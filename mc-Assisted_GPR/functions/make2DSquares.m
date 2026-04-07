function [avgRSRPs, sampleIndices, sampleCounts, boxCentersX, boxCentersY,...
    boxCentersLon, boxCentersLat] = make2DSquares(latitudes, longitudes, rsrps,...
    minLat, maxLat, minLon, maxLon, boxLength, boxWidth)
    
    latDeltaToMeters = 6371000*2*pi/360; %111139 ;
    lonDeltaToMeters = 6371000*2*pi/360*cosd(35.727451);

    % Number of samples
    numSamples = length(rsrps);

    % Calculate the number of boxes along each dimension
    numLatBoxes = ceil((maxLat - minLat) * latDeltaToMeters / boxLength);
    numLonBoxes = ceil((maxLon - minLon) * lonDeltaToMeters / boxWidth);

    % Initialize the 3D matrices for average RSRPs and sample indices
    avgRSRPs = NaN * ones(numLatBoxes, numLonBoxes);
    sampleIndices = NaN * ones(numLatBoxes, numLonBoxes);
    sampleCounts = zeros(numLatBoxes, numLonBoxes);

    boxCentersLon = NaN * ones(numLatBoxes, numLonBoxes);
    boxCentersLat = NaN * ones(numLatBoxes, numLonBoxes);

    boxCentersX = NaN * ones(numLatBoxes, numLonBoxes);
    boxCentersY = NaN * ones(numLatBoxes, numLonBoxes);

    for i = 1:numLatBoxes
        for j = 1:numLonBoxes
                boxCentersY(i,j) = (i - 1) * boxLength + boxLength / 2;
                boxCentersX(i,j) = (j - 1) * boxWidth + boxWidth / 2 ;

                boxCentersLat(i,j) = minLat + (i - 1) * boxLength / latDeltaToMeters + boxLength / 2 / latDeltaToMeters;
                boxCentersLon(i,j) = minLon + (j - 1) * boxWidth / lonDeltaToMeters + boxWidth / 2 / lonDeltaToMeters;
        end
    end

    % Loop through all the samples
    for i = 1:numSamples
        % Get the latitude, longitude, and altitude of the current sample
        lat = latitudes(i);
        lon = longitudes(i);
        
        % Calculate the box indices for this sample
        latIndex = min(max(floor((lat - minLat) * latDeltaToMeters / boxLength) + 1, 1), numLatBoxes);
        lonIndex = min(max(floor((lon - minLon) * lonDeltaToMeters / boxWidth) + 1, 1), numLonBoxes);
        
%         % Calculate the center of the current box
%         boxCenterY = (latIndex - 1) * boxLength + boxLength / 2;
%         boxCenterX = (lonIndex - 1) * boxWidth + boxWidth / 2 ;
%         boxCenterZ = (altIndex - 1) * boxHeight + boxHeight / 2;

        % Calculate the center of the current box
        boxCenterLat = minLat + (latIndex - 1) * boxLength / latDeltaToMeters + boxLength / 2 / latDeltaToMeters;
        boxCenterLon = minLon + (lonIndex - 1) * boxWidth / lonDeltaToMeters + boxWidth / 2 / lonDeltaToMeters;


        % Calculate the Euclidean distance from the sample to the center of the box
        distanceToCenter = sqrt((lat - boxCenterLat)^2 * latDeltaToMeters^2 + (lon - boxCenterLon)^2 * lonDeltaToMeters^2);
        
        % If this box has not been filled yet, or this sample is closer to the center, update the average RSRP and sample index
        if isnan(avgRSRPs(latIndex, lonIndex)) || distanceToCenter < sqrt((latitudes(sampleIndices(latIndex, lonIndex)) - boxCenterLat)^2 * latDeltaToMeters^2 + (longitudes(sampleIndices(latIndex, lonIndex)) - boxCenterLon)^2 * lonDeltaToMeters^2)
            % Store the index of the closest sample
            sampleIndices(latIndex, lonIndex) = i;
        end

        % Update the average RSRP for this box
        currentRSRP = rsrps(i);
        % Update the average RSRP for this box
        if isnan(avgRSRPs(latIndex, lonIndex))
            avgRSRPs(latIndex, lonIndex) = currentRSRP;
        else
            avgRSRPs(latIndex, lonIndex) = (currentRSRP + avgRSRPs(latIndex, lonIndex) * sampleCounts(latIndex, lonIndex)) / (1+sampleCounts(latIndex, lonIndex));
        end
        sampleCounts(latIndex, lonIndex) = sampleCounts(latIndex, lonIndex) + 1;

    end
end

