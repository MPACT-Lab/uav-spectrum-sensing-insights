function [avgRSRPs, sampleIndices, sampleCounts, boxCentersX, boxCentersY, boxCentersZ,...
    boxCentersLon, boxCentersLat, boxCentersAlt] = make3DBoxes(latitudes, longitudes, altitudes, rsrps,...
    minLat, maxLat, minLon, maxLon, minAlt, maxAlt, boxLength, boxWidth, boxHeight)
    
    latDeltaToMeters = 6371000*2*pi/360; %111139 ;
    lonDeltaToMeters = 6371000*2*pi/360*cosd(35.727451);

    % Number of samples
    numSamples = length(rsrps);

    % Calculate the number of boxes along each dimension
    numLatBoxes = ceil((maxLat - minLat) * latDeltaToMeters / boxLength);
    numLonBoxes = ceil((maxLon - minLon) * lonDeltaToMeters / boxWidth);
    numAltBoxes = ceil((maxAlt - minAlt) / boxHeight);

    % Initialize the 3D matrices for average RSRPs and sample indices
    avgRSRPs = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    sampleIndices = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    sampleCounts = zeros(numLatBoxes, numLonBoxes, numAltBoxes);

    boxCentersLon = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    boxCentersLat = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    boxCentersAlt = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);

    boxCentersX = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    boxCentersY = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);
    boxCentersZ = NaN * ones(numLatBoxes, numLonBoxes, numAltBoxes);

    for i = 1:numLatBoxes
        for j = 1:numLonBoxes
            for k = 1:numAltBoxes
                boxCentersY(i,j,k) = (i - 1) * boxLength + boxLength / 2;
                boxCentersX(i,j,k) = (j - 1) * boxWidth + boxWidth / 2 ;
                boxCentersZ(i,j,k) = (k - 1) * boxHeight + boxHeight / 2;

                boxCentersLat(i,j,k) = minLat + (i - 1) * boxLength / latDeltaToMeters + boxLength / 2 / latDeltaToMeters;
                boxCentersLon(i,j,k) = minLon + (j - 1) * boxWidth / lonDeltaToMeters + boxWidth / 2 / lonDeltaToMeters;
                boxCentersAlt(i,j,k) = minAlt + (k - 1) * boxHeight + boxHeight / 2;
            end
        end
    end

    % Loop through all the samples
    for i = 1:numSamples
        % Get the latitude, longitude, and altitude of the current sample
        lat = latitudes(i);
        lon = longitudes(i);
        alt = altitudes(i);
        
        % Calculate the box indices for this sample
        latIndex = min(max(floor((lat - minLat) * latDeltaToMeters / boxLength) + 1, 1), numLatBoxes);
        lonIndex = min(max(floor((lon - minLon) * lonDeltaToMeters / boxWidth) + 1, 1), numLonBoxes);
        altIndex = min(max(floor((alt - minAlt) / boxHeight) + 1, 1), numAltBoxes);
        
%         % Calculate the center of the current box
%         boxCenterY = (latIndex - 1) * boxLength + boxLength / 2;
%         boxCenterX = (lonIndex - 1) * boxWidth + boxWidth / 2 ;
%         boxCenterZ = (altIndex - 1) * boxHeight + boxHeight / 2;

        % Calculate the center of the current box
        boxCenterLat = minLat + (latIndex - 1) * boxLength / latDeltaToMeters + boxLength / 2 / latDeltaToMeters;
        boxCenterLon = minLon + (lonIndex - 1) * boxWidth / lonDeltaToMeters + boxWidth / 2 / lonDeltaToMeters;
        boxCenterAlt = minAlt + (altIndex - 1) * boxHeight + boxHeight / 2;

        % Calculate the Euclidean distance from the sample to the center of the box
        distanceToCenter = sqrt((lat - boxCenterLat)^2 * latDeltaToMeters^2 + (lon - boxCenterLon)^2 * lonDeltaToMeters^2 + (alt - boxCenterAlt)^2);
        
        % If this box has not been filled yet, or this sample is closer to the center, update the average RSRP and sample index
        if isnan(avgRSRPs(latIndex, lonIndex, altIndex)) || distanceToCenter < sqrt((latitudes(sampleIndices(latIndex, lonIndex, altIndex)) - boxCenterLat)^2 * latDeltaToMeters^2 + (longitudes(sampleIndices(latIndex, lonIndex, altIndex)) - boxCenterLon)^2 * lonDeltaToMeters^2 + (altitudes(sampleIndices(latIndex, lonIndex, altIndex)) - boxCenterAlt)^2)
            % Store the index of the closest sample
            sampleIndices(latIndex, lonIndex, altIndex) = i;
        end

        % Update the average RSRP for this box
        currentRSRP = rsrps(i);
        % Update the average RSRP for this box
        if isnan(avgRSRPs(latIndex, lonIndex, altIndex))
            avgRSRPs(latIndex, lonIndex, altIndex) = currentRSRP;
        else
            avgRSRPs(latIndex, lonIndex, altIndex) = (currentRSRP + avgRSRPs(latIndex, lonIndex, altIndex) * sampleCounts(latIndex, lonIndex, altIndex)) / (1+sampleCounts(latIndex, lonIndex, altIndex));
        end
        sampleCounts(latIndex, lonIndex, altIndex) = sampleCounts(latIndex, lonIndex, altIndex) + 1;

    end
end

