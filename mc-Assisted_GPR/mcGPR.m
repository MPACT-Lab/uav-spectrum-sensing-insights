clc;
clear;
close all;
% pre requisite threed_corr_gen_meas.m (without normal transformation)
addpath('functions\')
addpath('MatrixCompletion\')

all_heights = [70];
all_heights_train = [90];
offset = [-78.698 35.727] ;
scaler = 111139 ;

R_earth = 6378.137 * 10^3 ; % [m]

lat_all_array = cell(length(all_heights),1); 
lon_all_array = lat_all_array; sha_all_array = lat_all_array; h_all_array = lat_all_array;

%% load shadow fading data
for i = 1:length(all_heights)
    load(sprintf('../data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_%d_meas_rad.mat',all_heights(i))) % 'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s', 'elev_s', height_s
    lat_all_array{i} = lat_s;
    lon_all_array{i} = lon_s;
    sha_all_array{i} = sha_two_s;
    h_all_array{i} = height_s;
    clear lat_s lon_s height_s sha_two_s RSRP_s sha_free_s RSRP_PL_two_s RSRP_PL_free_s elev_s x_f y_f
end

latDeltaToMeters = 6371000*2*pi/360; %111139 ;
lonDeltaToMeters = 6371000*2*pi/360*cosd(35.727451);
minLat = 35.7232; maxLat = 35.7303; minLon = -78.7002; maxLon = -78.6961;
%% start
close all
num_repeat = 5000;
num_known_samples = [50, 75, 100, 125, 150, 175, 200, 225, 250];
radius_all = [1];
    
radius_colors = ["r", "g", "b", "c", "m", "y", "k"];
num_test_points = 25; % per height
threed = 0; % exp(-y) improves performance, although it makes variance bigger (no one close to 0)

plot_inter_points = 0;

if plot_inter_points>0
    num_known_samples = [100];
    radius_all = [100];
    num_test_points = 7;
    num_repeat = 5;
    analysis_tar_h = 70;
end

%% 
error_all_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
num_meas_used_2d = zeros(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points);

error_all_2d_gpr = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
error_all_2d_gpr_grid = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;

lat_calc = zeros(num_repeat, length(all_heights), num_test_points);
lon_calc = zeros(num_repeat, length(all_heights), num_test_points);
h_calc = zeros(num_repeat, length(all_heights), num_test_points);

plot_fig = 0;

rng(58)

for height_id = 1:length(all_heights)
    load(sprintf("./gpr_coeffs/gpr_coeff_%dm_rad.mat", all_heights_train(height_id))) 
    corr_coef = [gpr_values.a, gpr_values.b, gpr_values.c, 0];
    var_sha = gpr_values.v; mu_sha = 0;
    var_gp = gpr_values.sig;
    for repeat_id = 1:num_repeat
        if mod(repeat_id,10)==0
            fprintf('iter %d|\n',repeat_id);
        end

        sha_all = sha_all_array{height_id};
        lat_all = lat_all_array{height_id};
        lon_all = lon_all_array{height_id};
        h_all = h_all_array{height_id};

        boxLength = 5; boxWidth =5;

        % target point set
        tot_points = length(sha_all);
        target_ids = randi([1 tot_points],1,num_test_points);

        lat_tar = lat_all(target_ids) ;
        lon_tar = lon_all(target_ids) ;
        h_tar = h_all(target_ids) ;
        sha_tar = sha_all(target_ids) ;
        
        lat_calc(repeat_id, height_id, :) = lat_tar;
        lon_calc(repeat_id, height_id, :) = lon_tar;
        h_calc(repeat_id, height_id, :) = h_tar;
    
        meas_ids_probable_set = setdiff(1:tot_points, target_ids);
    
        for meas_i = 1:length(num_known_samples)
            % meas point set
            num_meas = num_known_samples(meas_i);
            meas_ids_temp = unique(randi([1, tot_points-num_test_points],1,ceil(num_meas*2)));
            randomize_idx = randperm(length(meas_ids_temp));
            meas_ids = meas_ids_probable_set(meas_ids_temp(randomize_idx(1:num_meas)));
            
            lat_meas = lat_all(meas_ids) ;
            lon_meas = lon_all(meas_ids) ;
            h_meas = h_all(meas_ids) ;
            sha_meas = sha_all(meas_ids) ;

            [avgSFs, sampleIndices, sampleCounts, centerX, centerY, centerLon, centerLat] = make2DSquares(lat_meas, lon_meas, sha_meas, minLat, maxLat, minLon, maxLon, boxLength, boxWidth);
            use_raw_meas = 1;

            meas_correlation_all = my_meas_correlation_2D(lat_meas,lon_meas, h_meas, corr_coef, var_sha); %(num_meas, num_meas)
            target_cross_correlation_all = my_target_cross_correlation_2D(lat_meas,lon_meas, h_meas, centerLat(:), centerLon(:), zeros(size(centerLon(:))), corr_coef, var_sha); %(num_tar, num_meas)

            dist_hor_tar = R_earth .* acos(sin(centerLat(:) * pi/180) .* sin(lat_meas' * pi/180) + cos(centerLat(:) * pi/180) .* cos(lat_meas' * pi/180) .* cos( (centerLon(:) - lon_meas')  * pi/180 )) ;

            index_valid_radi_2d = (dist_hor_tar < 70); % & (dist_ver_tar<10) ;

            % predict on the centers using GPR
            [centerSF, centerVar, ~] = my_gpr(meas_correlation_all,...
    target_cross_correlation_all, index_valid_radi_2d, sha_meas, var_sha, var_gp);
            avgSFs = reshape(centerSF, size(avgSFs,1), size(avgSFs,2));
            centerVars = reshape(centerVar, size(avgSFs,1), size(avgSFs,2));
            centerVars = 20+max(centerVars(:))*(centerVars-min(centerVars(:)))./(max(centerVars(:))-min(centerVars(:)));

            mask = ~isnan(avgSFs); % this should be all true, as Kriging prediction is everywhere
            [~, nearest_idx] = bwdist(mask);
            avgSFs_completed = avgSFs(nearest_idx);

            mask_all_true = ~isnan(avgSFs_completed);
            mask_all_false = false(size(mask));

            n_iter = 1200;
            matrix_change_tolerance = 1; % enforces to match with given values, decrease -> more close to original value (variance also matters)
            norm_stability_threshold = 20; % how much norm you want to reduce? 10-> almost 0 norm
            [compMat, status] = MatrixCompletionIntervalMatrix(avgSFs_completed, mask,...
                n_iter,'nuclear', norm_stability_threshold, sqrt(centerVars)*matrix_change_tolerance, 0);

            %% 
            prominent_trace = avgSFs_completed - compMat;
            
            prominent_trace_clear = prominent_trace;
            prominent_trace_clear(abs(prominent_trace)<1) = NaN;

            prominent_trace_clear(isnan(prominent_trace_clear)) = 0;
            smooth_data = avgSFs_completed - prominent_trace_clear;
            sigma = 5; % Standard deviation of the Gaussian
            filterSize = 7; % Size of the filter (e.g., 5x5)
            prominent_trace_clear_pos = prominent_trace_clear .* (prominent_trace_clear >=0);
            prominent_trace_clear_neg = prominent_trace_clear .* (prominent_trace_clear <0);
            prominent_trace_clear_pos = gaussianMaxDilation(prominent_trace_clear_pos, sigma, filterSize); % 'same' keeps output size
            prominent_trace_clear_neg = -gaussianMaxDilation(-prominent_trace_clear_neg, sigma, filterSize); % 'same' keeps output size
            prominent_trace_clear_res = prominent_trace_clear_pos + prominent_trace_clear_neg;


            sigma = 5; % Standard deviation of the Gaussian
            filterSize = 11; % Size of the filter (e.g., 5x5)
            gaussianKernel = fspecial('gaussian', filterSize, sigma);
            prominent_trace_clear_res = conv2(prominent_trace_clear_res, gaussianKernel, 'same'); % 'same' keeps output size
          
            
            compMat = smooth_data + prominent_trace_clear_res;



            %   linear', 'nearest', 'cubic', 'makima', or 'spline'
            sha_tar_est = interp2(centerLon,centerLat,compMat,lon_tar,lat_tar,'spline');
            error_all_2d(repeat_id,height_id,meas_i,1,:) = sha_tar_est - sha_tar;

            sha_tar_est_center_gpr = interp2(centerLon,centerLat,avgSFs_completed,lon_tar,lat_tar,'spline');
            error_all_2d_gpr_grid(repeat_id,height_id,meas_i,1,:) = sha_tar_est_center_gpr - sha_tar;

            % predict on the test set using GPR
            target_cross_correlation_all = my_target_cross_correlation_2D(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef, var_sha); %(num_tar, num_meas)

            [est_tar, num_points] = my_gpr2(meas_correlation_all,...
    target_cross_correlation_all, true(size(target_cross_correlation_all)), sha_meas, var_sha, var_gp);
            error_all_2d_gpr(repeat_id, height_id, meas_i, 1, :) = est_tar - sha_tar;

        end
    end
end

%% 
save results\res_mc4_lte_70m.mat error_all_2d error_all_2d_gpr error_all_2d_gpr_grid num_meas_used_2d lat_calc lon_calc h_calc radius_all num_known_samples all_heights threed