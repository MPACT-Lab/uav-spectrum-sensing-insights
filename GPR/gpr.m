clc;
clear;
close all;
addpath('functions\')

% pre-requisite threed_corr_gen_GPR.m

all_heights = [70];
all_heights_train = [90];
offset = [-78.698 35.727] ;
scaler = 111139 ;
%h_BS = 10 ; % [m]
%lat_BS = 35.727451 ; % [degree]
%lon_BS = -78.695974 ; % [degree]
R_earth = 6378.137 * 10^3 ; % [m]
%c = 3.0e8;
%f0 = 3.5e9;
%lamda = c/f0;

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

%% start
close all
num_repeat = 5000;
num_known_samples = [50, 150, 250, 350, 450]; % per height
num_known_samples = [50, 75, 100, 125, 150, 175, 200, 225, 250]; % per height
num_known_samples = [50, 75, 100, 125, 150, 175, 200, 225, 250]; % per height
radius_all = [ 70, 200];

radius_colors = ["r", "g", "b", "c", "m", "y", "k"];
num_test_points = 25; % per height

%% 

%error_all_3d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
%num_meas_used_3d = zeros(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points);
error_all_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
num_meas_used_2d = zeros(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points);

lat_calc = zeros(num_repeat, length(all_heights), num_test_points);
lon_calc = zeros(num_repeat, length(all_heights), num_test_points);
h_calc = zeros(num_repeat, length(all_heights), num_test_points);


plot_fig = 0;

for height_id = 1:length(all_heights)
    load(sprintf("./gpr_coeffs/gpr_coeff_%dm_rad.mat", all_heights_train(height_id))) 
    corr_coef = [gpr_values.a, gpr_values.b, gpr_values.c, 0];
    var_sha = gpr_values.v; mu_sha = 0;
    var_gp = gpr_values.sig;
    
    added_term = var_sha*2;
    for repeat_id = 1:num_repeat
        if mod(repeat_id,100)==0
            fprintf('iter %d|\n',repeat_id);
        end
        sha_all = sha_all_array{height_id};
        lat_all = lat_all_array{height_id};
        lon_all = lon_all_array{height_id};
        h_all = h_all_array{height_id};
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
            
            meas_correlation_all = my_meas_correlation_2D(lat_meas,lon_meas, h_meas, corr_coef, var_sha); %(num_meas, num_meas)
            target_cross_correlation_all = my_target_cross_correlation_2D(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef, var_sha); %(num_tar, num_meas)
        
            %dist_ver_tar = abs(h_tar - h_meas'); %(num_tar, num_meas)
            dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas' * pi/180) .* cos( (lon_tar - lon_meas')  * pi/180 )) ;
            %dist_target = sqrt(dist_ver_tar.^2 + dist_hor_tar.^2);
         
            for radi_i = 1:length(radius_all)
                % 3D interpolation
                %index_valid_radi_3d = (dist_target < radius_all(radi_i)) & (dist_ver_tar<50);
                index_valid_radi_2d = (dist_hor_tar < radius_all(radi_i)); % & (dist_ver_tar<10) ;
    
                % end of 3D interpolation
                % 2D interpolation
                index_valid_radi = index_valid_radi_2d;

               [est_tar, num_points] = my_gpr(meas_correlation_all,...
    target_cross_correlation_all, index_valid_radi, sha_meas, var_sha, var_gp);
               error_all_2d(repeat_id, height_id, meas_i, radi_i, :) = est_tar - sha_tar;
               num_meas_used_2d(repeat_id, height_id, meas_i, radi_i, :) = num_points;
            end
        end
    end
end

%% 
save results/res_gpr_lte_70m.mat error_all_2d num_meas_used_2d lat_calc lon_calc h_calc radius_all num_known_samples all_heights
