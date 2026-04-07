clc;
clear;
close all;
% pre requisite threed_corr_gen_trans_ok_sk_ant.m (without normal transformation)
addpath('functions\')
addpath("../effect_ant_patt/patts\")

is_smooth = "_smooth"; % "_smooth" or "" % smooth for Kriging

all_heights = [30,50,70,90,110];
all_heights_train = [50,30,90,110,90]; % use data from another altitude for correlation coefficients
offset = [-78.698 35.727] ;
scaler = 111139 ;

R_earth = 6378.137 * 10^3 ; % [m]

lat_all_array = cell(length(all_heights),1); 
lon_all_array = lat_all_array; sha_all_array = lat_all_array; h_all_array = lat_all_array;

%% load shadow fading data
for i = 1:length(all_heights)
    load(sprintf('../data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_%d_meas_rad.mat',all_heights(i))) % 'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s', 'elev_s', height_s
    mat_name = "flipx_antenna_tx_rx_updated_sung_"+num2str(all_heights_train(i))+"m"+ is_smooth +".mat"; % for antenna pattern
    [sha_two_s, sha_free_s] = update_PL_sung_new_rad_pat(sha_two_s, sha_free_s, azim_s, elev_s, ori_s, mat_name);
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

num_known_samples = [50, 100, 200]; % per height

radius_all = [ 70, 200];

radius_colors = ["r", "g", "b", "c", "m", "y", "k"];
num_test_points = 25; % per height

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
error_all_2d_biased = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
num_meas_used_2d = zeros(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points);

variance_all_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;

lat_calc = zeros(num_repeat, length(all_heights), num_test_points);
lon_calc = zeros(num_repeat, length(all_heights), num_test_points);
h_calc = zeros(num_repeat, length(all_heights), num_test_points);


plot_fig = 0;

for repeat_id = 1:num_repeat
    fprintf('iter %d|\n',repeat_id);
    for height_id = 1:length(all_heights)
        load(sprintf("../corr_prof/profiles/correlation_table_3D_gauss_a2g_smooth_%dm.mat", all_heights_train(height_id))) % 'hor_dist_list','ver_dist_list', 'corr_list', 'mu_sha', "var_sha", 'd_int', 'h_int', 'correlation_f'
        load(sprintf("../corr_prof/profiles/normal_transformation_ant_%dm.mat", all_heights_train(height_id))) % trans_gauss trans_raw second_der_value

        corr_coef = [correlation_f.a, correlation_f.b, correlation_f.c, 0];
        added_term = var_sha*2;

        sha_all = sha_all_array{height_id};
        lat_all = lat_all_array{height_id};
        lon_all = lon_all_array{height_id};
        h_all = h_all_array{height_id};
        % target point set
        tot_points = length(sha_all);
        target_ids = randi([1 tot_points],1,num_test_points);
        if plot_inter_points>0
            all_indx = 1:tot_points;
            valid_indx = all_indx((h_all > analysis_tar_h - 10) & (h_all < analysis_tar_h + 10));
            target_ids_temp = randi([1, length(valid_indx)],1,num_test_points);
            target_ids = valid_indx(target_ids_temp);
        end
        lat_tar = lat_all(target_ids) ;
        lon_tar = lon_all(target_ids) ;
        h_tar = h_all(target_ids) ;
        sha_tar_raw = sha_all(target_ids) ;
        sha_tar = trans_gauss(sha_tar_raw);
        
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
            sha_meas_raw = sha_all(meas_ids) ;
            sha_meas = trans_gauss(sha_meas_raw);

            meas_correlation_all = my_meas_correlation_3D(lat_meas,lon_meas, h_meas, corr_coef, var_sha); %(num_meas, num_meas)
            target_cross_correlation_all = my_target_cross_correlation_3D(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef, var_sha); %(num_tar, num_meas)
        
            dist_ver_tar = abs(h_tar - h_meas'); %(num_tar, num_meas)
            dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas' * pi/180) .* cos( (lon_tar - lon_meas')  * pi/180 )) ;
            dist_target = sqrt(dist_ver_tar.^2 + dist_hor_tar.^2);
    
            if plot_inter_points>0
                %figure(repeat_id)
                FigH = figure('Position', get(0, 'Screensize'), 'visible','off');
                hold on
                same_height_indx = (h_meas > analysis_tar_h - 10) & (h_meas < analysis_tar_h + 10);
                scatter((lon_all(valid_indx)-offset(1))*scaler, (lat_all(valid_indx)-offset(2))*scaler,'CData',[0.7,0.7,0.7], 'DisplayName','trajectory')
                scatter((lon_meas(same_height_indx)-offset(1))*scaler,(lat_meas(same_height_indx)-offset(2))*scaler,[],sha_meas(same_height_indx),'filled','DisplayName','known points at same height')
                scatter((lon_meas(~same_height_indx)-offset(1))*scaler,(lat_meas(~same_height_indx)-offset(2))*scaler,[],sha_meas(~same_height_indx),'o', 'DisplayName','known points at another height')
                colorbar
                colormap jet
                ylim([-500, 400])
                xlim([-700,700])
                clim([-7 7])
                fig_text = sprintf('Height %d m, Known points %d, Radius %d',analysis_tar_h, num_known_samples(meas_i),radius_all(1));
                title("Analysis at "+fig_text)
                xlabel('X [m]')
                ylabel('Y [m]')
                grid on
                legend show
            end
         
            for radi_i = 1:length(radius_all)
                index_valid_radi_2d = (dist_hor_tar < radius_all(radi_i)); % & (dist_ver_tar<10) ;
    
                % 2D interpolation
                index_valid_radi = index_valid_radi_2d;
               [est_tar, est_variance, num_points, lagrange_mul] = my_oridinary_kriging_lagrange_output(meas_correlation_all,...
    target_cross_correlation_all, index_valid_radi, sha_meas, var_sha);
               err_tar = trans_raw(est_tar) + real(second_der_value*(est_variance/2 - lagrange_mul)) - sha_tar_raw;
               error_all_2d_biased(repeat_id, height_id, meas_i, radi_i, :) = trans_raw(est_tar) - sha_tar_raw;
               error_all_2d(repeat_id, height_id, meas_i, radi_i, :) = err_tar;
               num_meas_used_2d(repeat_id, height_id, meas_i, radi_i, :) = num_points;
               variance_all_2d(repeat_id, height_id, meas_i, radi_i, :) = est_variance;
            end
        end
    end
end

%% 
save results/res_kriging_ordinary_gaussian_sung_joon_50_250_a2g_3_5k.mat error_all_2d error_all_2d_biased num_meas_used_2d lat_calc lon_calc h_calc radius_all num_known_samples variance_all_2d all_heights
