clc;clear; close all;
addpath('.\GPR\functions\')
addpath('.\kriging\results\')
radius_colors = ["r", "b", "k", "m", "g", "c" ];
radius_plot_type = ["-s", "-<", "-x", "-d", "-o", "-v", "->"];
all_heihgts = [30,50,70,90,110];
all_radius = [70, 100, 200];

font_size = 14;

h_BS = 10 ; % [m]
lat_BS = 35.727451 ; % [degree]
lon_BS = -78.695974 ; % [degree]
R_earth = 6378.137 * 10^3 ; % [m]

load res_kriging_ordinary_sung_joon_50_250_3_5k.mat % baseline

fig_id1 = "latex_fig_803_";
fig_id2 = "latex_fig_802_";
num_known_samples1 = num_known_samples;
data1_all_heihgts = error_all_2d; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
data1_legend = 'TRPL, OK'; % change it
% for calc elev
data1_lat_all_h = lat_calc;
data1_lon_all_h = lon_calc;
data1_h_all_h = h_calc;
shaded_plot = 0;
%%
load res_kriging_ordinary_sung_joon_50_250_3_a2g_5k.mat

% for calc elev
data2_lat_all_h = lat_calc;
data2_lon_all_h = lon_calc;
data2_h_all_h = h_calc;

data2_all_heihgts = error_all_2d;
data2_legend = 'Calibrated, OK';

y_label = 'Median of RMSE [dB]'; % change it
y_label2 = 'Valid estimation ratio'; %

elev_grids = (5:10:70)'*pi/180;
% futher smoothed on this code (idx_theta1==ii1)|(idx_theta1==ii1+1)|(idx_theta1==ii1-1)

data1_lat = data1_lat_all_h;
data1_lon = data1_lon_all_h;
data1_h = data1_h_all_h;
data2_lat = data2_lat_all_h;
data2_lon = data2_lon_all_h;
data2_h = data2_h_all_h;

% calc elevation
dist_2D = R_earth .* acos(sin(lat_BS * pi/180) .* sin(data1_lat * pi/180) + cos(lat_BS * pi/180) .* cos(data1_lat * pi/180) .* cos( (data1_lon - lon_BS)  * pi/180 )) ;
elev1 = atan( ( data1_h - h_BS )./dist_2D) ;
dist_2D = R_earth .* acos(sin(lat_BS * pi/180) .* sin(data2_lat * pi/180) + cos(lat_BS * pi/180) .* cos(data2_lat * pi/180) .* cos( (data2_lon - lon_BS)  * pi/180 )) ;
elev2 = atan( ( data2_h - h_BS )./dist_2D) ;
elev1_size = size(elev1);

elev1_lin = elev1(:);
[~,idx_theta] = min(abs(elev1_lin - elev_grids')') ;
figure
histogram(idx_theta)
idx_theta1 = reshape(idx_theta, elev1_size);

elev2_lin = elev2(:);
[~,idx_theta] = min(abs(elev2_lin - elev_grids')') ;
figure
histogram(idx_theta)
idx_theta2 = reshape(idx_theta, elev1_size);

elev_fspl = [];
fspl_sha_all = [];
trpl_sha_all = [];
% results for FSPL and TRPL
for i = 1:length(all_heihgts)
    load(sprintf('./data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_%d_meas_rad2.mat',all_heihgts(i))) % 'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s', 'elev_s', height_s
    elev_fspl = [elev_fspl; elev_s];
    fspl_sha_all = [fspl_sha_all; sha_free_s];
    trpl_sha_all = [trpl_sha_all; sha_two_s];
end
RMSE_fspl = [];
RMSE_trpl = [];
variable_smooth_range = [5,5,10,20,30,30,30];
for ii1 = 1:length(elev_grids)
    elev_filter = abs(elev_fspl - elev_grids(ii1)) < variable_smooth_range(ii1)*pi/180;
    RMSE_fspl = [RMSE_fspl, sqrt(var(fspl_sha_all(elev_filter)))];
    RMSE_trpl = [RMSE_trpl, sqrt(var(trpl_sha_all(elev_filter)))];
end
RMSE_fspl
RMSE_trpl
%% 
for height_i = 1:length(all_radius)
    current_heihgt = all_radius(height_i);
    data1 = squeeze(data1_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,radius,targets)
    data2 = squeeze(data2_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,num_fixed_points,targets)

    calc_rmse_points  = 20;


    res_rmse_elev_meas1 = zeros(length(elev_grids), size(data1, 3));
    res_rmse_elev_meas2 = zeros(length(elev_grids), size(data1, 3));
    

    for ii1 = 1:length(elev_grids)

        
        for ii2 = 1:size(data1, 3) % meas 50 100 200
            data1_sub = squeeze(data1(:,:,ii2,:)); % now matches shape with idx_theta1;
            data2_sub = squeeze(data2(:,:,ii2,:)); % now matches shape with idx_theta1;
            data1_sub_elev = data1_sub((idx_theta1==ii1)|(idx_theta1==ii1+1)|(idx_theta1==ii1-1));
            data2_sub_elev = data2_sub((idx_theta2==ii1)|(idx_theta2==ii1+1)|(idx_theta2==ii1-1));
            data1_sub_elev = data1_sub_elev(~isnan(data1_sub_elev));
            data2_sub_elev = data2_sub_elev(~isnan(data2_sub_elev));

            medain_points1 = floor(length(data1_sub_elev)/calc_rmse_points);
            medain_points2 = floor(length(data2_sub_elev)/calc_rmse_points);
            data_mad1 = zeros(medain_points1, 1);
            data_mad2 = zeros(medain_points2, 1);
            for ii3 = 1:medain_points1
                data_mad1(ii3,1) = sqrt(mean(data1_sub_elev((ii3-1)*calc_rmse_points+1:ii3*calc_rmse_points).^2));
            end
            for ii3 = 1:medain_points2
                data_mad2(ii3,1) = sqrt(mean(data2_sub_elev((ii3-1)*calc_rmse_points+1:ii3*calc_rmse_points).^2));
            end
    
            res_rmse_elev_meas1(ii1, ii2) = median(data_mad1);
            res_rmse_elev_meas2(ii1, ii2) = median(data_mad2);
        end

    end


    %FigH = figure('Position', get(0, 'Screensize'), 'visible','off');
    FigH = figure('Position', [100,100,680,610], 'visible','on');
    hax1=axes;
    hold on
    
    new_iter = 1;
    for radi_i = [1,2,3] %[1,3,7] % 50, 100, 200 measurements
        plot(elev_grids*180/pi, res_rmse_elev_meas1(:,radi_i),radius_colors(1)+'-'+radius_plot_type(new_iter),'DisplayName', sprintf('%s (M = %d)',data1_legend,num_known_samples(radi_i)),'LineWidth',2)
        plot(elev_grids*180/pi, res_rmse_elev_meas2(:,radi_i),radius_colors(2)+radius_plot_type(new_iter),'DisplayName', sprintf('%s (M = %d)',data2_legend,num_known_samples(radi_i)),'LineWidth',2)
        new_iter = new_iter + 1;
    end
    plot(elev_grids*180/pi, RMSE_trpl, '--','DisplayName', 'Shadow Fading $\sigma$ (TRPL)', 'LineWidth',2,'color',[0.0, 0.0, 0.0])
    plot(elev_grids*180/pi, RMSE_fspl, ':','DisplayName', 'Shadow Fading $\sigma$ (FSPL)', 'LineWidth',2,'color',[0.0, 0.0, 0.0])

    xlabel('Elevation Angle [degree]')
    ylabel(y_label)
    %ylim([5.4 7.4])
    %xticks(num_known_samples1)
    grid on
    axis tight
    legend show
    legend('Interpreter', 'latex');
    fig_text = fig_id1 + y_label + " (R = " + num2str(current_heihgt) + " m)";
    fontsize(FigH,font_size,"points")
    box on
    ylim([5, 7.75])

end

