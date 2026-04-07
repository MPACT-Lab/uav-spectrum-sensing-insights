clc;
clear;
close all;

addpath('functions\')
addpath('functions\Inpaint_nans\')

%load("filtererd_data_with_dist_ori_" + num2str(30) + "_meas_rad" + ".mat")
%% this is for rx antenna
altitudes = [110]; % change this and the next line as well {30, 50, 70, 90, 110}
altitude_string = num2str(altitudes(1))+"m"; 
% change _smoothing in ant_place_effect_sung.m as well ALLLEEERT!!!!!!!!!
is_smooth = "_smooth"; % "_smooth" or "" % smooth for Kriging

alpha = 1.0;

unknown_angle = 0;
mat_name = "patts/flipx_antenna_tx_rx_updated_sung_" + altitude_string + is_smooth + ".mat"; % for Kriging


prediction_err = zeros(1,0);

% with alpha = 1.0, i_iter=1 is enough
iter_end_at = 1;
for i_iter = 1:iter_end_at+1 %25
    % update ant patt for rx
    all_az = [];
    all_el = [];
    all_shadow = [];
    for ii=1:length(altitudes)
        altitude = altitudes(ii);
        load("../data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_" + num2str(altitude) + "_meas_rad" + ".mat")
        if i_iter>1
            [sha_two_s, sha_free_s] = update_PL_sung_new_rad_pat(sha_two_s, sha_free_s, azim_s, elev_s, ori_s, mat_name);
        end
        % -(phi_t + uav_yaw + 180)

        % ori_s is the direction following trajectory
        ori_s = ori_s - pi/2;
        ori_s = ori_s + (ori_s<-pi).*pi*2;

        azim_wrt_uav_ori = -azim_s - ori_s*pi/180;
        azim_wrt_uav_ori = azim_wrt_uav_ori - (azim_wrt_uav_ori>=0).*pi + (azim_wrt_uav_ori<0).*pi;
        all_az = [all_az; azim_wrt_uav_ori];
        all_el = [all_el; elev_s];
        all_shadow = [all_shadow; sha_free_s];
    end

    fig = figure(40);
    fig.Position = [20 + 500 * 2, 200, 500, 500];
    histogram(all_shadow);
    % Add labels and title
    xlabel('Shadow Fading [dB]');
    ylabel('Frequency');
    title('Histogram of residuals')
    xlim([-40, 40])

    prediction_err(1,i_iter) = sqrt(mean(all_shadow.^2)); % rmse

    %pause(1.0)

    if i_iter>iter_end_at
        break
    end
    % for i_iter=1, prev_patt and G_r_dB are same
    [Z, prev_patt, G_r_dB] = ant_place_effect_sung2(all_az,all_el,all_shadow*alpha,1,1,i_iter,101*(i_iter<2),102,1,1,mat_name);

end

fig = figure(70);
%fig.Position = [20 + 500 * 2, 200, 500, 500];
plot(1:length(prediction_err),prediction_err,'-o');
% Add labels and title
xlabel('Iteration Number');
ylabel('RMSE of Prediction (dB)');
% title('Histogram of residuals')
% xlim([-40, 40])

%% show rad pat before and after attachment (Azimuth)
% condition: iter_number must be 1, otherwise prev_patt will reflect
% intermidiate patt, Z also will be delta
% use [Z, prev_patt]
offset = 50;
figure;
placement_effect = zeros(61,120)*NaN;
placement_effect(1:31,:) = Z;
phi_grid = (180:-3:-177)*pi/180;
theta_grid = (90:-3:0)*pi/180;
theta_grid_full = (90:-3:-90)*pi/180;

% if altitude==30
%     elMax = 10;
% elseif altitude==70
%     elMax = 20;
% else
%     elMax = 30;
% end
elMax = 30;

theta_filter = (theta_grid_full>0*pi/180) & (theta_grid_full<elMax*pi/180); % 10 for 40m, 20 for 70m, 30 for 100m
prev_patt_subset = prev_patt(theta_filter,:);
orig_patt_subset = G_r_dB(theta_filter,:);
placement_effect_subset = placement_effect(theta_filter,:);

phi_filter = (phi_grid>-145*pi/180) & (phi_grid<50*pi/180);
% make this adaptive
phi_filter = ~isnan(nanmean(placement_effect_subset,1));

polarplot(phi_grid,max(mean(orig_patt_subset+offset,1),0),'--', 'DisplayName', 'Anechoic Chamber', 'LineWidth',2)
hold on;
polarplot(phi_grid(phi_filter),max(nanmean(prev_patt_subset(:,phi_filter)+offset+placement_effect_subset(:,phi_filter),1),0), 'DisplayName', 'UAV-to-ground environment', 'LineWidth',2)
hold on;
set(gcf,'color','w');
rlim([0 60])
%rticks([-35:10:5])
legend('Location', 'best');  % Adjust the label location as needed
% Increase figure text size
set(gca, 'FontSize', 14);  % Set the axis labels' font size

%% show rad pat before and after attachment (Elevation)
% condition: iter_number must be 1, otherwise prev_patt will reflect
% intermidiate patt, Z also will be delta
% use [Z, prev_patt]
offset = 50;
figure;
placement_effect = zeros(61,120)*NaN;
placement_effect(1:31,:) = Z;
phi_grid = (180:-3:-177)*pi/180;
theta_grid = (90:-3:0)*pi/180;
theta_grid_full = (90:-3:-90)*pi/180;

% if altitude=="40m"
%     elMax = 10;
% elseif altitude=="70m"
%     elMax = 20;
% else
%     elMax = 30;
% end

%azim_filter = phi_grid>?
prev_patt_subset = prev_patt;
orig_patt_subset = G_r_dB;
placement_effect_subset = placement_effect;

% make this adaptive
theta_filter = ~isnan(nanmean(placement_effect_subset,2));


polarplot(theta_grid_full,max(mean(orig_patt_subset+offset,2),0),'--', 'DisplayName', 'Anechoic Chamber', 'LineWidth',2)
hold on;
polarplot(theta_grid_full(theta_filter),max(nanmean(prev_patt_subset(theta_filter,:)+offset+placement_effect_subset(theta_filter,:),2),0), 'DisplayName', 'UAV-to-ground environment', 'LineWidth',2)
hold on;
set(gcf,'color','w');
rlim([0 70])
%rticks([-35:10:5])
legend('Location', 'best');  % Adjust the label location as needed
% Increase figure text size
set(gca, 'FontSize', 14);  % Set the axis labels' font size