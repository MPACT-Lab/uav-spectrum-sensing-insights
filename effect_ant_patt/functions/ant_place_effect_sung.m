function [Z, prev_patt, G_r_dB] = ant_place_effect_sung(all_az,all_el,error_all,save_patt,is_rx,n_iter,fig_prv, fig_new, fig_variance, verbose, save_name)
%ANT_PLACE_EFFECT Summary of this function goes here
%   Detailed explanation goes here
all_shadow = error_all;
%% plot
if verbose
    figure
    scatter(all_az, all_el, 50, all_shadow, 'filled'); % 'filled' option fills the markers
    %caxis([-50,20])
    
    % Adding labels and title
    xlabel('Azimuth [deg]');
    ylabel('Elevation [deg]');
    title('Shadow Fading (dB)');
    colormap(jet);
    % Adding a colorbar
    colorbar;
end

%% interpolation
phi_grid = (180:-3:-177)*pi/180;
theta_grid = (90:-3:0)*pi/180;
theta_grid_full = (90:-3:-90)*pi/180;
% phi_grid = (-177:3:180)*pi/180;
% theta_grid = (0:3:90)*pi/180;
% theta_grid_full = (-90:3:90)*pi/180;
d_grid = 3*pi/180;

%addpath("Inpaint_nans\")
[phi_mesh,theta_mesh] = meshgrid(phi_grid,theta_grid);

Z = zeros(size(phi_mesh));
Z_var = zeros(size(phi_mesh));
for ii=1:length(phi_grid)
    for jj=1:length(theta_grid)
        phi_this = phi_grid(ii);
        theta_this = theta_grid(jj);
        shadow_this = all_shadow((abs(all_az-phi_this)<d_grid) & (abs(all_el-theta_this)<d_grid));
        Z(jj,ii) = mean(shadow_this); % mean([]) will be nan
        Z_var(jj,ii) = var(shadow_this);
    end
end

% Z_filled = inpaint_nans(Z,1);
Z_filled = Z;
Z_filled(isnan(Z_filled)) = 0;
% Define Gaussian kernel parameters
sigma = 3; % Standard deviation of the Gaussian
filterSize = 15; % Size of the filter (e.g., 5x5)
gaussianKernel = fspecial('gaussian', filterSize, sigma);
Z_filled = conv2(Z_filled, gaussianKernel, 'same'); % 'same' keeps output size
Z_filled = conv2(Z_filled, gaussianKernel, 'same'); % 'same' keeps output size

if verbose
    figure
    mesh(phi_grid*180/pi,theta_grid*180/pi,Z,'FaceColor','flat')
    % Adding labels and title
    xlabel('Azimuth [deg]');
    ylabel('Elevation [deg]');
    title('Shadow Fading (dB)');
    colormap(jet);
    % Adding a colorbar
    colorbar
    %caxis([-50,20])
    view(2)
    figure
    mesh(phi_grid*180/pi,theta_grid*180/pi,Z_filled,'FaceColor','flat')
    % Adding labels and title
    xlabel('Azimuth [deg]');
    ylabel('Elevation [deg]');
    title('Shadow Fading (dB)');
    colormap(jet);
    % Adding a colorbar
    colorbar
    %caxis([-50,20])
    view(2)
    
    % figure
    % histogram(Z(:))
    % 
    % figure
    % histogram(Z_filled(:))
end

Z_filled(Z_filled<min(Z(:))) = min(Z(:));
Z_filled(Z_filled>max(Z(:))) = max(Z(:));

%% 
load rad_pat_orig.mat
if n_iter > 1
    eval("load " + save_name)
end

placement_effect = zeros(61,120);
placement_effect(1:31,:) = Z_filled;
%placement_effect(31:end,:) = Z_filled;

if is_rx
    if n_iter == 1
        rx_patt_with_placement_effect = G_r_dB + placement_effect; % in baseline.m
        prev_patt = G_r_dB;
        tx_patt_with_placement_effect = G_t_dB; % initial (default)
    else
        rx_patt_with_placement_effect = rx_patt_with_placement_effect + placement_effect; % in baseline_ant_placement.m
        prev_patt = rx_patt_with_placement_effect;
    end
    new_patt = rx_patt_with_placement_effect;
else
    if n_iter == 1
        tx_patt_with_placement_effect = SA3300; % in baseline.m
        prev_patt = SA3300;
        rx_patt_with_placement_effect = G_r_dB; % initial
    else
        tx_patt_with_placement_effect = tx_patt_with_placement_effect + placement_effect; % in baseline_ant_placement.m
        prev_patt = tx_patt_with_placement_effect;
    end
    new_patt = tx_patt_with_placement_effect;
end

if(fig_prv)
    fig = figure(fig_prv);
    fig.Position = [20 + (fig_prv>102)*500, 200, 500, 500];
    surf(phi_grid*180/pi,theta_grid_full*180/pi,prev_patt);
    colormap(jet);
    colorbar;
    %mesh(phi_grid,theta_grid_full,theta_grid_full,'FaceColor','flat')
    %colorbar
    %caxis([min(prev_patt(:)) max(prev_patt(:))])
    caxis([-15 15])
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Antenna gain [dBi]';
    set(colorTitleHandle ,'String',titleString);
    set(gcf,'color','w')
    xlabel('Azimuth angle [degree]'); ylabel('Elevation angle [degree]'); zlabel('Antenna gain [dB]');
    set(gcf,'Color','w');
    axis([-177 180 -90 90 min(prev_patt(:)) max(prev_patt(:)) ]);
    view([0,0,1])
end



if(fig_variance)
    fig = figure(fig_variance);
    summarize_variance = zeros(61,120);
    summarize_variance(1:31,:) = Z_var;
    %summarize_variance(31:end,:) = Z_var;
    surf(phi_grid*180/pi,theta_grid_full*180/pi,summarize_variance);
    colormap(jet);
    colorbar;
    %mesh(phi_grid,theta_grid_full,theta_grid_full,'FaceColor','flat')
    %colorbar
    caxis([min(summarize_variance(:)) 70])
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Variance [dBi]';
    set(colorTitleHandle ,'String',titleString);
    set(gcf,'color','w')
    xlabel('Azimuth angle [degree]'); ylabel('Elevation angle [degree]'); zlabel('Variance [dB]');
    set(gcf,'Color','w');
    axis([-177 180 -90 90 min(summarize_variance(:)) max(summarize_variance(:)) ]);
    view([0,0,1])
end

if(fig_new)
    fig = figure(fig_new);
    fig.Position = [20 + (fig_new>102)*500, 200, 500, 500];
    %mesh(phi_grid,theta_grid_full,rx_patt_with_placement_effect,'FaceColor','flat')
    surf(phi_grid*180/pi,theta_grid_full*180/pi,new_patt);
    colormap(jet);
    colorbar;
    %mesh(phi_grid,theta_grid_full,theta_grid_full,'FaceColor','flat')
    %colorbar
    %caxis([min(new_patt(:)) max(new_patt(:))])
    caxis([-15 15])
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Antenna gain [dBi]';
    set(colorTitleHandle ,'String',titleString);
    %set(gcf,'color','w')
    xlabel('Azimuth angle [degree]'); ylabel('Elevation angle [degree]'); zlabel('Antenna gain [dB]');
    set(gcf,'Color','w');
    axis([-177 180 -90 90 min(new_patt(:)) max(new_patt(:)) ]);
    view([0,0,100])
end

if save_patt
    eval("save " + save_name + " phi2 theta2 G_t_dB G_r_dB rx_patt_with_placement_effect tx_patt_with_placement_effect")
end
end

