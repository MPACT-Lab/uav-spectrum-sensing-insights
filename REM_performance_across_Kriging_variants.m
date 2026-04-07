clc;clear; close all;
addpath('myFunctions\')
radius_colors = ["r", "b", "k", "m", "g", "c" ];
radius_colors = ["r", "r", "b", "b", "g", "c" ];
%radius_colors = ["k", "k", "k", "k", "k", "k" ];
radius_plot_type = ["-s", "-<", "-x", "-d", "-o", "-v", "->"];
all_heihgts = [30,50,70,90,110];
all_radius = [70,200];
font_size = 14;

load res_kriging_ordinary_sung_joon_50_250_3_a2g_5k.mat % baseline
%make 70,100,200 to 70,200
error_all_2d = error_all_2d(:,:,:,[1,3],:);
num_meas_used_2d = num_meas_used_2d(:,:,:,[1,3],:);
variance_all_2d = variance_all_2d(:,:,:,[1,3],:);
radius_all = [70,200];

fig_id1 = "latex_fig_803_";
fig_id2 = "latex_fig_802_";
num_known_samples1 = num_known_samples;
data1_all_heihgts = error_all_2d; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
data1_legend = 'Calibrated, OK'; % change it
shaded_plot = 0;
%%
% not ready yet
load res_kriging_ordinary_gaussian_sung_joon_50_250_3_a2g_5k.mat
error_all_2d = error_all_2d(:,:,:,[1,3],:);
num_meas_used_2d = num_meas_used_2d(:,:,:,[1,3],:);
variance_all_2d = variance_all_2d(:,:,:,[1,3],:);
radius_all = [70,200];

data2_all_heihgts = error_all_2d;
data2_legend = 'Calibrated, OK, Gaussian';

y_label = 'Median of RMSE [dB]'; % change it
y_label2 = 'Valid estimation ratio'; %

load res_kriging_simple_sung_joon_50_250_2_a2g_5k.mat % baseline
data3_all_heihgts = error_all_2d; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
data3_legend = 'Calibrated, SK'; % change it

load res_kriging_simple_gaussian_sung_joon_50_250_2_a2g_5k.mat
data4_all_heihgts = error_all_2d;
data4_legend = 'Calibrated, SK, Gaussian';

for height_i = 1:length(all_radius)
    current_heihgt = all_radius(height_i);
    data1 = squeeze(data1_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,radius,targets)
    data2 = squeeze(data2_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,num_fixed_points,targets)
    data3 = squeeze(data3_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,radius,targets)
    data4 = squeeze(data4_all_heihgts(:,:,:,height_i,:)); % (iter,height,num_meas,num_fixed_points,targets)

    % previously we had (num_meas, radius)
    % now we will have (height, num_meas)

    data1 = data1(:,:,:,1:20);
    data2 = data2(:,:,:,1:20);
    data3 = data3(:,:,:,1:20);
    data4 = data4(:,:,:,1:20);


    [rmse_meas_radius_data1, rmse_meas_radius_data2, valid_est1, valid_est2] = calc_rmse_median_sung_joon(data1,data2);
    [rmse_meas_radius_data3, rmse_meas_radius_data4, valid_est3, valid_est4] = calc_rmse_median_sung_joon(data3,data4);
    %FigH = figure('Position', get(0, 'Screensize'), 'visible','off');
    FigH = figure('Position', [100,100,680,610], 'visible','on');
    hax1=axes;
    hold on
    
    new_iter = 1;
    for radi_i = [1,3] %[1,3,7] % 50, 100, 200 measurements
        if radi_i>1
        plotLineStyle = '-';
        else
            plotLineStyle = '';
        end
        plot(all_heihgts, rmse_meas_radius_data1(:,radi_i),radius_colors(1)+plotLineStyle+radius_plot_type(1),'DisplayName', sprintf('%s (M = %d)',data1_legend,num_known_samples(radi_i)),'LineWidth',2)
        plot(all_heihgts, rmse_meas_radius_data2(:,radi_i),radius_colors(2)+plotLineStyle+radius_plot_type(2),'DisplayName', sprintf('%s (M = %d)',data2_legend,num_known_samples(radi_i)),'LineWidth',2)
        plot(all_heihgts, rmse_meas_radius_data3(:,radi_i),radius_colors(3)+plotLineStyle+radius_plot_type(3),'DisplayName', sprintf('%s (M = %d)',data3_legend,num_known_samples(radi_i)),'LineWidth',2)
        plot(all_heihgts, rmse_meas_radius_data4(:,radi_i),radius_colors(4)+plotLineStyle+radius_plot_type(4),'DisplayName', sprintf('%s (M = %d)',data4_legend,num_known_samples(radi_i)),'LineWidth',2)
        new_iter = new_iter + 1;
    end

    plot(all_heihgts, [sqrt(45.7237);sqrt(39.8541);sqrt(41.5870);sqrt(46.7718);sqrt(44.9827)], '--','DisplayName', 'Shadow Fading $\sigma$ (TRPL)', 'LineWidth',2,'color',[0.0, 0.0, 0.0])
    plot(all_heihgts, [sqrt(44.6398);sqrt(39.7005);sqrt(41.4070);sqrt(46.7749);sqrt(44.7768)], ':','DisplayName', 'Shadow Fading $\sigma$ (FSPL)', 'LineWidth',2,'color',[0.0, 0.0, 0.0])
    xlabel('Experiment, Altitude [m]')
    ylabel(y_label)
    grid on
    axis tight
    legend show
    legend('Interpreter', 'latex');
    fig_text = fig_id1 + y_label + " (R = " + num2str(current_heihgt) + " m)";
    fontsize(FigH,font_size,"points")
    box on
    ylim([5.4, 8.8])
    xticks(all_heihgts);  % Set the x-tick positions
    xticklabels({'A1, 30', 'A2, 50', 'A3, 70', 'A4, 90', 'A5, 110'});  % Set the custom string labels
end

