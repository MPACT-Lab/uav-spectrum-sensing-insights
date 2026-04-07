clc;clear; close all;
addpath('functions\')
radius_colors = ["r", "b", "k", "m", "g", "c" ];
radius_plot_type = ["-s", "-<", "-x", "-d", "-o", "-v", "->"];
all_heights = [70];
all_radius = [70, 200];
font_size = 16;

load ./results/res_ok_lte_70m.mat % baseline

fig_id2 = "latex_fig_902_";
num_known_samples1 = num_known_samples;
data1_all_heights = error_all_2d; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
data1_legend = 'OK'; % change it
%%
load ./results/res_sk_lte_70m.mat

data2_all_heights = error_all_2d;
data2_legend = 'SK';

y_label = 'Median of RMSE [dB]'; % change it
y_label2 = 'Valid estimation ratio'; %


data1 = squeeze(data1_all_heights); % (iter,height=1,num_meas,radius,targets)
data2 = squeeze(data2_all_heights); % (iter,height=1,num_meas,radius,targets)

data1 = data1(:,:,:,1:20);
data2 = data2(:,:,:,1:20);

[rmse_meas_rad_data1, rmse_meas_rad_data2, valid_est1, valid_est2] = calc_rmse_median_sung_joon(data1,data2);

load ./results/res_gpr_lte_70m.mat

data3_all_heights = error_all_2d;
data3_legend = 'GPR';

load ../mc-Assisted_GPR/results/res_mc4_lte_70m.mat

data4_all_heights = error_all_2d;
data4_legend = 'MC-assisted GPR';
data5_all_heights = error_all_2d_gpr_grid;
data5_legend = 'GPR';

data3 = squeeze(data3_all_heights);
data3 = data3(:,:,:,1:20);


data4 = reshape(data4_all_heights, size(data4_all_heights,1), size(data4_all_heights,3), size(data4_all_heights,4), size(data4_all_heights,5));
data5 = reshape(data5_all_heights, size(data5_all_heights,1), size(data5_all_heights,3), size(data5_all_heights,4), size(data5_all_heights,5));


[rmse_meas_rad_data1, rmse_meas_rad_data3, valid_est1, valid_est3] = calc_rmse_median_sung_joon(data1,data3);
[rmse_meas_rad_data4, rmse_meas_rad_data5, valid_est1, valid_est3] = calc_rmse_median_sung_joon(data4,data5);

FigH = figure('Position', [100,100,680,610], 'visible','on');
hax1=axes;
hold on

new_iter = 1;
for rad_i = [1,2] % 70, 200 m
    if new_iter>1
        plot(num_known_samples1, rmse_meas_rad_data1(:,rad_i),radius_colors(1)+radius_plot_type(1),'DisplayName', sprintf('%s (R = %dm)',data1_legend,all_radius(rad_i)),'LineWidth',2)
        plot(num_known_samples1, rmse_meas_rad_data2(:,rad_i),radius_colors(2)+radius_plot_type(2),'DisplayName', sprintf('%s (R = %dm)',data2_legend,all_radius(rad_i)),'LineWidth',2)
    else
        plot(num_known_samples1, rmse_meas_rad_data1(:,rad_i),radius_colors(1)+'-'+radius_plot_type(1),'DisplayName', sprintf('%s (R = %dm)',data1_legend,all_radius(rad_i)),'LineWidth',2)
        plot(num_known_samples1, rmse_meas_rad_data2(:,rad_i),radius_colors(2)+'-'+radius_plot_type(2),'DisplayName', sprintf('%s (R = %dm)',data2_legend,all_radius(rad_i)),'LineWidth',2)
    end
    if new_iter>1
        plot(num_known_samples1, rmse_meas_rad_data5(:,1),radius_colors(5)+radius_plot_type(3),'DisplayName', sprintf('%s',data5_legend),'LineWidth',2)
        plot(num_known_samples1, rmse_meas_rad_data4(:,1),radius_colors(4)+radius_plot_type(4),'DisplayName', sprintf('%s',data4_legend),'LineWidth',2)
    end
    new_iter = new_iter + 1;
end

xlabel('Number of Sparse Samples, M')
ylabel(y_label)
grid on
axis tight
legend show
fontsize(FigH,font_size,"points")
box on
ylim([5.5, 6.2])

