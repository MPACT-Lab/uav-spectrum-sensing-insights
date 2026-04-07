clc;
clear;
close all;

addpath("functions\")
addpath("../effect_ant_patt/patts\")

all_heihgts = [110];
exp_name = "_" + num2str(all_heihgts(1))+"m"; % change the mat_name as well
% policy: 30m (train ant patt) -> 30 m (train corr profile)
% policy: 50m (train ant patt) -> 50 m (train corr profile)
% policy: 70m (train ant patt) -> 70 m (train corr profile)
% policy: 90m (train ant patt) -> 90 m (train corr profile)
% policy: 110m (train ant patt) -> 110 m (train corr profile)
offset = [-78.698 35.727] ;
scaler = 111139 ;
%h_BS = 10 ; % [m]
%lat_BS = 35.727451 ; % [degree]
%lon_BS = -78.695974 ; % [degree]
R_earth = 6378.137 * 10^3 ; % [m]
%c = 3.0e8;
%f0 = 3.5e9;
%lamda = c/f0;

lat_all = [];
lon_all = [];
h_all = [];
sha_all = [];

is_smooth = "_smooth"; % "_smooth" or "" % smooth for Kriging
mat_name = "flipx_antenna_tx_rx_updated_sung_"+num2str(all_heihgts(1))+"m"+ is_smooth +".mat"; % for Kriging

%% load shadow fading data
for i = 1:length(all_heihgts)
    load(sprintf('../data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_%d_meas_rad.mat',all_heihgts(i))) % 'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s', 'elev_s', height_s
    [sha_two_s, sha_free_s] = update_PL_sung_new_rad_pat(sha_two_s, sha_free_s, azim_s, elev_s, ori_s, mat_name);
    lat_all = [lat_all; lat_s];
    lon_all = [lon_all; lon_s];
    h_all = [h_all; height_s];
    sha_all = [sha_all; sha_two_s];
    if(0)
        figure(1)
        hold on
        scatter3( (lon_s - offset(1))*scaler,(lat_s - offset(2))*scaler,height_s,10,RSRP_s)
        if i==5
            colormap(jet);
            colorbar;
            caxislim=[-160 -50];
            caxis(caxislim)
            xlabel('X [m]')
            ylabel('Y [m]')
            zlabel('Height [m]')
            zlim([0 120]);
            hcb=colorbar;
            colorTitleHandle = get(hcb,'Title');
            titleString = 'RSRP (dBm)';
            set(colorTitleHandle ,'String',titleString);
            set(gcf,'color','w');
            grid on;
            view([45 45])
        end
    end
end
%% 
% finding the correlation in 3D
rng(100)
N_seed = 4000 ;
d_int = 2 ;
d_max = 500 ;
d_grid = [0:d_int:d_max]' ;

h_int = 20 ;
h_max = 10 ;
h_grid = [-10:h_int:h_max]' ;

s_idx = randi(length(sha_all),round(N_seed*1.8),1) ;
s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
s_idx = s_idx(1:N_seed) ;

mu_sha = mean(sha_all(s_idx)) ;
var_sha = var(sha_all(s_idx)) ;

dist_hor = zeros(N_seed, N_seed);
dist_ver = zeros(N_seed, N_seed);
corr_all = zeros(N_seed, N_seed);
for ii=1:N_seed
    dist_ver(:,ii) = abs(h_all(s_idx) - h_all(s_idx(ii)));
    dist_hor(:,ii) = R_earth .* acos(sin(lat_all(s_idx(ii)) * pi/180) .* sin(lat_all(s_idx) * pi/180) + cos(lat_all(s_idx(ii)) * pi/180) .* cos(lat_all(s_idx) * pi/180) .* cos( (lon_all(s_idx) - lon_all(s_idx(ii)))  * pi/180 )) ;
    corr_all(:,ii) = ( sha_all(s_idx(ii)) - mu_sha ) .* ( sha_all(s_idx) - mu_sha ) / var_sha ;
    %corr_s(:,ii) = 1 - (sha_s(s_idx(ii)) - sha_s(s_idx)).^2 / (2*var_sha);
end


dist_hor2 = dist_hor(:);
dist_ver2 = dist_ver(:);
corr_all2 = corr_all(:);
if(0)
    figure(2)
    hold on
    scatter3( dist_hor2, dist_ver2, corr_all2, 5, corr_all2)
end

%% 

corr_grid = zeros(length(h_grid)-1, length(d_grid)-1);
N_grid = zeros(length(h_grid)-1, length(d_grid)-1);
hor_dist_grid = zeros(length(h_grid)-1, length(d_grid)-1);
ver_dist_grid = zeros(length(h_grid)-1, length(d_grid)-1);
for jj = 1: length(h_grid)-1
    for ii = 1:length(d_grid)-1
        corr_tmp1 = corr_all2(dist_hor2>=d_grid(ii) & dist_hor2<d_grid(ii+1)...
            & dist_ver2>=h_grid(jj) & dist_ver2<h_grid(jj+1)) ;
        N_grid(jj,ii) = length(corr_tmp1) ;
        corr_grid(jj,ii) = mean(corr_tmp1) ;
        hor_dist_grid(jj,ii) = (d_grid(ii) + d_grid(ii+1))/2;
        ver_dist_grid(jj,ii) = (h_grid(jj) + h_grid(jj+1))/2;
    end
end
% 
% x_f = d_grid(2:end) ;
% y_f = corr_grid(1:end) ;

hor_dist_list = hor_dist_grid(:);
ver_dist_list = ver_dist_grid(:); 
corr_list = corr_grid(:);

disp('done')
%% 
ft = fittype( '(a*exp(-(b*x)) + (1-a)*exp(-(c*x)))', 'independent', {'x'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', [1,2,2], 'Lower',[0,0,0]);
correlation_f = fit([hor_dist_list],corr_list,ft,opts)
save(sprintf("profiles/correlation_table_3D_a2g"+is_smooth+exp_name+".mat"),'hor_dist_list','ver_dist_list', 'corr_list', 'mu_sha', "var_sha", 'd_int', 'h_int', 'correlation_f')
%% 

figure(3)
hold on
scatter( hor_dist_list, corr_list)
xlabel('Horizontal distance [m]')
ylabel('Correlation')
set(gcf,'color','w');
grid on;


fitted_f = @(x) (correlation_f.a*exp(-(correlation_f.b*x)) + (1-correlation_f.a)*exp(-(correlation_f.c*x))); % transformed sha_all
fplot(fitted_f,[0,500])

 function v=shuffle(v)
     v=v(randperm(length(v)));
 end