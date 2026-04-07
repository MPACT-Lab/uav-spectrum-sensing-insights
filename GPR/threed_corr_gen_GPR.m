clc;
clear;
close all;

all_heihgts = [90]; % change it
exp_name = "_" + num2str(all_heihgts(1))+"m";
offset = [-78.698 35.727] ;
scaler = 111139 ;
%h_BS = 10 ; % [m]f
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
% var (2-ray): 46.0687 39.9257 41.6243 46.7801 45.0187
% var (fspl): 44.6398 39.7005 41.4070 46.7749 44.7768
% load shadow fading data
for i = 1:length(all_heihgts)
    load(sprintf('../data_gen/LTE_dataset/processed_data/filtererd_data_with_dist_ori_%d_meas_rad.mat',all_heihgts(i))) % 'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s', 'elev_s', height_s
    lat_all = [lat_all; lat_s];
    lon_all = [lon_all; lon_s];
    h_all = [h_all; height_s];
    sha_all = [sha_all; sha_two_s]; % sha_two_s sha_free_s
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
var(sha_all)
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
%save(sprintf("profiles/correlation_table_3D_meas"+exp_name+"_rad.mat"),'hor_dist_list','ver_dist_list', 'corr_list', 'mu_sha', "var_sha", 'd_int', 'h_int', 'correlation_f')
%% 
%load('data/correlation_table_3D_meas.mat')


figure(3)
hold on
scatter( hor_dist_list, corr_list)
xlabel('Horizontal distance [m]')
ylabel('Correlation')
set(gcf,'color','w');
grid on;
%% 
N_seed = 100 ;

cum_dis = zeros(size(lat_all));
for jj=2:length(lat_all)
    cum_dis(jj) = cum_dis(jj-1) + sqrt((lat_all(jj)-lat_all(jj-1))^2 + (0.64*lon_all(jj)-lon_all(jj-1))^2);
end
equal_dis = linspace(cum_dis(1),cum_dis(end-1),N_seed);
s_idx = min(max(round(interp1(cum_dis,1:length(lat_all),equal_dis)),1),length(lat_all));

gap = s_idx(2) - s_idx(1);
random_starts = 0:10:gap-5;

dist_sam_multiple = zeros(length(random_starts),N_seed, N_seed);
sha_sam_multiple = zeros(length(random_starts),N_seed);
for i_random_start = 1:length(random_starts)
    s_idx = min(max(round(interp1(cum_dis,1:length(lat_all),equal_dis)+random_starts(i_random_start)),1),length(lat_all));
    sha_sam = sha_all(s_idx);
    dist_sam = zeros(N_seed, N_seed);
    for ii=1:N_seed
        dist_sam(:,ii) = R_earth .* acos(sin(lat_all(s_idx(ii)) * pi/180) .* sin(lat_all(s_idx) * pi/180) + cos(lat_all(s_idx(ii)) * pi/180) .* cos(lat_all(s_idx) * pi/180) .* cos( (lon_all(s_idx) - lon_all(s_idx(ii)))  * pi/180 )) ;
    end
    dist_sam = real(dist_sam);
    dist_sam_multiple(i_random_start,:,:) = dist_sam;
    sha_sam_multiple(i_random_start,:) = sha_sam';
end

% Define the residual function
% fun_covar = @(a,b,c,v,sig_gp) (v*a*exp(-(b*dist_sam)) + v*(1-a)*exp(-(c*dist_sam)) + sig_gp*eye(N_seed));
% fun = @(a,b,c,v,sig_gp) sha_sam'*((fun_covar(a,b,c,v,sig_gp))\sha_sam) + log(det(fun_covar(a,b,c,v,sig_gp)));
fun_covar = @(x,rand_i) (x(4)*x(1)*exp(-(x(2)*squeeze(dist_sam_multiple(rand_i,:,:)))) + x(4)*(1-x(1))*exp(-(x(3)*squeeze(dist_sam_multiple(rand_i,:,:)))) + x(5)*eye(N_seed));
fun = @(x,rand_i) sha_sam_multiple(rand_i,:)*((fun_covar(x,rand_i))\sha_sam_multiple(rand_i,:)') + log(det(fun_covar(x,rand_i)));

% first iteration is to bypass local minimums, from second iteration it is
% to make the prediction precise
n_elements = 11; % odd number to keep the center
a_hop = 0.3; % expected a_res = 0.2/10 = 0.02
b_hop = 0.15; % expected b_res = 0.2/10 = 0.02
c_hop = 0.15; % expected b_res = 0.2/10 = 0.02
v_hop = 20; % expected v_res = 40/10 = 4
sig_hop = 20; % expected v_res = 40/10 = 4
current_a = correlation_f.a; current_b = correlation_f.b; current_c = correlation_f.c; current_v = var_sha; current_sig = 0;
current_guess = [current_a, current_b, current_c, current_v, current_sig];
%fprintf("initial residual: %.2f\n",fun(current_guess))
disp(current_guess)
res_current = 0;
for i_random_start = 1:length(random_starts)
    res_current = res_current + fun(current_guess,i_random_start);
end
disp(res_current)

for i_outer=1:7
    a_min = max(0,current_a-a_hop); a_max = min(1,current_a+a_hop); a_res = (a_max-a_min)/(n_elements-1);
    a_list = linspace(a_min,a_max,n_elements);
    b_min = max(0,current_b-b_hop); b_max = min(1,current_b+b_hop); b_res = (b_max-b_min)/(n_elements-1);
    b_list = linspace(b_min,b_max,n_elements);
    c_min = max(0,current_c-c_hop); c_max = min(1,current_c+c_hop); c_res = (c_max-c_min)/(n_elements-1);
    c_list = linspace(c_min,c_max,n_elements);
    v_min = max(0,current_v-v_hop); v_max = current_v+v_hop; v_res = (v_max-v_min)/(n_elements-1);
    v_list = linspace(v_min,v_max,n_elements);
    sig_min = max(0,current_sig-sig_hop); sig_max = current_sig+sig_hop; sig_res = (sig_max-sig_min)/(n_elements-1);
    sig_list = linspace(sig_min,sig_max,n_elements);
    
    res_array = zeros(n_elements, n_elements, n_elements, n_elements, n_elements);
    for i1=1:n_elements
        for i2=1:n_elements
            for i3=1:n_elements
                for i4=1:n_elements
                    for i5=1:n_elements
                        res_array(i1,i2,i3,i4,i5) = 0;
                        for i_random_start = 1:length(random_starts)
                            res_array(i1,i2,i3,i4,i5) = res_array(i1,i2,i3,i4,i5) + fun([a_list(i1),b_list(i2),c_list(i3),v_list(i4),sig_list(i5)],i_random_start);
                        end
                    end
                end
            end
        end
    end
    [min_val, linear_idx] = min(res_array, [], 'all');
    [idx1, idx2, idx3, idx4, idx5] = ind2sub(size(res_array), linear_idx);
    current_a = a_list(idx1); current_b =  b_list(idx2); current_c = c_list(idx3); current_v = v_list(idx4); current_sig = sig_list(idx5);
    
    current_guess = [current_a, current_b, current_c, current_v, current_sig];
    disp(current_guess)
    disp(min_val)
    %fprintf("updated residual: %.2f\n",fun(current_guess))
    a_hop = a_res/2; b_hop = b_res/2; c_hop = c_res/2; v_hop = v_res/2; sig_hop = sig_res/2;
end
%% 
gpr_values = {};
gpr_values.a = current_a; gpr_values.b = current_b; gpr_values.c = current_c; gpr_values.v = current_v; gpr_values.sig = current_sig;
save(sprintf("gpr_coeffs/gpr_coeff"+exp_name+"_rad.mat"),'gpr_values')
%% 
N_seed = 100 ;
random_start = 40;
cum_dis = zeros(size(lat_all));
for jj=2:length(lat_all)
    cum_dis(jj) = cum_dis(jj-1) + sqrt((lat_all(jj)-lat_all(jj-1))^2 + (0.64*lon_all(jj)-lon_all(jj-1))^2);
end
equal_dis = linspace(cum_dis(1),cum_dis(end-1),N_seed);
s_idx = min(max(round(interp1(cum_dis,1:length(lat_all),equal_dis)+ random_start),1),length(lat_all));

sha_sam = sha_all(s_idx);
dist_sam = zeros(N_seed, N_seed);
for ii=1:N_seed
    dist_sam(:,ii) = R_earth .* acos(sin(lat_all(s_idx(ii)) * pi/180) .* sin(lat_all(s_idx) * pi/180) + cos(lat_all(s_idx(ii)) * pi/180) .* cos(lat_all(s_idx) * pi/180) .* cos( (lon_all(s_idx) - lon_all(s_idx(ii)))  * pi/180 )) ;
end
dist_sam = real(dist_sam);

% Define the residual function
% fun_covar = @(a,b,c,v,sig_gp) (v*a*exp(-(b*dist_sam)) + v*(1-a)*exp(-(c*dist_sam)) + sig_gp*eye(N_seed));
% fun = @(a,b,c,v,sig_gp) sha_sam'*((fun_covar(a,b,c,v,sig_gp))\sha_sam) + log(det(fun_covar(a,b,c,v,sig_gp)));
fun_covar = @(x) (x(4)*x(1)*exp(-(x(2)*dist_sam)) + x(4)*(1-x(1))*exp(-(x(3)*dist_sam)) + x(5)*eye(N_seed));
fun = @(x) sha_sam'*((fun_covar(x))\sha_sam) + log(det(fun_covar(x)));

% first iteration is to bypass local minimums, from second iteration it is
% to make the prediction precise
n_elements = 3; % odd number to keep the center
a_hop = 0.1; % expected a_res = 0.2/10 = 0.02
b_hop = 0.1; % expected b_res = 0.2/10 = 0.02
c_hop = 0.1; % expected b_res = 0.2/10 = 0.02
v_hop = 20; % expected v_res = 40/10 = 4
sig_hop = 20; % expected v_res = 40/10 = 4
current_a = correlation_f.a; current_b = correlation_f.b; current_c = correlation_f.c; current_v = var_sha; current_sig = 0;
current_guess = [current_a, current_b, current_c, current_v, current_sig];
fprintf("initial residual: %.2f\n",fun(current_guess))

for i_outer=1:5
a_min = max(0,current_a-a_hop); a_max = min(1,current_a+a_hop); a_res = (a_max-a_min)/(n_elements-1);
a_list = linspace(a_min,a_max,n_elements);
b_min = max(0,current_b-b_hop); b_max = min(1,current_b+b_hop); b_res = (b_max-b_min)/(n_elements-1);
b_list = linspace(b_min,b_max,n_elements);
c_min = max(0,current_c-c_hop); c_max = min(1,current_c+c_hop); c_res = (c_max-c_min)/(n_elements-1);
c_list = linspace(c_min,c_max,n_elements);
v_min = max(0,current_v-v_hop); v_max = current_v+v_hop; v_res = (v_max-v_min)/(n_elements-1);
v_list = linspace(v_min,v_max,n_elements);
sig_min = max(0,current_sig-sig_hop); sig_max = current_sig+sig_hop; sig_res = (sig_max-sig_min)/(n_elements-1);
sig_list = linspace(sig_min,sig_max,n_elements);

res_array = zeros(n_elements, n_elements, n_elements, n_elements, n_elements);
for i1=1:n_elements
    for i2=1:n_elements
        for i3=1:n_elements
            for i4=1:n_elements
                for i5=1:n_elements
                    res_array(i1,i2,i3,i4,i5) = fun([a_list(i1),b_list(i2),c_list(i3),v_list(i4),sig_list(i5)]);
                end
            end
        end
    end
end
[min_val, linear_idx] = min(res_array, [], 'all');
[idx1, idx2, idx3, idx4, idx5] = ind2sub(size(res_array), linear_idx);
current_a = a_list(idx1); current_b =  b_list(idx2); current_c = c_list(idx3); current_v = v_list(idx4); current_sig = sig_list(idx5);

current_guess
current_guess = [current_a, current_b, current_c, current_v, current_sig];
current_guess
fprintf("updated residual: %.2f\n",fun(current_guess))
a_hop = a_res/2; b_hop = b_res/2; c_hop = c_res/2; v_hop = v_res/2; sig_hop = sig_res/2;
end

%% 
N_seed = 100 ;
cum_dis = zeros(size(lat_all));
for jj=2:length(lat_all)
    cum_dis(jj) = cum_dis(jj-1) + sqrt((lat_all(jj)-lat_all(jj-1))^2 + (0.64*lon_all(jj)-lon_all(jj-1))^2);
end
equal_dis = linspace(cum_dis(1),cum_dis(end),N_seed);
s_idx = min(max(round(interp1(cum_dis,1:length(lat_all),equal_dis)),1),length(lat_all));

sha_sam = sha_all(s_idx);
dist_sam = zeros(N_seed, N_seed);
for ii=1:N_seed
    dist_sam(:,ii) = R_earth .* acos(sin(lat_all(s_idx(ii)) * pi/180) .* sin(lat_all(s_idx) * pi/180) + cos(lat_all(s_idx(ii)) * pi/180) .* cos(lat_all(s_idx) * pi/180) .* cos( (lon_all(s_idx) - lon_all(s_idx(ii)))  * pi/180 )) ;
end
dist_sam = real(dist_sam);

% Define the residual function
% fun_covar = @(a,b,c,v,sig_gp) (v*a*exp(-(b*dist_sam)) + v*(1-a)*exp(-(c*dist_sam)) + sig_gp*eye(N_seed));
% fun = @(a,b,c,v,sig_gp) sha_sam'*((fun_covar(a,b,c,v,sig_gp))\sha_sam) + log(det(fun_covar(a,b,c,v,sig_gp)));
% fun_covar = @(x) (x(4)*x(1)*exp(-(x(2)*dist_sam)) + x(4)*(1-x(1))*exp(-(x(3)*dist_sam)) + x(5)*eye(N_seed));
% fun = @(x) sha_sam'*((fun_covar(x))\sha_sam) + log(det(fun_covar(x)));
fun_covar = @(x) (x(1)*correlation_f.a*exp(-(correlation_f.b*dist_sam)) + x(1)*(1-correlation_f.a)*exp(-(correlation_f.c*dist_sam)) + x(2)*eye(N_seed));
fun = @(x) sha_sam'*((fun_covar(x))\sha_sam) + log(det(fun_covar(x)));

% Optimization options
options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt'...
    ,'OutputFcn', @(x, optimValues, state) disp_values(x, optimValues, state)...
    ); % 'Display', 'iter'

%current_guess = [correlation_f.a, correlation_f.b, correlation_f.c, var_sha-10, 1e-4];
current_guess = [var_sha-10, 20];
% lsqnonlin(func, initialization, lower_bound, upper bound)

optimal_soln = lsqnonlin(fun, current_guess, [], [], options);
current_guess
optimal_soln



fitted_f = @(x) (correlation_f.a*exp(-(correlation_f.b*x)) + (1-correlation_f.a)*exp(-(correlation_f.c*x))); % transformed sha_all
fplot(fitted_f,[0,500])

%% 
fun = @(x) ((x-4)^2);

% Optimization options
options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt',...
    'OutputFcn', @(x, optimValues, state) disp_values(x, optimValues, state)); % 'Display', 'iter'

current_guess = [0];
% lsqnonlin(func, initialization, lower_bound, upper bound)

optimal_soln = lsqnonlin(fun, current_guess, [], [], options);
current_guess
optimal_soln

 function v=shuffle(v)
     v=v(randperm(length(v)));
 end

 % This function runs at every iteration
function stop = disp_values(x, optimValues, state)
    stop = false; % Tells MATLAB to keep going
    fprintf('Iteration: %d | Current x: [%s]\n', optimValues.iteration, num2str(x));
end