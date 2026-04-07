clear; 
close all; 
addpath("raw_data\")
%%
UAV_height = 110; % [m]

switch UAV_height
    case(30)
        load plotdata1.mat ;
        h_UAV = 30 ; % [m]
    case(50)
        load plotdata2.mat ; 
        h_UAV = 50 ; % [m]
    case(70)
        load plotdata3.mat ;
        h_UAV = 70 ; % [m]
    case(90)
        load plotdata4.mat ;
        h_UAV = 90 ; % [m]
    case(110)
        load plotdata5.mat ;
        h_UAV = 110 ; % [m]
end

%% parameters
h_BS = 10 ; % [m]
lat_BS = 35.727451 ; % [degree]
lon_BS = -78.695974 ; % [degree]
R_earth = 6378.137 * 10^3 ; % [m]
c = 3.0e8;
f0 = 3.5e9;
lamda = c/f0;

%% Height correction
Z_dif = mZ(end) - mZ(1) ;
N_sam = length(mZ) ;
Z_del = Z_dif / N_sam ;
h_UAV_new = mZ - Z_del*[1:N_sam]' ;

% orientation calculation
[ori_all, speed_all, acceleration_all] = calc_orientation_and_speed(mX, mY, timestamp_meas);

if(1)
    figure(1)
    hold on
    plot((timestamp_meas-timestamp_gps(1)),mZ,'bx','MarkerSize',2,'DisplayName', sprintf(['GPS']))
    plot((timestamp_meas-timestamp_gps(1)),h_UAV_new,'ro','MarkerSize',2,'DisplayName', sprintf(['corrected']))
    grid on;
    xlabel('Time (s)')
    ylabel('Altitude(m)')
    set(gcf,'color','w');
    legend show
end

%% RSRP correction
RSRP_offset = 98 ; % [dB]
RSRP_new = RSRP2' - RSRP_offset ; % [dBm]
RSRP_idx = RSRP_new>-140 ;

%% flight time
time = timestamp_meas-timestamp_gps(1) ;
idx_tmp = find(mZ > UAV_height - 0.0) ; % flight time without departure and landing
if (UAV_height == 110)
    idx_tmp = find(mZ > UAV_height - 0.5) ;
end
time_idx_tmp = [min(idx_tmp):1: max(idx_tmp)] ;
time_idx = zeros(size(time)) ;
time_idx(time_idx_tmp) = 1 ;
time_idx = logical(time_idx) ;

%% effective data index
tot_idx = RSRP_idx & time_idx ;

%% Distance
dist_2D_new = R_earth .* acos(sin(lat_BS * pi/180) .* sin(mY * pi/180) + cos(lat_BS * pi/180) .* cos(mY * pi/180) .* cos( (mX - lon_BS)  * pi/180 )) ;
dist_3D = sqrt( (dist_2D_new).^2 + (h_UAV_new - h_BS).^2 ) ;

if(1)  
    figure(2)
    hold on
    plot(time,dist_2D_new,'bx','MarkerSize',2,'DisplayName', sprintf(['2D distance']))
    plot(time,dist_3D,'ro','MarkerSize',2,'DisplayName', sprintf(['3D distance']))   
    grid on;
    xlabel('Time (s)')
    ylabel('3D distance (m)')
    set(gcf,'color','w');
    legend show
end

%% Elevation angle
elev = atan( ( h_UAV_new - h_BS )./dist_2D_new) ;

%% Azimuth angle

for ii=1:N_sam
    azim(ii,1) = atan2( sin( pi/180*(lon_BS - mX(ii)) ) * cos( pi/180*lat_BS ), cos( pi/180*mY(ii) ) * sin( pi/180*lat_BS ) - sin( pi/180*mY(ii) ) * cos( pi/180*lat_BS ) * cos( pi/180*(lon_BS - mX(ii)) ) ) ;
end


radiation_pat;
load rad_pat.mat


P_Tx_dBm = 10 ; % [dBm]

G_t_max = 0;
G_r_max = 1 ;

phi_os = 88; % 0 ~ 119 % 90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ant_opt = 'measured' ; % 'measured', 'constant', 'dipole'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(ant_opt)
    case('measured')       
        % Antenna pattern from the new measurement
        [~,idx_theta] = min(abs(elev - theta2'*pi/180)') ;
        [~,idx_phi] = min(abs(azim - circshift(phi2, phi_os)*pi/180)') ;

        azim_plus_pi = azim - (azim>=0).*pi + (azim<0).*pi;
        [~,idx_phi2] = min(abs(-azim_plus_pi  - phi2*pi/180)') ; % eqn (58) on latex document
        
        for ii=1:N_sam
            G_t(ii,1)= G_t_dB(idx_theta(ii), idx_phi(ii)) ; % [dBi]
        end
       G_r = G_r_dB(idx_theta,1) + G_r_max ; % [dBi]
    case('constant')
        % Antenna pattern constant gain
        G_t = 0 ;
        G_r = 0 ;
    case('dipole')
        % Antenna pattern dipole shape
        const_gain = 0 ;
        d_len = 0.5*lamda;
        
        G_t_lin = (cos(pi*d_len/lamda*cos(pi/2 - elev))-cos(pi*d_len/lamda))./sin(pi/2 - elev);
        G_r_lin = (cos(pi*d_len/lamda*cos(pi/2 - elev))-cos(pi*d_len/lamda))./sin(pi/2 - elev);
        
        G_t = 20*log10(G_t_lin) + const_gain ;
        G_r = 20*log10(G_r_lin) ;
end

%% Two-ray ground reflection pathloss model
rho = 0.005 ;
eps = 15 ;
eps_0 = eps - j*60*rho*lamda ;

th_r = atan( ( h_BS + h_UAV_new ) ./ (dist_2D_new) ) ;
th_r_2 = 2*pi - th_r ; % not used
th_r_22 = - th_r ;

r_1 = h_BS ./ sin(th_r) ;
r_2 = h_UAV_new ./ sin(th_r) ;
r_3 = dist_3D ;

phi_del = 2*pi*( r_1 + r_2 - r_3)/lamda ;

z_V = sqrt(eps_0 - cos(th_r).^2) / eps_0 ;
R_V = (sin(th_r) - z_V) ./ (sin(th_r) + z_V) ;

G_t_l = 10.^(G_t/10) ;
G_r_l = 10.^(G_r/10) ;

[~,idx_Gt_r] = min(abs(th_r_2 - G_t_grid')') ;
[~,idx_Gr_r] = min(abs(th_r - G_r_grid')') ;


G_t_r_dB = G_t_table(idx_Gt_r) ; % [dBi]
G_r_r_dB = G_r_table(idx_Gr_r) + G_r_max ; % [dBi]

G_t_r = 10.^(G_t_r_dB/10) ;
G_r_r = 10.^(G_r_r_dB/10) ;

unknown_offset = -5 ; % previously it was -5 (2/4/2026)

P_Tx_2_dBm = P_Tx_dBm + unknown_offset ;
P_Tx_2 = 10.^(P_Tx_2_dBm/10) ;

RSRP_PL_two_lin = P_Tx_2 * (lamda/(4*pi))^2 * abs( sqrt(G_t_l) .* sqrt(G_r_l) ./ ( r_3 ) + R_V .* sqrt(G_t_r) .* sqrt(G_r_r) ./ ( r_1 + r_2) .* exp( -j * phi_del ) ).^2 ;
RSRP_PL_free_lin = P_Tx_2 * (lamda/(4*pi))^2 * abs( sqrt(G_t_l) .* sqrt(G_r_l) ./ ( r_3 ) ).^2 ;

RSRP_PL_two = 10*log10(RSRP_PL_two_lin) ;
RSRP_PL_free = 10*log10(RSRP_PL_free_lin) ;

sha_two = RSRP_new - RSRP_PL_two ; 
sha_free = RSRP_new - RSRP_PL_free ; 

if(1)
    figure(3)
    hold on
    plot(time(tot_idx),RSRP_new(tot_idx),'bx','MarkerSize',3, 'DisplayName', sprintf(['Measurement (' num2str(h_UAV) 'm)']))
    plot(time(tot_idx),RSRP_PL_two(tot_idx),'cx','MarkerSize',3,'DisplayName', sprintf(['Two-ray (' num2str(h_UAV) 'm)']))
    plot(time(tot_idx),RSRP_PL_free(tot_idx),'rx','MarkerSize',3,'DisplayName', sprintf(['Free space (' num2str(h_UAV) 'm)']))   
    grid on;
    xlabel('Time (s)')
    ylabel('RSRP (dBm)')
    set(gcf,'color','w');
    legend show
    
    % elev grid to model shadow fading based on elevation angle
    d_elev_me = 10 ;
    elev_grid_me = (5:d_elev_me:65)*pi/180 ;
    sha_s = sha_two(tot_idx);
    elev_me = elev(tot_idx) ;
    % for later selecting the elev grid for all valid data
    elev_dif_table = abs(elev_me - elev_grid_me);
    [~,choose_grid_idx] = min(elev_dif_table');
    
    for ii=1:length(elev_grid_me)
        idx_tmp1 = (elev_me>=elev_grid_me(ii)-d_elev_me*pi/360)&...
            (elev_me<=elev_grid_me(ii)+d_elev_me*pi/360) ;
        mu_sha_grid(ii,1) = mean(sha_s(idx_tmp1)) ;
        var_sha_grid(ii,1) = var(sha_s(idx_tmp1)) ;
    end
    
    figure(4)
    hold on;
    semilogx(dist_3D(tot_idx),RSRP_PL_two(tot_idx),'cx','MarkerSize',3,'DisplayName', sprintf(['Two-ray (' num2str(h_UAV) 'm)']))
    semilogx(dist_3D(tot_idx),RSRP_PL_free(tot_idx),'rx','MarkerSize',3,'DisplayName', sprintf(['free space (' num2str(h_UAV) 'm)']))
    grid on;
    %xlim([0 500])
    axis tight
    xlabel('3D distance [m]')
    ylabel('RSRP [dBm]')
    set(gcf,'color','w');
    legend show

    figure(40)
    hold on;
    semilogx(dist_3D(tot_idx),RSRP_new(tot_idx),'bx','MarkerSize',3, 'DisplayName', sprintf(['Measurement (' num2str(h_UAV) 'm)']))
    semilogx(dist_3D(tot_idx),RSRP_PL_two(tot_idx),'cx','MarkerSize',3,'DisplayName', sprintf(['Two-ray (' num2str(h_UAV) 'm)']))
    semilogx(dist_3D(tot_idx),RSRP_PL_free(tot_idx),'rx','MarkerSize',3,'DisplayName', sprintf(['free space (' num2str(h_UAV) 'm)']))
    grid on;
    %xlim([0 500])
    axis tight
    xlabel('3D distance [m]')
    ylabel('RSRP [dBm]')
    set(gcf,'color','w');
    legend show

    figure(6)
    hold on;
    semilogx(dist_3D(tot_idx),phi_del(tot_idx),'cx','MarkerSize',3,'DisplayName', sprintf(['Two-ray (' num2str(h_UAV) 'm)']))
    grid on;
    %xlim([0 500])
    axis tight
    xlabel('3D distance [m]')
    ylabel('Phase Offset [radian]')
    set(gcf,'color','w');
    legend show

end

%% saving the data

lat_s = mY(tot_idx) ;
lon_s = mX(tot_idx) ;
RSRP_s = RSRP_new(tot_idx);
elev_s = elev(tot_idx);
azim_s = azim(tot_idx);
ori_s = ori_all(tot_idx);
speed_s = speed_all(tot_idx);
acceler_s = acceleration_all(tot_idx);

sha_two_s = sha_two(tot_idx);
sha_free_s = sha_free(tot_idx);
RSRP_PL_two_s = RSRP_PL_two(tot_idx);
RSRP_PL_free_s = RSRP_PL_free(tot_idx);

%save(sprintf('filtererd_data_%d.mat',UAV_height),'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s')
%% Shadowing componet

x = sha_two(tot_idx);
pd = fitdist(x,'Normal') ; % 'ExtremeValue', 'Kernel', 'Normal'
x_values = -40:1:40;
y = pdf(pd,x_values);

gaussian = @(x,mu,sigma) (1/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma).^2/2))
skewedgaussian = @(x,mu,sigma,alpha) 2*gaussian(x,mu,sigma).*normcdf(alpha*x,alpha*mu,sigma)

if(1)
    figure(5)
    hold on
    histogram(x,'Normalization','pdf','DisplayName', sprintf(['Shadowing (' num2str(h_UAV) 'm)']))
    plot(x_values,y,'LineWidth',2,'DisplayName', sprintf(['Gaussian (' num2str(h_UAV) 'm)']))
%     plot(x_values, skewedgaussian(x_values,pd.mu+8,pd.sigma*1.5, -2),'LineWidth',2,'DisplayName', sprintf(['Left-skewed Gaussian (' num2str(h_UAV) 'm)'])) % 30 m
    plot(x_values, skewedgaussian(x_values,pd.mu+7,pd.sigma*1.4, -2),'LineWidth',2,'DisplayName', sprintf(['Left-skewed Gaussian (' num2str(h_UAV) 'm)'])) % 50 m | 70 m | 90 m | 110 m
    xlabel('Power (dB)')
    ylabel('Probability density function')
    xlim([-40 40])
    grid on;
    set(gcf,'color','w');
    legend show
end

%% UAV speed by GPS logs
if(0)
    lat_GPS_1 = GPSy(1:end-1) ;
    lat_GPS_2 = GPSy(2:end) ;
    
    lon_GPS_1 = GPSx(1:end-1) ;
    lon_GPS_2 = GPSx(2:end) ;
    
    t_del = timestamp_gps(2:end) - timestamp_gps(1:end-1) ;
    
    
    dist_2D_new = R_earth .* acos(sin(lat_GPS_2 * pi/180) .* sin(lat_GPS_1 * pi/180) + cos(lat_GPS_2 * pi/180) .* cos(lat_GPS_1 * pi/180) .* cos( (lon_GPS_1 - lon_GPS_2)  * pi/180 )) ;
    
    UAV_speed = dist_2D_new./t_del ;
    time_GPS = timestamp_gps(1:end-1)-timestamp_gps(1) ;
    
    v = 5;
    f_d_norm = v*f0/c/15000 ;
    
    figure(6)
    hold on
    plot(time_GPS,UAV_speed,'bx','MarkerSize',3, 'DisplayName', sprintf(['Measurement (' num2str(h_UAV) 'm)']))
    xlabel('Time (s)')
    ylabel('UAV speed (m/s)')
    grid on;
    set(gcf,'color','w');
    legend show
end
[];

%% Trajectory
if(1)
    figure(7)
    offset = [-78.698 35.727] ;
    scaler = 111139 ;
    hold on
    scatter3( (mX - offset(1))*scaler,(mY - offset(2))*scaler,h_UAV_new,10,RSRP_new)
    viscircles([0,0],[40],'Color','k')
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
if(1)
    figure(70)
    offset = [-78.698 35.727] ;
    scaler = 111139 ;
    hold on
    scatter3( (mX - offset(1))*scaler,(mY - offset(2))*scaler,h_UAV_new,10,ori_all*180/pi)
    viscircles([0,0],[40],'Color','k')
    colormap(jet);
    colorbar;
    caxislim=[-160 -50];
    %caxis(caxislim)
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Height [m]')
    zlim([0 120]);
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Yaw (degree)';
    set(colorTitleHandle ,'String',titleString);
    set(gcf,'color','w');
    grid on;
    view([45 45])
end 
if(1)
    figure(71)
    offset = [-78.698 35.727] ;
    scaler = 111139 ;
    hold on
    scatter3( (mX - offset(1))*scaler,(mY - offset(2))*scaler,h_UAV_new,10,speed_all)
    viscircles([0,0],[40],'Color','k')
    colormap(jet);
    colorbar;
    caxislim=[-160 -50];
    %caxis(caxislim)
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Height [m]')
    zlim([0 120]);
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Spped (m/s)';
    set(colorTitleHandle ,'String',titleString);
    set(gcf,'color','w');
    grid on;
    view([45 45])
end 
%% Spatial correlation (2D)
if(1)
    rng(100)
    N_seed = 2000 ;
    d_int = 2 ;
    d_max = 500 ;
    d_grid = [0:d_int:d_max]' ;
    
    sha_s = sha_two(tot_idx) ;
    lat_s = mY(tot_idx) ;
    lon_s = mX(tot_idx) ;
    height_s = h_UAV_new(tot_idx);
    dist_2D_new_s = dist_2D_new(tot_idx);
    dist_3D_s = dist_3D(tot_idx);
    
    s_idx = randi(length(sha_s),round(N_seed*1.8),1) ;
    s_idx = unique(s_idx) ;
    s_idx = s_idx(1:N_seed) ;
    
    mu_sha = mean(sha_s(s_idx)) ;
    var_sha = var(sha_s(s_idx)) ;

    
    for ii=1:N_seed
        dist_s5(:,ii) = 111139.*sqrt((lon_s(s_idx(ii))-lon_s(s_idx)).^2 + (lat_s(s_idx(ii))-lat_s(s_idx)).^2);
        dist_s(:,ii) = R_earth .* acos(sin(lat_s(s_idx(ii)) * pi/180) .* sin(lat_s(s_idx) * pi/180) + cos(lat_s(s_idx(ii)) * pi/180) .* cos(lat_s(s_idx) * pi/180) .* cos( (lon_s(s_idx) - lon_s(s_idx(ii)))  * pi/180 )) ;
        corr_s(:,ii) = ( sha_s(s_idx(ii)) - mu_sha ) .* ( sha_s(s_idx) - mu_sha ) / var_sha ;
        %corr_s(:,ii) = ( sha_s(s_idx(ii)) - mu_sha_grid(choose_grid_idx(s_idx(ii))) ) .* ( sha_s(s_idx) - mu_sha_grid(choose_grid_idx(s_idx)) ) ./ sqrt(var_sha_grid(choose_grid_idx(s_idx(ii))).*var_sha_grid(choose_grid_idx(s_idx))) ;
        %corr_s(:,ii) = 1 - (sha_s(s_idx(ii)) - sha_s(s_idx)).^2 / (2*var_sha);
    end
    
    [dist_s2,idx_ds2] = sort(dist_s(:)) ;
    dist_s2 = real(dist_s2) ;
    dist_s6 = real(dist_s5(idx_ds2)) ;
    corr_s2 = real(corr_s(idx_ds2)) ;
    figure(10)
    scatter(dist_s2,dist_s6)
    
    for ii=1:length(d_grid)-1
        corr_tmp1 = corr_s2(dist_s2>=d_grid(ii)&dist_s2<d_grid(ii+1)) ;
        N_grid(ii,1) = length(corr_tmp1) ;
        corr_grid(ii,1) = mean(corr_tmp1) ;
    end
    
    x_f = d_grid(2:end) ;
    y_f = corr_grid(1:end) ;
%     x_f_filtered = x_f(y_f > 0);
%     y_f_filtered = y_f(y_f > 0);
%     x_f = x_f_filtered;
%     y_f = y_f_filtered;
    % save data
    save(sprintf('processed_data/filtererd_data_with_dist_ori_%d_meas_rad.mat',UAV_height),'lat_s','lon_s','RSRP_s','sha_two_s','sha_free_s','RSRP_PL_two_s','RSRP_PL_free_s','x_f','y_f','elev_s','azim_s','ori_s','speed_s','acceler_s','height_s', 'dist_2D_new_s', 'dist_3D_s')
    % original
    ft = fittype( 'exp(-a*x)', 'independent', 'x', 'dependent', 'y' );
    %ft = fittype( 'exp(-(a*x+c))+d', 'independent', 'x', 'dependent', 'y' );
    ft = fittype( 'a*exp(-(b*x)) + (1-a)*exp(-(c*x))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', [1,2,2], 'Lower',[0,0,0] );
    %f = fit(x_f,y_f,ft,opts)
    
    % Plot
    figure(8)
    hold on;
    %plot(f,x_f, y_f)
    scatter(x_f, y_f)
    set(gcf,'color','w');
    grid on;
    ylim([-0.2 1])
    xlabel('Distance [m]')
    ylabel('Corrleation')
end

function [ori_all, speed_all, acceleration_all] = calc_orientation_and_speed(lon_all, lat_all, timestmp)
    %CALC_ORIENTATION Summary of this function goes here
    %   Detailed explanation goes here
    space = 40; % around 1 second
    prev_weight = 0.0;
    R_earth = 6378.137 * 10^3 ; % [m]
    ori_all = zeros(size(lon_all));
    speed_all = zeros(size(lon_all));
    acceleration_all = zeros(size(lon_all));
    ori_all(1:space) = pi/2;
    dist_all = zeros(1,length(ori_all) - space);
    for i = space+1: length(lon_all)
        %distX = cos( pi/180*lat_all(i) ) * sin( pi/180*lat_all(i-space) ) - sin( pi/180*lat_all(i) ) * cos( pi/180*lat_all(i-space) ) * cos( pi/180*(lon_all(i-space) - lon_all(i)) );
        %distY = sin( pi/180*(lon_all(i-space) - lon_all(i)) ) * cos( pi/180*lat_all(i-space) );
        distX = cos(pi/180*(lon_all(i))).*cos(pi/180*(lat_all(i))) - cos(pi/180*(lon_all(i-space))).*cos(pi/180*(lat_all(i-space)));
        distY = sin(pi/180*(lon_all(i))).*cos(pi/180*(lat_all(i))) - sin(pi/180*(lon_all(i-space))).*cos(pi/180*(lat_all(i-space)));
        dist_all(i) = R_earth .* acos(sin(lat_all(i-space) * pi/180) .* sin(lat_all(i) * pi/180) + cos(lat_all(i-space) * pi/180) .* cos(lat_all(i) * pi/180) .* cos( (lon_all(i) - lon_all(i-space))  * pi/180 )) ;
        current_ori = ori_all(i-1);
        %if dist_all(i) > 1
            current_ori = atan2(distY, distX);
            %current_ori = atan2(lat_all(i)-lat_all(i-space), lon_all(i)-lon_all(i-space));
        %end
        % get rid of round of error: discard angles near -180
%         if current_ori < -175*(pi/180)
%             current_ori = 177*(pi/180);
%         end
        ori_all(i) = ori_all(i-1) * prev_weight + current_ori * (1-prev_weight);
        speed_all(i) = dist_all(i)/ (timestmp(i) - timestmp(i-space));
        acceleration_all(i) = (speed_all(i) - speed_all(i-space))/ (timestmp(i) - timestmp(i-space));
    end
end

