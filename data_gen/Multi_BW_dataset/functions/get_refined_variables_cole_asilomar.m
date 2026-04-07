function [lat_all,lon_all,h_all,power_all,power_all2,speed_allX,speed_allY,speed_allZ,...
    ori_all,roll_all,pitch_all,dist_2D_new,dist_3D,elev,azim,pitch_towards_src] = get_refined_variables_cole_asilomar(mX, mY, mZ, mSpeedX, mSpeedY, mSpeedZ, mpower1,...
    mpower2, timestamp_power, mYaw, mRoll, mPitch)
%GET_REFINED_VARIABLES_COLE Summary of this function goes here
%   Detailed explanation goes here
    scaler = 111139 ;
    R_earth = 6378.137 * 10^3 ; % [m]

    ori_all = mYaw - 360 .* (mYaw>=180) + 360 .* (mYaw<-180);
    roll_all = mRoll - 360 .* (mRoll>=180) + 360 .* (mRoll<-180);
    pitch_all = mPitch - 360 .* (mPitch>=180) + 360 .* (mPitch<-180);
    
    %%
    height_avg = mode(round(mZ));
    
    valid_height = ((mZ>height_avg-5) & (mZ<height_avg+5));
    lat_all = mY(valid_height);
    lon_all = mX(valid_height);
    h_all = mZ(valid_height);
    power_all = mpower1(valid_height);
    power_all2 = mpower2(valid_height);
    speed_allX = mSpeedX(valid_height);
    speed_allY = mSpeedY(valid_height);
    speed_allZ = mSpeedZ(valid_height);
    ori_all = ori_all(valid_height);
    roll_all = roll_all(valid_height);
    pitch_all = pitch_all(valid_height);
    
    N_sam = length(power_all);
    
    origin_y=35.727451; % LW1
    origin_x=-78.695974;  % LW1
    
    lat_BS = origin_y;
    lon_BS = origin_x;
    
    %% 
    h_UAV_new = h_all;
    h_BS = 10;
    P_Tx_dBm = 100 ; % [dBm]
    c = 3.0e8;
    f0 = 3.32e9; % what is the actual f0 here
    lamda = c/f0;
    dist_2D_new = R_earth .* acos(sin(lat_BS * pi/180) .* sin(lat_all * pi/180) + cos(lat_BS * pi/180) .* cos(lat_all * pi/180) .* cos( (lon_all - lon_BS)  * pi/180 )) ;
    dist_3D = sqrt( (dist_2D_new).^2 + (h_UAV_new - h_BS).^2 ) ;
    unknown_offset = 0;
    P_Tx_2_dBm = P_Tx_dBm + unknown_offset ;
    P_Tx_2 = 10.^(P_Tx_2_dBm/10) ;
    
    %% Elevation angle
    elev = atan( ( h_UAV_new - h_BS )./dist_2D_new) ;
    
    %% Azimuth angle
    % cos*sin(lat_BS) - cos * sin(mY(ii)) % so it is respect to receiver i!! no
    % the data shows it is respect to the BS
    % guess
    % azim calculation was incorrect i think
    for ii=1:length(lon_all)
        azim(ii,1) = great_circle_azimuth(lat_BS,lon_BS,lat_all(ii),lon_all(ii));
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ant_opt = 'measured' ; % 'measured', 'constant', 'dipole'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hRx = 0.35;
    hBs = 10;
    pitch_towards_src = get_uav_body_reflection_asilomar(azim,elev,ori_all,pitch_all,roll_all,dist_2D_new,h_all);
end

