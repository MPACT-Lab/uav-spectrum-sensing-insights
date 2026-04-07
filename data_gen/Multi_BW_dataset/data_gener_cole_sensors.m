clc;
clear;
close all;

addpath('functions\')
addpath('raw_data\')

%%
close all

fig_rows = 1;
fig_cols = 1;
i_run = 1;

use_ground_relfec = 0; % no use so far, will produce both, PL_free and two way

smooth = ""; % "" for Propagation modeling, "_smooth" for Kriging

altitude = "100m"; % test
fs = "fs20MHz"; % test
load(altitude + "_" + fs + ".mat") % basic data
[lat_all,lon_all,h_all,power_all,power_all2,speed_allX,speed_allY,speed_allZ,...
ori_all,roll_all,pitch_all,dist_2D_new,dist_3D,elev,azim,pitch_towards_src] = get_refined_variables_cole_asilomar(mX, mY, mZ, mSpeedX, mSpeedY, mSpeedZ, mpower1,...
mpower2, timestamp_power, mYaw, mRoll, mPitch);
speed_all = speed_allX;

%plot_fspl = 1;
%[RSRP_PL_free, RSRP_PL_two] = baseline_flipx_cole(altitude, fs, use_ground_relfec, 0, "", 0, 0);
% ugv_angle = 120*pi/180;
[RSRP_PL_free, RSRP_PL_two] = baseline_flipx_cole_unknown_orientation_fix2(altitude, fs, use_ground_relfec, 0, "", 0, 0);
%[legends1, legends2] = plot_elements_cole(fig_rows,fig_cols,i_run,power_all,RSRP_PL_free,RSRP_PL_two,dist_3D,ori_all, roll_all, pitch_all, azim, elev,lat_all, lon_all, h_all, speed_all,"k",'--', "",[], [],1,plot_fspl,rcv_azimuth_arm_len, pitch_towards_src,arm_len);

%save("processed_data/results_cole_with_tworay_tilt_src_alt_"+altitude+"_"+fs+".mat", 'power_all', 'ori_all', 'roll_all', 'pitch_all', 'azim', 'elev', 'dist_3D', 'dist_2D_new', 'RSRP_PL_free', 'RSRP_PL_two', 'lat_all', 'lon_all', 'h_all', 'speed_all',...
            %'pitch_towards_src');