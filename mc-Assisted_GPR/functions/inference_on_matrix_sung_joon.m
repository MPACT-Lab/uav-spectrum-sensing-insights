function [sha_tar_est] = inference_on_matrix_sung_joon(RSRP_grid,lon_meshgrid,lat_meshgrid,lat_tar,lon_tar,method)
%UNTITLED Summary of this function goes here
%   linear', 'nearest', 'cubic', 'makima', or 'spline'
    sha_tar_est = interp2(lon_meshgrid,lat_meshgrid,RSRP_grid,lon_tar,lat_tar,method);
end
