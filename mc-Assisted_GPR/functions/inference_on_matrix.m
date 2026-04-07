function [sha_tar_est] = inference_on_matrix(RSRP_grid_all_heights,lon_meshgrid,lat_meshgrid,lat_tar,lon_tar,h_tar,method,all_heights)
%UNTITLED Summary of this function goes here
%   linear', 'nearest', 'cubic', 'makima', or 'spline'
sha_tar_est = zeros(size(lat_tar))*NaN;
    for hi = 1:length(all_heights)
        RSRP_grid = squeeze(RSRP_grid_all_heights(hi,:,:));
        valid_h = h_tar > all_heights(hi) - 10 & h_tar < all_heights(hi) + 10;
        sha_tar_est(valid_h) = interp2(lon_meshgrid,lat_meshgrid,RSRP_grid,lon_tar(valid_h),lat_tar(valid_h),method);
    end
end
