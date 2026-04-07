function [RSRP_grid, Variance_grid] = get_almost_completed_matrix_sung_joon_efficient(lat_meas_h, lon_meas_h, h_meas_h, sha_meas_h, ...
    max_allowed_est_variance, kriging_method, var_sha, mu_sha, fixed_num_points_all, corr_coef, fixed_height, lon_tar, lat_tar, num_y, num_x, threed)
    
    % 2D field constants
    R_earth = 6378.137 * 10^3 ; % [m]

    minVal = extend_boundary(min(sha_meas_h), 1.3);
    maxVal = extend_boundary(max(sha_meas_h), 1.3);

    h_tar = ones(size(lon_tar)) * fixed_height;
    
    if threed
        meas_correlation_all = my_meas_correlation_3D(lat_meas_h,lon_meas_h, h_meas_h, corr_coef, var_sha); %(num_meas, num_meas)
        target_cross_correlation_all = my_target_cross_correlation_3D(lat_meas_h,lon_meas_h, h_meas_h, lat_tar, lon_tar, h_tar, corr_coef, var_sha); %(num_tar, num_meas)
    else
        meas_correlation_all = my_meas_correlation_2D(lat_meas_h,lon_meas_h, h_meas_h, corr_coef, var_sha); %(num_meas, num_meas)
        target_cross_correlation_all = my_target_cross_correlation_2D(lat_meas_h,lon_meas_h, h_meas_h, lat_tar, lon_tar, h_tar, corr_coef, var_sha); %(num_tar, num_meas)
    end
    dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas_h' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas_h' * pi/180) .* cos( (lon_tar - lon_meas_h')  * pi/180 )) ;
    %dist_ver_tar = abs(h_tar - h_meas_h'); %(num_tar, num_meas)

    % 2D interpolation
    [~,index_sorted_dist] = sort(dist_hor_tar,2); % (num_tar, num_meas)
    
    est_sha_target = zeros(length(fixed_num_points_all), length(lon_tar))*nan;
    est_sha_target_variance = zeros(length(fixed_num_points_all), length(lon_tar))*nan;

    for num_points_idx = 1 : length(fixed_num_points_all)
        fixed_num_points = fixed_num_points_all(num_points_idx);

        index_valid_closest = index_sorted_dist(:, 1:fixed_num_points);

        if kriging_method=="ordinary"
                [est_tar, est_variance, ~] = my_oridinary_kriging(meas_correlation_all,...
        target_cross_correlation_all, index_valid_closest, sha_meas_h, var_sha);
        elseif kriging_method=="simple"
                [est_tar, est_variance, ~] = my_simple_kriging(meas_correlation_all,...
        target_cross_correlation_all, index_valid_closest, sha_meas_h, var_sha, mu_sha);
        end
        
        is_est_logical = (est_tar>minVal) & (est_tar<maxVal) & (est_variance<max_allowed_est_variance) & (est_variance>0);
        est_sha_target(num_points_idx, is_est_logical) = est_tar(is_est_logical);
        est_sha_target_variance(num_points_idx, is_est_logical) = est_variance(is_est_logical);
    end

    RSRP_grid = reshape(est_sha_target, length(fixed_num_points_all), num_y, num_x);

    Variance_grid = reshape(est_sha_target_variance, length(fixed_num_points_all), num_y, num_x);

end

