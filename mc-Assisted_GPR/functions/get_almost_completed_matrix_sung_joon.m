function [RSRP_grid,Variance_grid,lon_meshgrid,lat_meshgrid] = get_almost_completed_matrix_sung_joon(lat_meas_h, lon_meas_h, h_meas_h, sha_meas_h, ...
    grid_dist, max_allowed_est_variance, kriging_method, var_sha, mu_sha, fixed_num_points, corr_coef )
    
    % 2D field constants
    R_earth = 6378.137 * 10^3 ; % [m]
    scaler = 111139 ;

    minVal = extend_boundary(min(sha_meas_h), 1.3);
    maxVal = extend_boundary(max(sha_meas_h), 1.3);

    added_term = var_sha*2;

    grid_res = grid_dist/scaler; % in latitude or longitude
    
    mX_range = -78.7002: grid_res: -78.6961;
    mY_range = 35.7232: grid_res: 35.7303;

    [lon_meshgrid, lat_meshgrid] = meshgrid(mX_range, mY_range);
    lon_tar = lon_meshgrid(:) ;
    lat_tar = lat_meshgrid(:) ;

    plot_inter_points = 1;
    plotted1 = 0;
    plotted2 = 0;

    dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas_h' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas_h' * pi/180) .* cos( (lon_tar - lon_meas_h')  * pi/180 )) ;
    %dist_ver_tar = abs(h_tar - h_meas_h'); %(num_tar, num_meas)

    est_sha_target = zeros(1, length(lon_tar))*nan;
    est_sha_target_variance = zeros(1, length(lon_tar))*nan;

    %grid_radius = grid_dist*sqrt(2);

    % 2D interpolation
    [~,index_sorted_dist] = sort(dist_hor_tar,2);

    if plot_inter_points>0 && plotted1<1
        %figure(repeat_id)
        show_scattered_samples_aerpaw_black_white(lon_meas_h, lat_meas_h, sha_meas_h, 192, 'o', 'Measurements (M=50)', "inter result")
        plotted1 = 1;
    end
    random_tar = randi([1,length(lat_meas_h)], 1, 1);

    for tar_idx = 1 : length(lat_tar)
        sha_meas_valid = sha_meas_h(index_sorted_dist(tar_idx, 1:fixed_num_points));
        num_meas_valid = fixed_num_points;

        %disp(num_meas_valid)
        lat_meas_valid = lat_meas_h(index_sorted_dist(tar_idx, 1:fixed_num_points));
        lon_meas_valid = lon_meas_h(index_sorted_dist(tar_idx, 1:fixed_num_points));
        h_meas_valid = h_meas_h(index_sorted_dist(tar_idx, 1:fixed_num_points));
        dist_hor_tar_valid = dist_hor_tar(tar_idx, index_sorted_dist(tar_idx, 1:fixed_num_points));
        %dist_ver_tar_valid = dist_ver_tar(tar_idx, index_sorted_dist(tar_idx, 1:fixed_num_points));

        if plot_inter_points>0 && tar_idx==random_tar && plotted2<1
            show_scattered_samples_aerpaw_black_white(lon_tar(tar_idx), lat_tar(tar_idx), 0, 192, '*', 'Unknown', "inter result")
            show_scattered_samples_aerpaw_black_white(lon_meas_valid, lat_meas_valid, sha_meas_valid, 192, 'filled', 'Neighbours (N=10)', "inter result")
            plotted2 = 1;
        end
    
        %num_meas_used(repeat_id, meas_i, radi_i, tar_idx) = num_meas_valid;
    
        % covariance matrix generation
        if kriging_method=="ordinary"
            left_mat = ones(num_meas_valid+1);
            left_mat(num_meas_valid+1, num_meas_valid+1) = 0;
        elseif kriging_method=="simple"
            left_mat = ones(num_meas_valid);
        else
            error("unknown kriging interpolation type. valid: ordinary or simple")
        end

        % auto correlation model
        for i=1:num_meas_valid
            for j=1:num_meas_valid
                dist_ver_ij = abs(h_meas_valid(i) - h_meas_valid(j));
                dist_hor_ij = R_earth * acos(sin(lat_meas_valid(i) * pi/180) * sin(lat_meas_valid(j) * pi/180)...
                    + cos(lat_meas_valid(i) * pi/180) * cos(lat_meas_valid(j) * pi/180) * cos(...
                    (lon_meas_valid(i) - lon_meas_valid(j))  * pi/180 )) ;
                left_mat(i,j) = (1-get_correlation_from_table(dist_hor_ij, dist_ver_ij, corr_coef)) * (var_sha); 
            end
        end
    
        tar_semivariogram = zeros(num_meas_valid, 1);
        for i=1:num_meas_valid
            tar_semivariogram(i, 1) = (1-get_correlation_from_table(dist_hor_tar_valid(i), 0, corr_coef)) * (var_sha);
        end
        if kriging_method=="ordinary"
            left_mat(1:num_meas_valid, 1:num_meas_valid) = left_mat(1:num_meas_valid, 1:num_meas_valid) + added_term;
            right_mat = [tar_semivariogram + added_term;1];
            lammda0 = left_mat \ right_mat ; %left_mat \ right_mat ;
            lammda = lammda0(1:end-1) ;
            %lagrange_multi = lammda0(end) ; % lagrange multiplier
            w_kri_est = real(lammda')*sha_meas_valid ;
%                 if(isnan(w_kri_est))
%                     disp(lammda0)
%                     disp(tar_semivariogram)
%                     disp(left_mat(1:num_meas_valid, 1:num_meas_valid))
%                     disp(left_mat)
%                     
% 
%                     for i=1:num_meas_valid
%                         for j=1:num_meas_valid
%                             %dist_hor_ij = R_earth * acos(sin(lat_meas_valid(i) * pi/180) * sin(lat_meas_valid(j) * pi/180)...
%                                 %+ cos(lat_meas_valid(i) * pi/180) * cos(lat_meas_valid(j) * pi/180) * cos(...
%                                 %(lon_meas_valid(i) - lon_meas_valid(j))  * pi/180 )) ;
%                             disp(dist_hor_ij)
%                             %left_mat(i,j) = (1-get_correlation_from_table(dist_hor_ij,0,corr_coef)) * (var_sha); 
%                         end
%                     end
%                     show_scattered_samples_aerpaw(lon_tar(tar_idx), lat_tar(tar_idx), 0, 192, '*', 'target', "inter result")
%                     show_scattered_samples_aerpaw(lon_meas_valid, lat_meas_valid, sha_meas_valid, 192, 'filled', 'neighbours', "inter result")
%                     sdgklf = kjdgjo+1
%                 end
            est_variance = real(lammda0')*[tar_semivariogram;1] ;
            %disp(est_variance)
            if (w_kri_est>minVal) && (w_kri_est<maxVal) && (est_variance<max_allowed_est_variance) && (est_variance>0)
                est_sha_target(tar_idx) = w_kri_est;
                est_sha_target_variance(tar_idx) = est_variance;
            end
        elseif kriging_method=="simple"
            left_mat(1:num_meas_valid, 1:num_meas_valid) = var_sha - left_mat(1:num_meas_valid, 1:num_meas_valid);
            tar_covar = var_sha - tar_semivariogram;
            lammda = left_mat \ tar_covar ; %left_mat \ right_mat ;
            w_kri_est = real(lammda')*(sha_meas_valid - mu_sha) + mu_sha ;
            est_variance = var_sha - real(lammda')*tar_covar ;
            if (w_kri_est>minVal) && (w_kri_est<maxVal) && (est_variance<max_allowed_est_variance) && (est_variance>0)
                est_sha_target(tar_idx) = w_kri_est;
                est_sha_target_variance(tar_idx) = est_variance;
            end
        else
            error("unknown kriging interpolation type. valid: ordinary or simple")
        end
    end

    RSRP_grid = reshape(est_sha_target, size(lon_meshgrid));

    Variance_grid = reshape(est_sha_target_variance, size(lon_meshgrid));

    function r=get_correlation_from_table(dis_hor, dis_ver, corr_coef)
        % correlation_f
        r = (corr_coef(1)*exp(-(corr_coef(2)*dis_hor)) + (1-corr_coef(1))*exp(-(corr_coef(3)*dis_hor)))*exp(-corr_coef(4)*dis_ver); %using data points hor 250 ver 50
    end
end

