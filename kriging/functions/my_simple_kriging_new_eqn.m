function [est_tar, est_tar_biased, est_variance, num_points] = my_simple_kriging_new_eqn(meas_correlation_all,...
    target_cross_correlation_all, index_valid_radi, meas_sha, var_sha, mu_sha, second_der_value, trans_raw)

    % input
    % meas_correlation_all (num_meas, num_meas)
    % target_cross_correlation_all (num_tar, num_meas)
    % index_valid_radi (num_tar, num_meas)

    num_test_points = size(target_cross_correlation_all,1);

    est_tar = ones(num_test_points,1) * NaN;
    est_tar_biased = ones(num_test_points,1) * NaN;
    est_variance = ones(num_test_points,1) * NaN;
    num_points = zeros(num_test_points,1);

    for tar_idx = 1 : num_test_points

        index_activated = squeeze(index_valid_radi(tar_idx, :));
        num_meas_valid = sum(index_activated>0);

        if num_meas_valid > 0
            meas_sha_activated = meas_sha(index_activated); % column matrix
    
            left_mat = ones(num_meas_valid);

            % covariance
            left_mat(1:num_meas_valid, 1:num_meas_valid) = meas_correlation_all(index_activated,index_activated);
    
            tar_correlation = target_cross_correlation_all(tar_idx,index_activated)'; % column matrix
    
            lammda = left_mat \ tar_correlation ; %left_mat \ right_mat ;
            
            % convert to raw
            est_tar_biased(tar_idx, 1) = trans_raw(real(lammda')*(meas_sha_activated - mu_sha) + mu_sha) ;
            est_variance(tar_idx, 1) = var_sha - real(lammda')*tar_correlation ;
            est_tar(tar_idx, 1) = est_tar_biased(tar_idx, 1) + (second_der_value/2) * (est_variance(tar_idx, 1) - real(lammda')*tar_correlation);
            
            num_points(tar_idx, 1) = num_meas_valid;

        end
    end
end



