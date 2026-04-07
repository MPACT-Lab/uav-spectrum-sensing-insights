function [sha_two_s, sha_free_s] = update_PL_sung_new_rad_pat(sha_two_s, sha_free_s, azim_s, elev_s, ori_s, mat_name)
    eval("load " + mat_name)

    % ori_s is the direction following trajectory
    ori_s = ori_s - pi/2;
    ori_s = ori_s + (ori_s<-pi).*pi*2;

    azim_wrt_uav_ori = -azim_s - ori_s*pi/180;
    azim_wrt_uav_ori = azim_wrt_uav_ori - (azim_wrt_uav_ori>=0).*pi + (azim_wrt_uav_ori<0).*pi;
    elev_wrt_uav_ori = elev_s;

    induced_rad_pat = rx_patt_with_placement_effect - G_r_dB;

    [~,idx_theta] = min(abs(elev_wrt_uav_ori - theta2'*pi/180)') ;
    [~,idx_phi] = min(abs(azim_wrt_uav_ori - phi2*pi/180)') ;
    for ii=1:length(sha_two_s)
        G_induced(ii,1)= induced_rad_pat(idx_theta(ii), idx_phi(ii)) ; % [dBi]
    end
    
    % we are ignoring the effect of changed radiation pattern for the
    % reflected rays

    sha_two_s = sha_two_s - G_induced;
    sha_free_s = sha_free_s - G_induced;

end