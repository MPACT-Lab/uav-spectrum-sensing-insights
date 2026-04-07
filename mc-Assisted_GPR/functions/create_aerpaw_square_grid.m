function [lon_tar,lat_tar,lon_meshgrid,lat_meshgrid] = create_aerpaw_square_grid(grid_dist)
% grid_dist in meter
    scaler = 111139 ;
    grid_res = grid_dist/scaler; % in latitude or longitude
    
    mX_range = -78.7002: grid_res: -78.6961;
    mY_range = 35.7232: grid_res: 35.7303;

    [lon_meshgrid, lat_meshgrid] = meshgrid(mX_range, mY_range);
    lon_tar = lon_meshgrid(:) ;
    lat_tar = lat_meshgrid(:) ;
end

