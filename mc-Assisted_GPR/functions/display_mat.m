function display_mat(mX_grid,mY_grid,mPower2D)
    figure;
    %surf(mX_grid,mY_grid,mPower2D);
    %mesh(mX_grid,mY_grid,mPower2D,'FaceColor','flat')
    %imagesc(mX_grid, mY_grid, mPower2D);
    h = imagesc(mX_grid, mY_grid, mPower2D);
    set(gca, 'YDir', 'normal');
    colorbar;
    
    % Create an alpha map: 1 for valid data, 0 for NaN
    alpha_map = ~isnan(mPower2D);
    
    % Apply the transparency
    set(h, 'AlphaData', alpha_map);
    colormap(jet);
    colorbar;
    caxis([-15 15])
    hcb=colorbar;
    %colorTitleHandle = get(hcb,'Title');
    %titleString = 'Antenna gain [dBi]';
    %set(colorTitleHandle ,'String',titleString);
    %set(gcf,'color','w')
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('None');
    %set(gcf,'Color','w');
    %axis([-177 180 -90 90 min(new_patt(:)) max(new_patt(:)) ]);
    %view([0,0,100])
end