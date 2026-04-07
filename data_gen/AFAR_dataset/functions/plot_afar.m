function plot_afar(lats, lons, power_all)
    location =2;
    scaler = 111139 ;
    if location ==1
    origin_y=35.72806709;
    origin_x=-78.69730398;
    end
    
    if location==2
    origin_y=35.72911779;
    origin_x=-78.69918128;
    end
    
    if location==3
    origin_y=35.72985129;
    origin_x=-78.69711002;
    end

    cmap = jet;
    x1=-78.69953825817279;
    y1= 35.72688213193035;
    x2=-78.69621514941473;
    y2=35.72931030026633;

    figure
    hold on
    grid on
    box on
    %scatter(scaler*(lons-x1),scaler*(lats-y1),4,active_shadow);
    scatter(scaler*(lons-x1),scaler*(lats-y1),4,power_all);
    %scatter3(lons,lats,hs,8,floor((1:length(lats))/dur_len_step)) % 
    scatter(scaler*(origin_x-x1),scaler*(origin_y-y1),10,0,'MarkerEdgeColor', 'red')
    scatter(scaler*(origin_x-x1),scaler*(origin_y-y1),100,0, 'MarkerEdgeColor', 'red')
    colormap("jet");
    colorbar;
    zlabel('Altitude')
    xlabel('X [m]')
    ylabel('Y [m]')
    xlim([-10, 400])
    ylim([-10, 350])
    clim([-59,29])
    hcb=colorbar;
    colorTitleHandle = get(hcb,'Title');
    %titleString = 'Shadow (dB)';
    titleString = 'RSRP (dBm)';
    set(colorTitleHandle ,'String',titleString);
    set(gca, 'FontSize', 18); % Change font size for the current axis
    % x1=-78.69953825817279;
    % y1= 35.72688213193035;
    % x2=-78.69621514941473;
    % y2=35.72931030026633;
    %rectangle('Position',[x1 y1 (x2-x1) (y2-y1)])
    %title("exp "+num2str(exp_no)+" loc "+num2str(location))
end