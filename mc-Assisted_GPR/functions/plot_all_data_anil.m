function plot_all_data_anil(mX_all,mY_all,mZ_all,mPower_all)

lat_BS = 35.727451 ; % [degree]
lon_BS = -78.695974 ; % [degree]
offset = [-78.698 35.727] ;
scaler = 6371000*2*pi/360; %111139 ;
scalerX = 6371000*2*pi/360*cosd(35.727451);
figure
hold on
scatter3( (mX_all - offset(1))*scalerX,(mY_all - offset(2))*scaler,mZ_all,10,mPower_all)
%viscircles([0,0],[40],'Color','k')
colormap(jet);
colorbar;
caxislim=[-160 -50];
%caxis(caxislim)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Height [m]')
%zlim([0 120]);
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'SF (dB)';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
fontsize(gcf,24,"points")
grid on;
scatter3( (lon_BS - offset(1))*scaler,(lat_BS - offset(2))*scaler,12,50,"red",'filled','marker','^')
end