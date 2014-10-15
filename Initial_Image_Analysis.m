%% Load in Data

data.('name') = 'Cropped_fiber_matrix_S_400_400'
data.('local') = '_data/Cropped_fiber_matrix_400_400.mat'
data.('header') = '_data/Cropped_fiber_matrix_400_400.mat.json'

% Mat variable name called cropped
load( data.local);
cropped = double( cropped );


data.('voxels').('min') = min( cropped(:)) ;
data.('voxels').('max') = max( cropped(:)) ;
data.('voxels').('mean') =  mean( cropped(:)) ;
data.('voxels').('std') =  std( cropped(:)) ;
data.('voxels').('size') =  size( cropped) ;

s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
%% Visualization Tools
% Vol3d
% Sliceomatic

%% Visualize Faces
pcolor( squeeze( cropped(:,:,1) )); 
shading flat;
axis equal;
colorbar;
colormap gray;
saveas(gcf, './assets/S_topslice-xy.png' );
data.('images').('top')(1) = struct( 'src', '/assets/S_topslice-xy.png', 'description', 'Top XY Slice of Sample Image' );

pcolor( squeeze( cropped(:,1,:) )); 
shading flat;
axis equal;
colorbar;
colormap gray;
saveas(gcf, './assets/S_topslice-xz.png' );
data.('images').('top')(2) = struct( 'src', '/assets/S_topslice-xz.png', 'description', 'Top XZ Slice of Sample Image' );


pcolor( squeeze( cropped(1,:,:) )); 
shading flat;
axis equal;
colorbar;
colormap gray;
saveas(gcf, './assets/S_topslice-yz.png' );
data.('images').('top')(3) = struct( 'src', '/assets/S_topslice-yz.png', 'description', 'Top YZ Slice of Sample Image' );

s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Visualize Slices
pcolor( squeeze( cropped(:,:,round(size(cropped,3)./2)) )); 
shading flat
axis equal;
colorbar
colormap gray
saveas(gcf, './assets/midslice-xy.png' )
data.('images').('mid')(1) = struct( 'src', '/assets/midslice-xy.png', 'description', 'Middle XY Slice of Sample Image' )

pcolor( squeeze( cropped(:,round(size(cropped,2)./2),:) )); 
shading flat
axis equal;
colorbar
colormap gray
saveas(gcf, './assets/midslice-xz.png' )
data.('images').('mid')(2) = struct( 'src', '/assets/midslice-xz.png', 'description', 'Middle XZ Slice of Sample Image' )


pcolor( squeeze( cropped(round(size(cropped,1)./2),:,:) )); 
shading flat
axis equal;
colorbar
colormap gray
saveas(gcf, './assets/midslice-yz.png' )
data.('images').('mid')(3) = struct( 'src', '/assets/midslice-yz.png', 'description', 'Middle YZ Slice of Sample Image' )


s = savejson( [], data )
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Image Histogram
cropped = double(cropped);
[yy,xx] = hist( cropped(:),151);

[ax, h1,h2] = plotyy( xx, yy,... Plot distribution
    xx, [gradient(yy); gradient(gradient(yy))] ... Plot gradients of distribution
    );
figure(gcf)
% set(gca,'Yscale','log')
legend([h1;h2], 'Raw','First Derivative','Second Derivative');
grid( ax(1), 'on');
set( [h1;h2], 'LineWidth',3);
set( ax,'Fontsize',16);
ylabel(ax(1), 'Raw Data');
ylabel(ax(2), 'Gradients');
saveas( gcf,'./assets/raw-image-stats.png');
data.('images').('distribution') = struct( 'src', '/assets/raw-image-stats.png', 'description', 'This distribution of pixels values and their gradients.' )
s = savejson( [], data )
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
%% Peak fitting
figure
[yy,xx] = hist( cropped(:),151)
[ p,e ]= peakfit( [xx;yy], 0,0,2,1 );

saveas( gcf,'./assets/2_peak_fit_N.png')
data.('images').('Peak_fit2') = struct( 'src', '/assets/2_peak_fit_S.png', 'description', 'This is the peak fit using 2 peaks.' )

figure
[yy,xx] = hist( cropped(:),151)
[ p,e ]= peakfit( [xx;yy], 0,0,3,1 );

saveas( gcf,'./assets/3_peak_fi-N.png')
data.('images').('Peak_fit3') = struct( 'src', '/assets/3_peak_fit_S.png', 'description', 'This is the peak fit using 3 peaks.' )
s = savejson( [], data )
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
%% Using imadjust


%%
[yy,xx] = hist( cropped(:),151)
[ p,e ]= peakfit( [xx;yy], 0,0,3,1 );

figure(gcf);
grid on
hold( gca, 'on')
GRADE = [ gradient(yy); gradient(gradient(yy)) ];


%%%%%%%%%% NORMALIZE GRADIENTS FOR VISUALIZATION ONLY!
% The range of error values in the peak fit plot is orders of magnitude
% less than the gradients.  The normalization allows them to be viewed
% simultaneously.  We are not relying on the magnitudes of the gradients in
% this plot.

dGradient = diff([min(GRADE(:)),max(GRADE(:))])
dScale = diff( ylim( gca ) );
GRADEPLOT = GRADE .* dScale./dGradient;

%%%%%%%%%%
% h = plot( gca, xx, GRADEPLOT, ... Plot gradients of distribution
%     'LineWidth',3)
hold off

% legend( h, 'First Gradient', 'Second Gradient' )
set( gcf, 'Position', get(0,'ScreenSize')./[ 1 1 2 1] )
saveas( gcf,'./assets/peak-fit-test.png')
data.('images').('peaks') = struct( 'src', '/assets/peak-fit-test.png', 'description', 'Peak fitting of the histogram generated from a 3-D volume.' )



%% 
return
bnds = 6000
for ll = [100 + [-1:1]; [1:3]];
subplot(1,3,ll(2))
pcolor( double( cropped(:,:,ll(1)) > bnds(1) & cropped(:,:,ll(1)) > bnds(1) ) .* cropped(:,:,ll(1)) );
shading flat
axis equal
colorbar
end
figure(gcf)