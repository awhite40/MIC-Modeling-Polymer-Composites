
normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) )
%% Spacial Stats on Original sample 2

data.('name') = 'Spatial_Stats_Original2';
data.('local') = '_data/Spatial_Stats_Original2.mat';
data.('header') = '_data/Spatial_Stats_Original2.mat.json';

data.('source') = '_data/Ncropped.mat';

load( data.source);
cropped = double( Ncropped );
org2.phase =(cropped);
org2.stats_xy150 = SpatialStatsFFT( org2.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_xy.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_xy.png', 'description', 'Spatial statistics on original image N in X-Y direction.' );
org2.stats_xz150 = SpatialStatsFFT( squeeze(org2.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_xz.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_xz.png', 'description', 'Spatial statistics on original image N in X-Z direction.' );
org2.stats_yz150 = SpatialStatsFFT( squeeze(org2.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_yz.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_yz.png', 'description', 'Spatial statistics on original image N in Y-Z direction.' );


s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on Original sample 1

data.('name') = 'Spatial_Stats_Original';
data.('local') = '_data/Spatial_Stats_Original.mat';
data.('header') = '_data/Spatial_Stats_Original.mat.json';

data.('source') = '_data/Cropped_fiber_matrix_400_400.mat';

load( data.source);
cropped = double( Scropped );
org.phase = (cropped);
org.stats_xy150 = SpatialStatsFFT( org.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_xy.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_xy.png', 'description', 'Spatial statistics on original image S in X-Y direction.' );
org.stats_xz150 = SpatialStatsFFT( squeeze(org.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_xz.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_xz.png', 'description', 'Spatial statistics on original image in X-Z direction.' );
org.stats_yz150 = SpatialStatsFFT( squeeze(org.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_yz.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_yz.png', 'description', 'Spatial statistics on original image in Y-Z direction.' );

s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on segmented sample 2

data.('name') = 'Spatial_Stats_Segmented2';
data.('local') = '_data/Spatial_Stats_Segmented2.mat';
data.('header') = '_data/Spatial_Stats_Segmented2.mat.json';

data.('source') = '_data/SegSample2P2.mat';

load( data.source);
cropped = double( I);
seg2.phase = (cropped);
seg2.stats_xy150 = SpatialStatsFFT( seg2.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_xy.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_xy.png', 'description', 'Spatial statistics on segmented image N in X-Y direction.' );
seg2.stats_xz150 = SpatialStatsFFT( squeeze(seg2.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_xz.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_xz.png', 'description', 'Spatial statistics on segmented image N in X-Z direction.' );
seg2.stats_yz150 = SpatialStatsFFT( squeeze(seg2.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_yz.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_yz.png', 'description', 'Spatial statistics on segmented image N in Y-Z direction.' );


s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on segmented sample 1

data.('name') = 'Spatial_Stats_Segmented';
data.('local') = '_data/Spatial_Stats_Segmented.mat';
data.('header') = '_data/Spatial_Stats_Segmented.mat.json';

data.('source') = '_data/SegSample1P2.mat';

load( data.source);
cropped = double( I);
seg.phase = (cropped);
seg.stats_xy150 = SpatialStatsFFT( seg.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_xy.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_xy.png', 'description', 'Spatial statistics on segmented image S in X-Y direction.' );
seg.stats_xz150 = SpatialStatsFFT( squeeze(seg.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_xz.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_xz.png', 'description', 'Spatial statistics on segmented image S in X-Z direction.' );
seg.stats_yz150 = SpatialStatsFFT( squeeze(seg.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_yz.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_yz.png', 'description', 'Spatial statistics on segmented image S in Y-Z direction.' );



s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
