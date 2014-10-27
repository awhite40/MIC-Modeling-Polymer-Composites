
normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) )
%% Spacial Stats on Original sample 2

data.('name') = 'Spatial_Stats_Original2_normalized';
data.('local') = '_data/Spatial_Stats_Original2_normalized.mat';
data.('header') = '_data/Spatial_Stats_Original2_normalized.mat.json';

data.('source') = '_data/Ncropped.mat';

load( data.source);
cropped = double( Ncropped );
org2.phase = normalize(cropped);
org2.stats_xy150 = SpatialStatsFFT( org2.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_xy_normalized.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_xy_normalized.png', 'description', 'Spatial statistics on original image N in X-Y direction.' );
org2.stats_xz150 = SpatialStatsFFT( squeeze(org2.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_xz_normalized.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_xz_normalized.png', 'description', 'Spatial statistics on original image N in X-Z direction.' );
org2.stats_yz150 = SpatialStatsFFT( squeeze(org2.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample2_yz_normalized.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Original_sample2_yz_normalized.png', 'description', 'Spatial statistics on original image N in Y-Z direction.' );


s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on Original sample 1

data.('name') = 'Spatial_Stats_Original_normalized';
data.('local') = '_data/Spatial_Stats_Original_normalized.mat';
data.('header') = '_data/Spatial_Stats_Original_normalized.mat.json';

data.('source') = '_data/Cropped_fiber_matrix_400_400.mat';

load( data.source);
cropped = double( Scropped );
org.phase = normalize(cropped);
org.stats_xy150 = SpatialStatsFFT( org.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_xy_normalized.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_xy_normalized.png', 'description', 'Spatial statistics on original image S in X-Y direction.' );
org.stats_xz150 = SpatialStatsFFT( squeeze(org.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_xz_normalized.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_xz_normalized.png', 'description', 'Spatial statistics on original image in X-Z direction.' );
org.stats_yz150 = SpatialStatsFFT( squeeze(org.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Original_sample1_yz_normalized.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Original_sample1_yz_normalized.png', 'description', 'Spatial statistics on original image in Y-Z direction.' );

s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on segmented sample 2

data.('name') = 'Spatial_Stats_Segmented2_normalized';
data.('local') = '_data/Spatial_Stats_Segmented2_normalized.mat';
data.('header') = '_data/Spatial_Stats_Segmented2_normalized.mat.json';

data.('source') = '_data/SegSample2P2.mat';

load( data.source);
cropped = double( I);
seg2.phase = normalize(cropped);
seg2.stats_xy150 = SpatialStatsFFT( seg2.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_xy_normalized.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_xy_normalized.png', 'description', 'Spatial statistics on segmented image N in X-Y direction.' );
seg2.stats_xz150 = SpatialStatsFFT( squeeze(seg2.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_xz_normalized.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_xz_normalized.png', 'description', 'Spatial statistics on segmented image N in X-Z direction.' );
seg2.stats_yz150 = SpatialStatsFFT( squeeze(seg2.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample2_yz_normalized.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample2_yz_normalized.png', 'description', 'Spatial statistics on segmented image N in Y-Z direction.' );


s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%% Spacial Stats on segmented sample 1

data.('name') = 'Spatial_Stats_Segmented_normalized';
data.('local') = '_data/Spatial_Stats_Segmented_normalized.mat';
data.('header') = '_data/Spatial_Stats_Segmented_normalized.mat.json';

data.('source') = '_data/SegSample1P2.mat';

load( data.source);
cropped = double( I);
seg.phase = normalize(cropped);
seg.stats_xy150 = SpatialStatsFFT( seg.phase(:,:,150));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_xy_normalized.png' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_xy_normalized.png', 'description', 'Spatial statistics on segmented image S in X-Y direction.' );
seg.stats_xz150 = SpatialStatsFFT( squeeze(seg.phase(:,150,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_xz_normalized.png' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_xz_normalized.png', 'description', 'Spatial statistics on segmented image S in X-Z direction.' );
seg.stats_yz150 = SpatialStatsFFT( squeeze(seg.phase(150,:,:)));
saveas(gcf, './assets/Spacial_Stats_Segmented_sample1_yz_normalized.png' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Segmented_sample1_yz_normalized.png', 'description', 'Spatial statistics on segmented image S in Y-Z direction.' );



s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
