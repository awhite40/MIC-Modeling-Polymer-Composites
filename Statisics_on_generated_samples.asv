
normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) )
%% Spacial Stats on Elongated Sample

data.('name') = 'Generated microstructure with elongated fibers';
data.('local') = '_data/ElongFib.mat';
data.('header') = '_data/ElongFib.mat.json';

data.('source') = '_data/ElongFib.mat';

load( data.source);
cropped = double( ElongFib );
Elong.phase =normalize(cropped);
Elong.stats = SpatialStatsFFT(Elong.phase,[],'periodic',true);
Elong.cstats = SpatialStatsFFT(Elong.phase==1,Elong.phase==0,'periodic',true);
Elong.stats_xy10 = SpatialStatsFFT( Elong.phase(:,:,10),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Elong_xy.jpg' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Elong_xy.jpg', 'description', 'Spatial statistics on generated microstructure with elongated fibers in the XY plane.' );
Elong.stats_xz10 = SpatialStatsFFT( squeeze(Elong.phase(:,10,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Elong_xz.jpg' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Elong_xy.jpg', 'description', 'Spatial statistics on generated microstructure with elongated fibers in the XZ plane.' );
Elong.stats_yz10 = SpatialStatsFFT( squeeze(Elong.phase(10,:,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Elong_yz.jpg' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Elong_xy.jpg', 'description', 'Spatial statistics on generated microstructure with elongated fibers in the YZ plane.'  );

s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);
%% Statistics on Diagonal fibers
data.('name') = 'Generated microstructure with diagonal fibers';
data.('local') = '_data/DiagFib.mat';
data.('header') = '_data/DiagFib.mat.json';

data.('source') = '_data/DiagFib.mat';

load( data.source);
cropped = double( DiagFib );
Diag.phase =normalize(cropped);
Diag.stats = SpatialStatsFFT(Diag.phase,[],'periodic',true);
Diag.cstats = SpatialStatsFFT(Diag.phase==1,Diag.phase==0,'periodic',true);
Diag.stats_xy10 = SpatialStatsFFT( Diag.phase(:,:,10),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Diag_xy.jpg' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Diag_xy.jpg', 'description', 'Spatial statistics on Generated microstructures with diagonal fibers X-Y direction.' );
Diag.stats_xz10 = SpatialStatsFFT( squeeze(Diag.phase(:,10,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Diag_xz.jpg' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Diag_xz.jpg', 'description', 'Spatial statistics on Generated microstructures with diagonal fibers in X-Z direction.' );
Diag.stats_yz10 = SpatialStatsFFT( squeeze(Diag.phase(10,:,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Diag_yz.jpg' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Diag_yz.jpg', 'description', 'Spatial statistics on Generated microstructures with diagonal fibers in Y-Z direction.' );
s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);


%% Statistics on Random fibers

data.('name') = 'Generated microstructure with random fibers';
data.('local') = '_data/RandFib.mat';
data.('header') = '_data/RandFib.mat.json';

data.('source') = '_data/RandFib.mat';

load( data.source);
cropped = double( RandFib );
Rand.phase =normalize(cropped);
Rand.stats = SpatialStatsFFT(Rand.phase,[],'periodic',true);
Rand.cstats = SpatialStatsFFT(Rand.phase==1,Rand.phase==0,'periodic',true);
Rand.stats_xy10 = SpatialStatsFFT( Rand.phase(:,:,10),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Rand_xy.jpg' );
data.('images').('stats')(1) = struct( 'src', '/assets/Spacial_Stats_Rand_xy.jpg', 'description', 'Spatial statistics on Generated microstructures with Random fibers X-Y direction.' );
Rand.stats_xz10 = SpatialStatsFFT( squeeze(Rand.phase(:,10,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Rand_xz.jpg' );
data.('images').('stats')(2) = struct( 'src', '/assets/Spacial_Stats_Rand_xz.jpg', 'description', 'Spatial statistics on Generated microstructures with random fibers in X-Z direction.' );
Rand.stats_yz10 = SpatialStatsFFT( squeeze(Rand.phase(10,:,:)),[],'periodic',true);
saveas(gcf, './assets/Spacial_Stats_Rand_yz.jpg' );
data.('images').('stats')(3) = struct( 'src', '/assets/Spacial_Stats_Rand_yz.jpg', 'description', 'Spatial statistics on Generated microstructures with random fibers in Y-Z direction.' );
s = savejson( [], data );
fo = fopen( data.header, 'w'); fwrite( fo, s ); fclose(fo);

%%
close all
hsp = surf(linspace(0,21,21),linspace(0,21,21),...
      zeros(21) +7.5);
rotate(hsp,[-1,0,0],45)
xd = get(hsp,'XData');
yd = get(hsp,'YData');
zd = get(hsp,'ZData');

slice(fftshift(Rand.stats),xd,yd,zd)
hold on
slice(fftshift(Rand.stats),11,11,1)
colorbar
hold off

figure

slice(fftshift(Elong.stats),11,11,1)
colorbar
