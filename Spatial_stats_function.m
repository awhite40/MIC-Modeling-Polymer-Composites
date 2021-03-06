function [ data ] = Spatial_stats_function(matr)
% Calculates spacial statistics on input matrix which is binary and 3
% dimensional.  Also calculates 2D statistics on 3 orthonginal slices
%   Detailed explanation goes here

normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
cropped = double( (matr) );
data.phase =normalize(cropped);
data.stats = SpatialStatsFFT(data.phase, [],'periodic',true);
data.cstats = SpatialStatsFFT(data.phase==1,data.phase==0,'periodic',true);

% data.stats_xy_slice = SpatialStatsFFT( data.phase(:,:,slice),[],'periodic',true);
% data.stats_xz_slice = SpatialStatsFFT( squeeze(data.phase(:,slice,:)),[],'periodic',true);
% data.stats_yz_slice = SpatialStatsFFT( squeeze(data.phase(slice,:,:)),[],'periodic',true);

end

