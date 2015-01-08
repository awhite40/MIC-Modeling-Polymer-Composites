%% Load in Dependencies
%
% * Gaussian Likelihood segmentation containing an external peak fitting
% code
% * Some image adjustment tools
% * REQUIRES, BUT NOT INCLUDED:
% <http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave JSONlab>


%% Get Gaussian Likelihood Segmentation Codes

loc = 'segmentprob'; 
if ~exist( loc, 'dir' );
    try
        system('sudo git clone https://gist.github.com/tonyfast/de3d7b841cea38dc40fa.git');
    catch
        
        apicall = loadjson( urlread('https://api.github.com/gists/de3d7b841cea38dc40fa') );
        
        if ~exist( loc,'dir');  mkdir( loc ); end
        
        for ufilename = fieldnames(apicall.files)'
            if strcmp( apicall.files.(ufilename{1}).filename( end - [1:-1:0] ), '.m' )
                tofile = fullfile( loc, apicall.files.(ufilename{1}).filename );
                urlwrite( apicall.files.(ufilename{1}).raw_url, tofile );
            end
        end
    end
end
addpath( loc )

%% Get Image adjustment codes
%
% Steal these from segmentprob.m
% Adjust and normalize image data

normalize = @(A)( double(A)-double(min(A(:))) ) ./ (double( max(A(:))) - double(min(A(:))) );
adjust = @(A)reshape(  ... back to original shape
                imadjust( ... adjust image
                        reshape( ... flatten to 2-D image
                                normalize(A), ... normalize from zero to one
                                size(A,1), numel(A)./size(A,1))), ... Reshape to 2-D array
                                size(A) ... Reshape back to original size
                                );


%% Load in a sample
% Load in the <http://awhite40.github.io/MIC-Modeling-Polymer-Composites/2014/09/29/Sample-Analytics.html Original Polymer>
% composite sepcimen for testing.

load( '_data/Cropped_fiber_matrix_400_400.mat' )
cropped = double( Scropped );

%% Segment the Image
%
% The function takes in two numeric parameters.
% 
% * nbins  - Number of bins in the histogram. 
% * npeaks  - Number of peaks to find in the histogram. 
%
% Ideas to keep in mind:
%
% * _too many bins_ - the histogram is too noisy to fit
% * _too few bins_ - distinct peaks will not be resolved
%
% * _too many peaks_ - over fit the curve
% * Add more peaks than phases in case there are transition regions that
%   can be found.

A = adjust( normalize( cropped ) ) -  normalize( cropped ) ;
out = segmentprob( adjust( normalize(cropped) ), 151,3);

%%
%
% ABOVE: The Distribution of Voxel intensities in the image and the peak fit.

%% Plot the Segmentation
ll = [ 51  151 251 ];

for ii = [ ll ; ...
           1 : 3 ]
    ax(1) = subplot(3,3, (ii(2)-1) * 3 + 1);
    pcolor( cropped(:,:,ii(1) ) ); axis equal; axis tight; shading flat
    ylabel( sprintf('Layer %i', ii(1)), 'Fontsize', 16)
    title('Original')
    ax(2) = subplot(3,3, (ii(2)-1) * 3 + 2);
    pcolor( out.out.phase(:,:,ii(1) ) ); axis equal; axis tight; shading flat
    title('Phase')
    ax(3) = subplot(3,3, (ii(2)-1) * 3 + 3);
    pcolor( out.out.ci(:,:,ii(1) ) ); axis equal; axis tight; shading flat; colorbar
    title('Relative Confidence of the phase')
    linkaxes( ax )
end
figure(gcf)