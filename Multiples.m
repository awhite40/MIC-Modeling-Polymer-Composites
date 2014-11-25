
clear
clc
close all
load('RandFib.mat');
tic
for i = 1:size(RandFib,4)
    
    data(:,i) = Spatial_stats_function(RandFib(:,:,:,i));
    
end
toc



