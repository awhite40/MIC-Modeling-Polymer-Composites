
clear
clc
close all
load('DiagFib.mat');
tic
for i = 1:size(DiagFib,4)
    
    data(:,i) = Spatial_stats_function(DiagFib(:,:,:,i));
    
end
toc



