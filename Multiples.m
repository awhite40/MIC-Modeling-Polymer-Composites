
clear
clc
close all
load('DiagFib.mat');
tic
for i = 1:size(DiagFib,4)
    
    data(:,i) = Spatial_stats_function(DiagFib(:,:,:,i));
    
end
toc


%%
s = 364;
Savg = zeros(3,3,s);
Eavg = zeros(3,3,s);
Cavg = zeros(3,3,3,3,s);
for p=1:s
    
    filename = sprintf('ElongFib%d.dat',p);
    [Savg(:,:,p),Eavg(:,:,p),Cavg(:,:,:,:,p)] = My_Abaqus_file_read( filename,21 );
end

%%
C11 = zeros(100,1);
C11(1:100) = squeeze(Cavg(1,1,1,1,1:100));
