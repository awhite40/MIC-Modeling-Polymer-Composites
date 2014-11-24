
clear
clc
close all
load('DiagFib.mat');
tic
for i = 1:size(DiagFib,4)
    
    data(:,i) = Spatial_stats_function(DiagFib(:,:,:,i));
    
end
toc




% tmp = fft(a,[],1);
% tmp = fft(tmp, [],2);
% tmp = fft(tmp, [],3);
% for ii = 1:600
%     fftn(a(:,:,:,ii));
% end



% conj(tmp).*tmp
