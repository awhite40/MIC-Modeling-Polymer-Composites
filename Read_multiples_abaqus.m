clear
clc
close all
s = 100;
%Savg = zeros(3,3,s);
%Eavg = zeros(3,3,s);
% Cavg = zeros(3,3,3,3,s);
% Cs_avg = zeros(3,3,3,3,s);
CavgY = NaN(3,3,3,3,s);
Cs_avgY = NaN(3,3,3,3,s);
% CavgD = zeros(3,3,3,3,s);
% Cs_avgD = zeros(3,3,3,3,s);
% C11D = zeros(s,1);
% C11X = zeros(s,1);
%%
% for p=101:288
%     filename = sprintf('Abaqus_files/ElongFib%d.dat',p);
%     [Cavg(:,:,:,:,p),Cs_avg(:,:,:,:,p)] = My_Abaqus_file_read( filename );
% end

%%
% 
% C11X(1:s) = squeeze(Cavg(1,1,1,1,1:s));
% 
% 
% indX = find(isnan(C11X));

%%
tic
for p=101:200

    filename = sprintf('ElongYFib%d.dat',p);
    if exist(filename, 'file') == 2
        [CavgY(:,:,:,:,p),Cs_avgY(:,:,:,:,p)] = My_Abaqus_file_read( filename );
    end
    
end
toc

%%
C11Y = zeros(364,1);
C11Y(1:364) = squeeze(CavgY(1,1,1,1,1:364));

% %%
% for p=4:5
%     filename = sprintf('Abaqus_files/DiagFib%d.dat',p);
%     [CavgD(:,:,:,:,p),Cs_avgD(:,:,:,:,p)] = My_Abaqus_file_read( filename );
% end
% 
% 
% C11D(1:s) = squeeze(CavgD(1,1,1,1,1:s));
% 
% indD = find(isnan(C11D));

%%
filename = 'Ninp1.dat';
[CavgN,Cs_avgN] = My_Abaqus_file_read(filename,21);


%%
filename = 'Sinp1.dat';
[CavgS,Cs_avgS] = My_Abaqus_file_read(filename,21);

