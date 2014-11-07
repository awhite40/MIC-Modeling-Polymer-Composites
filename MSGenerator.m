clear
clc

%% Input - Debugging

x=11; y=11; z=11;
D = rand(x,y,z);

H = 2; % Number of Phases
% Control volume fraction of phases
fr = 0.15;
n = log(fr)/log(1/2);
r = 2; % radial Grain size
Iso1 = 3; % Change in Grain size (stretch in 1-Dir)
Iso2 = 5;

% P = MS(D,H,n,r,Iso1,Iso2)
%MS( A,H,n,r,Iso1,Iso2)
%% Matrix and Phase volume fraction check

Dim = size(D);
len = length(Dim); % number of dimension of A
if r >= min(Dim)/3
    error('The value for the parameter r, must be an integer smaller than 1/3 the smallest dimenion of the input matrix, A.')
end

if len == 2;
    DMS = rand(Dim(1),Dim(2),1);
    Dim = [Dim,1];
    f = zeros(ceil(2*Dim(1)/3),ceil(2*Dim(2)/3),1);%Will be used as a filter later, assumed to be 1/2 size of actual structure.
    fDim = [size(f),1];
    if Iso2 > 1
        warning('The parameter Iso2 is not used for 2-D datasets.')
    end
    
elseif len == 3;
    DMS = rand(Dim(1),Dim(2),Dim(3));
    f = zeros(ceil(2*Dim(1)/3),ceil(2*Dim(2)/3),ceil(2*Dim(3)/3));%Will be used as a filter later, assumed to be 1/2 size of actual structure.
    fDim = size(f);
    
else error('The input matrix A must must be either a 2-D or 3-d');
end

if H < 2;
    error('The number of phases, H, must be greater than 1');
end

P = zeros(1,H);%Percentage of volume fractions

for ii = 1:H-1
    P(ii+1)=(1/H)^(n^(sqrt(ii)));
end

SumP = cumsum(P);
Ptot = SumP(H);
P(1)=1-Ptot;
Pval = [0,cumsum(P)]
%% Filter and convulution

[X,Y,Z] = meshgrid(-ceil(fDim(2)/2):ceil(fDim(2)/2),-ceil(fDim(1)/2):ceil(fDim(1)/2),-ceil(fDim(3)/2):ceil(fDim(3)/2));

f = sqrt(X.^2+Y.^2+Z.^2);

%% Filter 1
cf = zeros(Dim);
cf(1,1:size(f,2),1) = f(1,:,1);

%% Filter 2
% cf = zeros(Dim);
% dimf = size(f);
% linInd = 1:prod(dimf);
% [a,b,c] = ind2sub(dimf,linInd');
% ind = find(a~=b | a~=c);
% ff = f(:);
% ff(ind)=0;
% tmp1 = reshape(ff,dimf);
% cf(1:dimf(1),1:dimf(2),1:dimf(3)) = tmp1;




%Convolution
DMS1 = ifftn(fftn(cf).*fftn(DMS));
DMS2 = ifftn(fftn(cf).*fftn(DMS1));
D = min(min(min(DMS2)));
E = ones(Dim)*D;
DMS2 = DMS2-E;
C = max(max(max(DMS2)));
DMS3 = DMS2/C;

%% Phases
%Assign different phases.             

DMStemp = zeros(Dim);

%To DO: make sure the voxel originally set equal to 1.
for ll = 1:H; %
    temp = (DMS3(:,:,:) >= Pval(ll)) &  (Pval(ll+1) >= DMS3(:,:,:) );
    DMStemp = DMStemp + temp*(ll+1)
end

DMS = DMStemp - ones(Dim)
[i,j,k] = size(DMS);
tempDMS = DMS-ones(Dim);
frac = sum(tempDMS(:))/(i*j*k)
