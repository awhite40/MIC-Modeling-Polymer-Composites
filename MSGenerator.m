function [ DMS, Vf, f ] = MsGenerator( A,H,n,Iso1,Iso2,Iso3,filter)
% Microstructure Generating Funcion
%   This function take a 2-D or 3-D matrix - A, and the number of phases -
%   H, the relative grain size and the distribtion of the phases - n, and
%   filter size modifiers - r, Iso1, Iso2, and returns an eigen
%   microstructure DMS, the the dimensions of A with H phases.
%
%   The parameter n can be thought of as volume fraction distribution
%   factor. When n = 1 the volume fraction for all H phases is equal. As
%   n increases, the volume fraction for the first phase increases and the
%   volume fraction for all other phases decreases.
%
%   The parameter r changes the size of the filter radious. The parameters
%   Iso1 and Iso2 make the filter more ellipsoidal.

%% Input - Debugging

% x=101; y=101; z=101;
% D = rand(x,y,z);
% A = D;
% H = 2; % Number of Phases
% n = 1.5; % Control volume fraction of phases
% % r = 2; % radial Grain size
% Iso1 = 20; % Dimension 1
% Iso2 = 2; % Dimension 2
% Iso3 = 2; % Dimension 3

% P = MS(D,H,n,r,Iso1,Iso2)
% MS( A,H,n,r,Iso1,Iso2)
%% Matrix and Phase volume fraction check
Dim = size(A);
len = length(Dim); % number of dimension of A
x = Dim(1);
y = Dim(2);
z = Dim(3);
% if r >= min(Dim)/3
%     error('The value for the parameter r, must be an integer smaller than 1/3 the smallest dimenion of the input matrix, A.')
% end

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
Pval = [0,cumsum(P)];

[X,Y,Z] = meshgrid(-ceil(fDim(2)/2):ceil(fDim(2)/2),-ceil(fDim(1)/2):ceil(fDim(1)/2),-ceil(fDim(3)/2):ceil(fDim(3)/2));

f = sqrt(X.^2+Y.^2+Z.^2);

%% Filters
if filter == 1
    % Filter 1
    cf = zeros(Dim);
    % cf1(1,1:size(f,2),1) = f(1,:,1);
    cf(1:Iso2,1:Iso1,1:Iso3) = f(1:Iso2,1:Iso1,1:Iso3);
else if filter == 2
        % Filter 2
        cf = zeros(Dim);
        dimf = size(f);
        linInd = 1:prod(dimf);
        [a,b,c] = ind2sub(dimf,linInd');
        ind = find(a~=b | a~=c);
        ff = f(:);
        ff(ind)=0;
        tmp1 = reshape(ff,dimf);
        % cf2(1:dimf(1),1:dimf(2),1:dimf(3)) = tmp1;
        cf(1:dimf(1),dimf(2):-1:1,dimf(3):-1:1) = tmp1;
    else if filter == 3
            % Filter 3
            slope1 = 1;
            slope2 = 1;
            normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
            
            fake = zeros(Dim);
            x = linspace(1,x,1000)';  %# x values at a higher resolution
            y = slope1*(x);
            z = slope2*(y);
            index = sub2ind(size(fake),round(y),round(x),round(z));  %# Indices of the line
            fake(index) = 5;
            cf = normalize(fake);
        else if filter == 4
                % Filter 1
                cf1 = zeros(Dim);
                % cf1(1,1:size(f,2),1) = f(1,:,1);
                cf1(1:Iso2,1:Iso1,1:Iso3) = f(1:Iso2,1:Iso1,1:Iso3);
                
                % Filter 2
                cf2 = zeros(Dim);
                dimf = size(f);
                linInd = 1:prod(dimf);
                [a,b,c] = ind2sub(dimf,linInd');
                ind = find(a~=b | a~=c);
                ff = f(:);
                ff(ind)=0;
                tmp1 = reshape(ff,dimf);
                % cf2(1:dimf(1),1:dimf(2),1:dimf(3)) = tmp1;
                cf2(1:dimf(1),dimf(2):-1:1,dimf(3):-1:1) = tmp1;
                
                % Filter 3
                slope1 = 1;
                slope2 = 1;
                normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
                fake = zeros(Dim);
                x = linspace(1,x,1000)';  %# x values at a higher resolution
                y = slope1*(x);
                z = slope2*(y);
                index = sub2ind(size(fake),round(y),round(x),round(z));  %# Indices of the line
                fake(index) = 5;
                cf3 = normalize(fake);
                
                cf = cf1 + cf2 + cf3;
            else if filter == 5
                    cf = zeros(Dim);
                    % cf1(1,1:size(f,2),1) = f(1,:,1);
                    cf(1:Iso1,1:Iso2,1:Iso3) = f(1:Iso1,1:Iso2,1:Iso3);                    
                end
            end
        end
    end
end

    % cf = cf1+cf2+cf3;
    % keyboard
    
    %% Convolution
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
    
    %To DO: make sure the voxel is originally set as 1.
    for ll = 1:H; %
        temp = (DMS3(:,:,:) >= Pval(ll)) &  (Pval(ll+1) >= DMS3(:,:,:) );
        DMStemp = DMStemp + temp*(ll+1);
    end
    
    DMS = DMStemp - ones(Dim);
    
    %% Volume Fraction
    
    a = find(DMS(:)==1);
    Vf = size(a,1)/(prod(Dim));
