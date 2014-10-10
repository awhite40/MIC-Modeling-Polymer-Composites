function [ values ] = Thresh3D_whole( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[nRows, nColumns, nFrames] = size(X); 
d= repmat(double(0), [nRows, nColumns, nFrames]);
minc = min(min(min((X(:,:,:)))));
maxc = max(max(max((X(:,:,:)))));
inte = zeros([(maxc-minc)/100, 1]);%edit denominator for finer thresholding
values = int16(zeros([(maxc-minc)/10,1]));%edit denominator for finer thresholding
b = 1;
for value = minc:10:maxc
    %Scanning through threshold values
    x = 1;
    while x <= nRows
        y = 1;
        while y <=nColumns
            z=1;
            while z <= nFrames
                if X(x,y,z) < value % Values less than the threshold are set to zero
                    d(x,y,z)= 0;
                else
                    d(x,y,z)= 1; %Everything else is set to 1
                end
                z=z+1;
            end
            y=y+1;
        end
        x = x+1;
    end
    
inte(b) = sum(sum(sum(d)))/(nRows*nColumns*nFrames);
values(b)=value;
b = b+1;
end

plot(values,inte,'red');
end

