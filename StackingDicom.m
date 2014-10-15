function [X] = StackingDicom(filename,short,nFrames,start)
%Read in file information
info = dicominfo(filename);
%nRows = info.Rows;
%nColumns = info.Columns;
%nPlanes = info.SamplesPerPixel;
xsize = 300;
ysize = 300;
X = repmat(int16(0), [xsize+1, ysize+1, nFrames]); %create empty array
%Stacking slices
n = 1;
for p=start:(nFrames-1+start)
    filename = sprintf('%s%05d.DCM', short,p);
    
    J= dicomread(filename);%add files to array
    % J is your original matrix obtained from a DICOM Imag
    J1 = imrotate(J,-4,'bilinear','crop');
    % Cropping
    J11=imcrop(J1, [400 400 ysize xsize]);
    %figure
    %imshow(J11)
%     [xs,ys]=size(J11)
%     x= 1;
%     while x<=300;
%         y = 1;
%         while y<=300;
%             bin = 0;
%             while bin<255;
%                 if (J11(x,y) > bin*248 && J11(x,y)< bin*248 + 248);
%                     J11 = int8(bin);
%                 end 
%                 bin = bin+1;
%             end
%             y = y+1;
%         end
%         x = x+1;
%     end
    X(:,:,n) = J11;
    n= n+1;
end
%figure
%imshow(X(:,:,10));
end



