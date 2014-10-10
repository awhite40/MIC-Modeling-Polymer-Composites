function slices_view( original, threshold, slice,slice2 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
subplot(2,2,1)
imshow(original(:,:,slice));
subplot(2,2,2)
imshow (threshold(:,:,slice));

subplot(2,2,3)
imshow(original(:,:,slice2));
subplot(2,2,4)
imshow (threshold(:,:,slice2));
end

