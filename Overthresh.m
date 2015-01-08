tic
close all
load('_data/Ncropped');


normalize = @(A)( double(A)-double(min(A(:))) ) ./ (double( max(A(:))) - double(min(A(:))) );
small = Ncropped(1:150,1:150,1:150);
imshow(small(:,:,1));
norms = normalize(small);
%%
for i=1:size(small,1)
    for j=1:size(small,2)
        for k=1:size(small,3)
            
            if small(i,j,k)<6500
                Thresh6500(i,j,k) = 0;
            else
                Thresh6500(i,j,k) =1;
            end
            if small(i,j,k)<7000
                Thresh7000(i,j,k) = 0;
            else
                Thresh7000(i,j,k) =1;
            end
            if small(i,j,k)<8500
                Thresh8500(i,j,k) = 0;
            else
                Thresh8500(i,j,k) =1;
            end
            if small(i,j,k)<8000
                Thresh8000(i,j,k) = 0;
            else
                Thresh8000(i,j,k) =1;
            end
        end
    end
end

%%

figure
imshow(Thresh6500(:,:,1))
figure

imshow(Thresh7000(:,:,1))

figure
imshow(Thresh8000(:,:,1))
figure
imshow(Thresh8500(:,:,1))
% diff = Thresh6000-Thresh9000;
% figure 
% imshow(diff(:,:,1))
toc