
close all
clear
load('_data/Cropped_fiber_matrix_400_400.mat')

normalize = @(A)(A-min(A(:)))./(max(A(:))-min(A(:)));

small = (Scropped(1:50,1:50,1:50));
Thresh = NaN(size(small));
imshow(small(:,:,1));
for i=1:size(small,1)
    for j=1:size(small,2)
        for k=1:size(small,3)
            if small(i,j,k)<6500
                Thresh(i,j,k) = 0;
            else
                Thresh(i,j,k) =1;
            end
            
        end
    end
end
figure
imshow(Thresh(:,:,1))
%%
%skel = Skeleton3D(small);
labels = bwlabeln(Thresh);
figure
pcolor(double(labels(:,:,1)));
%%
dy = NaN(1,max(labels(:)));
dx = NaN(1,max(labels(:)));
dz = NaN(1,max(labels(:)));
angle = NaN(size(labels));
j = 0;
[r2.(sprintf('l%d',j)),c2.(sprintf('l%d',j)),v2.(sprintf('l%d',j))] = ind2sub(size(labels),find(labels == j));
for i = 1:size(r2.(sprintf('l%d',j)))
            angle(r2.(sprintf('l%d',j))(i),c2.(sprintf('l%d',j))(i),v2.(sprintf('l%d',j))(i)) =NaN;
end
for j=1:max(labels(:))
    [r2.(sprintf('l%d',j)),c2.(sprintf('l%d',j)),v2.(sprintf('l%d',j))] = ind2sub(size(labels),find(labels == j));
    if (size(r2.(sprintf('l%d',j)),1)<5 && size(c2.(sprintf('l%d',j)),1)<5 && size(v2.(sprintf('l%d',j)),1) <5)
        for i = 1:size(r2.(sprintf('l%d',j)))
            labels(r2.(sprintf('l%d',j))(i),c2.(sprintf('l%d',j))(i),v2.(sprintf('l%d',j))(i)) =0;
            angle(r2.(sprintf('l%d',j))(i),c2.(sprintf('l%d',j))(i),v2.(sprintf('l%d',j))(i)) =NaN;
        end
    else
        dy(j) = max(r2.(sprintf('l%d',j)))-min(r2.(sprintf('l%d',j)));
        dx(j) = max(c2.(sprintf('l%d',j)))-min(c2.(sprintf('l%d',j)));
        dz(j) = max(v2.(sprintf('l%d',j)))-min(v2.(sprintf('l%d',j)));
        for i = 1:size(r2.(sprintf('l%d',j)))
            angle(r2.(sprintf('l%d',j))(i),c2.(sprintf('l%d',j))(i),v2.(sprintf('l%d',j))(i)) =atan(dz(j)/dy(j))*180/pi();
        end
        
    end
    
end
%%
figure

pcolor(angle(:,:,1));
colorbar

