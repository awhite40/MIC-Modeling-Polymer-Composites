load('SegSample2P2.mat')

normalize = @(A)(A-min(A(:)))./(max(A(:))-min(A(:)));

small = normalize(I(1:50,1:50,1:50));
%skel = Skeleton3D(small);
labels = bwlabeln(small);

[r,c,v] = ind2sub(size(labels),find(labels == 1));

[r2,c2,v2] = ind2sub(size(labels),find(labels == 2));

for i = 1:size(r)
    
    if labels(r(i),c(i),v(i)) ~=1
        x = 2;
    end
end

for i = 1:size(r2)
    
    if labels(r2(i),c2(i),v2(i)) ~=2
        x = 2;
    end
end

dy2 = max(r2)-min(r2)
dx2 = max(c2)-min(c2)
dz2 = max(v2)-min(v2)


dx = max(r)-min(r)
dy = max(c)-min(c)
dz = max(v)-min(v)

