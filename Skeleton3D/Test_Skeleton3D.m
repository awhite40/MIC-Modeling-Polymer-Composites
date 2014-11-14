clear all;
close all;

load('segSample2P2.mat')
normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
I2 = normalize(I);
test = logical(I2(1:5,1:5,1:5));


skel = Skeleton3D(test);

figure();
col=[.7 .7 .8];
hiso = patch(isosurface(test,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(test,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(test,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)
