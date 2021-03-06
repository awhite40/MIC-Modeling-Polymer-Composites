close all
[x,y,z] = size(I);
figure
ims = I(:,:,10);
pcolor( I(:,:,10) ); axis equal; axis tight; shading flat
%% If matrix is in 1s and 2s 
I2 = repmat(double(0), [x, y, z]);
i = 1;
while i< x
    j = 1;
    while j <y
        k = 1;
        while k<z
            if  I(i,j,k) == 1
                I2(i,j,k) = 0;
            else 
                I2(i,j,k) =1;
            end
            k =k+1;
        end
        j=j+1;
    end
    i = i+1;
end
figure
pcolor( I2(:,:,10) ); axis equal; axis tight; shading flat
%% Smoothing

sm = smooth3(I);
[x,y,z] = size(sm);
sm2 = repmat(double(0), [x, y, z]);
i = 1;
while i< x
    j = 1;
    while j <y
        k = 1;
        while k<z
            if  sm(i,j,k) >1.5
                sm2(i,j,k) = 1;
            else 
                sm2(i,j,k) =0;
            end
            k =k+1;
        end
        j=j+1;
    end
    i = i+1;
end
figure
pcolor( sm2(:,:,10) ); axis equal; axis tight; shading flat

%% Shifting
I_s= repmat(double(0), [x, y, z]);
I_s2= repmat(double(0), [x, y, z]);
I_s3= repmat(double(0), [x, y, z]);
i = 1;
while i< x-1
    I_s(i,:,:) = I(i+1,:,:);
    i = i +1;
end
I_s(x,:,:) = I(1,:,:);
i = 1;
while i< y-1
    I_s2(:,i,:) = I_s(:,i+1,:);
    i = i +1;
end
I_s2(:,y,:) = I_s(:,1,:);
while i< z-1
    I_s3(:,:,i) = I_s2(:,:,i+1);
    i = i +1;
end
I_s3(:,:,z) = I_s2(:,:,1);
[xs,ys,zs] = size(I_s);
figure
pcolor( I_s3(:,:,10) ); axis equal; axis tight; shading flat
shift = (I-I_s3).*I;
figure
pcolor( shift(:,:,10) ); axis equal; axis tight; shading flat
%% Labeling 
L = bwlabeln(I==2);

%%
L2 = bwlabeln(I2,6);
%%
L3 = bwlabeln(sm2,6);
L4 = bwlabeln(shift,6);
%L2 = bwlabeln(sm==2
maximum_L = max(L(:))
maximum_L2 = max(L2(:))
maximum_L3 = max(L3(:))
maximum_L4 = max(L4(:))
pcolor( I2(:,:,10 ) ); axis equal; axis tight; shading flat
figure
pcolor( shift(:,:,10) ); axis equal; axis tight; shading flat
%%


m1 = mode(L(L~=0)) 

%%
m2 = mode(L2(L2~=0))
m3 = mode(L3(L3~=0))
m4 = mode(L4(L4~=0))

%% Fibers that transverse the volume
i = 1;
through= zeros(maximum_L4,1);
j = 1;
while i <= maximum_L4/100
    % if a fiber can be found on both sides
    if any(any(L4(:,:,1) == i)) && any(any(L4(:,:,z) == i))
        through(j) = i;
        j=j+1;
    end
%     if any(any(L4(:,1,:) == i)) && any(any(L4(:,y,:) == i))
%         through(j) = i;
%         j=j+1;
%     end
%     if any(any(L4(1,:,:) == i)) && any(any(L4(x,:,:) == i))
%         through(j) = i;
%         j=j+1;
%     end
    i = i +1;
end
through
%% Segmented
figure
[x1,y1,z1] = ind2sub( size(L), find( L == 0));
plot3(x,y,z,'o')
axis equal
grid on
%% Segmented Adjusted
figure
[x2,y2,z2] = ind2sub( size(L2), find( L2 == m2));
plot3(x2,y2,z2,'o')
axis equal
grid on

%% Smoothed
figure
[x3,y3,z3] = ind2sub( size(L3), find( L3 == m3));
plot3(x3,y3,z3,'o')
axis equal
grid on

%% Shifted

xx = 1;
while xx <maximum_L4
    sz = size(find(L4 == xx));
    if sz(1)>90
        figure
        [x4,y4,z4] = ind2sub( size(L4), find( L4 == xx));
        plot3(x4,y4,z4,'.')
        axis equal
        grid on
    end
    
    xx = xx +1;
end
    
figure
[x4,y4,z4] = ind2sub( size(L4), find( L4 == m4));
plot3(x4,y4,z4,'o')
axis equal
grid on

%% Filtering
H = repmat(double(1), [3, 3, 10]);
I_F = imfilter(I2,H);
figure, imshow(I_F(:,:,10))

