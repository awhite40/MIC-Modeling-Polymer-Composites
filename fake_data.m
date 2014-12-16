close all
clear
clc

f =4 ;
slope =sqrt(2);

% Creating a fake microstructure in 2D
if f==1
    fake = zeros(30,30);

    fake(20:22,1:end) = 1;
    fake(1:20,20:22) = 1;
    fake(23:25,18:24) =1;
    fake(18:20,18:24) = 1;
    imshow(fake);
    colorbar;

% Fake microstructure 2
elseif f==2
    normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );

    fake = zeros(300,300,'uint8');
    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 30+sqrt(2)*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255; 
    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 31+sqrt(2)*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;%# Set the values to the max value of 255 for uint8


    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 29+sqrt(2)*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    y = 28+sqrt(2)*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    y = 31+sqrt(2)*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    fake = double(normalize(fake));
    imshow(fake);
    colorbar;

   

% Fake microstructure 3
elseif f==3
    normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
    fake = zeros(300,300,'uint8');
    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 30+slope*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255; 
    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 31+slope*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;%# Set the values to the max value of 255 for uint8


    x = linspace(30,180,10000);  %# x values at a higher resolution
    y = 29+slope*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    y = 28+slope*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    y = 31+slope*x;
    index = sub2ind(size(fake),round(y),round(x));  %# Indices of the line
    fake(index) = 255;
    fake = double(normalize(fake));

    fake(100:102,1:end) = 1;
    fake(1:100,50:52) = 1;
    fake(103:105,48:54) =1;
    fake(98:100,48:54) = 1;

    imshow(fake);
    colorbar;

     
elseif f==4
    normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
    load _data/SegSample2P2.mat
    fake = normalize(I(1:50,1:50,1:50));
    imshow(fake(:,:,10));
end




%% Isolating straight sections
% Run smaller window over image
[x, y] = size(fake);
new = NaN(x,y);
wind = 30; % Window size
i =wind;%108;
n= 1;
while i <= x%118
    j=wind;%90;
    while j<=y%100
        % Determine the white region in the window
        region = fake(i-wind+1:i,j-wind+1:j);
%         figure
%         imshow(region);
%         colorbar;
        % Determine the orientation and centroid of the white section in the window
        % assuming the possibility of unconnected sections in window
        % misses the edges without some sort of periodicity some extra
        % pixels found
        bw = bwlabel(region);
        lb = 1;
        while lb <= max(bw(:))
            reg = bw==lb;
            [row, col] = find(reg);
            row = row + i-wind;
            col = col + j-wind;
            cent(1) = sum(row)/sum(sum(reg));
            cent(2) = sum(col)/sum(sum(reg));
%             [V,D] = eigs(double(reg),3);
%             if V(:,1)~= V(:,2) & V(:,1)~= -V(:,2) 
%                 or = V(1:2,1);
%                 g=1;
%                 while g<=size(row)
%                 new2(row(g),col(g)) = norm(V(:,1));
%                 g=g+1;
%                 end
%             end
            [row2, col2] = find(reg(:,floor(cent(2))-j+wind));
            dX = max(row2) - min(row2);
            X = reg(:,floor(cent(2)-j +wind));
            sdx = std(X);
            sdy = std(reg(floor(cent(1)-i + wind),:));
            [row2, col2] = find(reg(floor(cent(1))- i+wind,:));
            dY = max(col2) - min(col2);
            rem.vec(n,:) = [sdx, sdy]/(sqrt(sdx^2+sdy^2));
            rem.centroid(n,:) = cent;
            rem.theta(n,:) = atan(rem.vec(n,2)/rem.vec(n,1));
            e = sqrt(1-((min(sdx,sdy))^2/(max(sdx,sdy)^2)));
            % if has an orientation add it to the seperated matrix
            g=1;
            while g<=size(row)
                new(row(g),col(g)) = 180*atan(rem.vec(n,2)/rem.vec(n,1))/pi();
                new2(row(g),col(g)) = e;
                g=g+1;
            end
            lb =lb+1;
            n= n+1;
        end
        j=j+1;
    end
    i=i+1;
end
%%
figure
pcolor((new));
colorbar
new3 = new2.*255;
figure
pcolor((new2).*255);
grid off
axis off
colorbar

%%
range =10;
new3 = NaN(size(new));
sz = range;
while sz <=max(new(:))+10
    
    new2 = new>=sz-range & new<=sz;
    bw = bwlabel(new2);
    b=1;
    while b<= max(bw(:))
    %     figure
    %     imshow(bw==b);
        if sum(sum(bw==b)) > 10
           [rows, cols] =find(bw==b);
           g=1;
           size(rows);
           while g<=max(size(rows))
           new3(rows(g),cols(g)) = 1;
           g=g+1;
           end
        end

        b=b+1;
    end
    sz = sz +range;
end

figure
imshow(new3);
% figure
% imshow(new>=0 & new <=20);
% colorbar
% 
% figure
% imshow(new>20 & new <=40);
% 
% figure
% imshow(new>=40 & new <=60);
% 
% figure
% imshow(new>=60 & new <=80);
% 
% figure
% imshow(new >80);

