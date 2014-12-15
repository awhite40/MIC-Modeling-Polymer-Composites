function [value] =  Thresh_3D( Matrix )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a = size(Matrix)
X = Matrix
%Thresholding in X-Y direction (X = rows Y = Columns)
for z = 50:2:60
    c = X(:,:,z);
    minc = min(min(c));
    maxc = max(max(c));
    %Preallocating the size of the vectors
    b=1;
    inte = (zeros([(maxc-minc)/1000, 1]));%edit denominator for finer thresholding
    values = int16(zeros([(maxc-minc)/1000,1]));%edit denominator for finer thresholding
    d = int16(zeros([nRows, nColumns]));
    value = minc;
    %Scanning through threshold values
    while value < maxc
        i = 1;
        while i <=nRows
            j=1;
            while j <= nColumns
                if c(i,j,1) < value % Values less than the threshold are set to zero
                    d(i,j)= 0;
                else
                    d(i,j)= 1; %Everything else is set to 1
                end
                j=j+1;
            end
            i=i+1;
        end
        %Adding values into vectors
        inte(b) = sum(sum(d));
        values(b)=value;
        b = b+1;
        %Save new array to workspace
        assignin('base','slice_edit',d);
    value = value + 1000;%Edit this for finer thresholding    
    end
%     assignin('base','intensity',inte); 
%     assignin('base','threshold_values',values);

%Display plot of intensity vs. threshold value
    assignin('base','rows',nRows*nColumns);
    assignin('base','inte',inte);
    inte = inte/(double(nRows)*double(nColumns));
    plot(values,inte); hold on;
    xlabel('Threshold values')
    ylabel('Percentage of pixels')
end
for x = 50:2:60
    c = X(x,:,:);
    %Preallocating the size of the vectors
    b=1;
    intex = zeros([(maxc-minc)/1000, 1]);%edit denominator for finer thresholding
    values = int16(zeros([(maxc-minc)/1000,1]));%edit denominator for finer thresholding
    dx = zeros([nColumns,nFrames]);
    value = minc;
    %Scanning through threshold values
    while value < maxc
        i = 1;
        while i <=nColumns
            j=1;
            while j <= nFrames
                if c(1,i,j) < value % Values less than the threshold are set to zero
                    dx(i,j)= 0;
                else
                    dx(i,j)= 1; %Everything else is set to 1
                end
                j=j+1;
            end
            i=i+1;
        end
        %Adding values into vectors
        intex(b) = sum(sum(dx));
        values(b)=value;
        b = b+1;
        %Save new array to workspace
        assignin('base','slice_editx',dx);
        
        
    value = value + 1000;%Edit this for finer thresholding    
    end
%     assignin('base','intensity',inte); 
%     assignin('base','threshold_values',values);

%Display plot of intensity vs. threshold value
    intex = intex/(double(nFrames)*double(nColumns));
    plot(values,intex,'red'); hold on;
end
for y = 50:2:60
    c = X(:,y,:);
    Preallocating the size of the vectors
    b=1;
    intey = zeros([(maxc-minc)/1000, 1]);%edit denominator for finer thresholding
    values = int16(zeros([(maxc-minc)/1000,1]));%edit denominator for finer thresholding
    dy = zeros([nRows,nFrames]);
    value = 250;
    Scanning through threshold values
    while value < 300
        i = 1;
        while i <=nRows
            j=1;
            while j <= nFrames
                if c(i,1,j) < value % Values less than the threshold are set to zero
                    dy(i,j)= 0;
                else
                    dy(i,j)= 1; %Everything else is set to 1
                end
                j=j+1;
            end
            i=i+1;
        end
        Adding values into vectors
        intey(b) = sum(sum(dy));
        values(b)=value;
        b = b+1;
        Save new array to workspace
        assignin('base','slice_editx',dx);
        
        
    value = value + 1000;%Edit this for finer thresholding    
    end
    assignin('base','intensity',inte); 
    assignin('base','threshold_values',values);

%Display plot of intensity vs. threshold value
    
    intey = intey/(double(nFrames)*double(nColumns));
    plot(values,intey,'green');
    legend('X-Y Plane','X-Z Plane','Y-Z Plane')
end

figure
imshow(dy)
end


