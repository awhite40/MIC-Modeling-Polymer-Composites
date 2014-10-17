matrix=load('C:/Users/Lenovo/Documents/Gatech/OneDrive for Business/Gatech 4th Sem/Coursework/MI/Project/Dicom images/Scropped.mat');
cropped=matrix.cropped;
cropped1=double(cropped);
cropped2=cropped1(:);
[yy,xx] = hist(cropped2,151);

% peakfit(signal,center,window,NumPeaks,peakshape); 
% "peakshape" specifies the peak shape of the model: (1=Gaussian
% (default), 2=Lorentzian, 3=logistic distribution, 4=Pearson,
% 5=exponentionally broadened Gaussian; 6=equal-width Gaussians;
% 7=Equal-width Lorentzians; 8=exponentionally broadened equal-width
% Gaussian, 9=exponential pulse, 10=up-sigmoid (logistic function),
% 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian; 13=Gaussian/
% Lorentzian blend; 14=Bifurcated Gaussian, 15=Breit-Wigner-Fano,
% 16=Fixed-position Gaussians; 17=Fixed-position Lorentzians;
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile; 21=triangular; 22=multiple shapes; 23=down-sigmoid;
% 25=lognormal; 26=sine wave; 27=Gaussian first derivative.

% [FitResults,FitError]=peakfit(signal,center,window,NumPeaks,peakshape); 
% Returns the FitResults vector in the order 
% peak number, peak position, peak height,peak width, and peak area), 
% and the FitError (the percent RMS
% difference between the data and the model in the selected segment of that
% data) of the best fit.

figure
[FitResults,FitError]= peakfit([xx;yy],0,0,2,1);

[r,c]=size(FitResults);

Y=zeros(size(cropped1,1),size(cropped1,2),size(cropped1,3),r);

for i=1:r
    Y(1:size(cropped1,1),1:size(cropped1,2),1:size(cropped1,3),i) = FitResults(i,3)*exp(-((cropped1-FitResults(i,2)).^2/(2*(FitResults(i,4)).^2)));
end

[M,I]=max(Y,[],4);

figure
imshow(I(:,:,1),[1 2])

save 'SegSample1P2.mat' I

% figure
% s  = regionprops(I(:,:,1), 'centroid');
% centroids = cat(1, s.Centroid);
% imshow(I(:,:,1),[1 2])
% hold on
% plot(centroids(:,1), centroids(:,2), 'b*')
% hold off
