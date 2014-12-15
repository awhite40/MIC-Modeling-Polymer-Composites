cropped = double(Ncropped_300);
[yy,xx] = hist( cropped(:),151);
%plot(xx,yy);
[ax, h1,h2] = plotyy( xx, yy,... Plot distribution
    xx, [gradient(yy); gradient(gradient(yy))] ... Plot gradients of distribution
    )
figure(gcf);
legend([h1;h2], 'Raw','First Derivative','Second Derivative')
grid( ax(1), 'on')
set( [h1;h2], 'LineWidth',3)
set( ax,'Fontsize',16)
ylabel(ax(1), 'Raw Data')
ylabel(ax(2), 'Gradients')
figure
[ p,e ]= peakfit( [xx;yy], 0,0,3,1 );
