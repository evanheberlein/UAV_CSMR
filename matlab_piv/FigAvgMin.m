function FigAvgMin(Bckgrnd1, Bckgrnd2, iMean)
%FigAvgMin(Bckgrnd1, Bckgrnd2, iMean) - function to plot Bkgrnd images
%and the mean value of each image

figure(1)
set(gcf,'position',[20 400 1000 250])
subplot(1,3,1)
imagesc(Bckgrnd1)
axis image
title('Bckgrnd 1 (minimum image)')
colorbar
   
subplot(1,3,2)
imagesc(Bckgrnd2)
axis image
title('Bckgrnd 2 (minimum image)')
colorbar

subplot(1,3,3)
plot(iMean)
title('Global median values')
xlabel('Image pair')
ylabel('Global Median (counts)')
legend('odd files','even files')


  