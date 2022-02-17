function FigImagesAndDifferences(I1,I2,Cmin, Cmax)
%FigImagesAndDifferences(I1,I2) - function to plot images (I1 and I2)
%in an image pair and their differences (I2 - I1)

figure(3)
set(gcf,'position',[20 200 1000 250])
subplot(1,3,1)
imagesc(I1,[Cmin Cmax])
axis image
title('Image 1')
colorbar
   
subplot(1,3,2)
imagesc(I2,[Cmin Cmax])
axis image
title('Image 2')
colorbar

subplot(1,3,3)
imagesc(I2-I1,[-Cmax Cmax])
axis image
title('Image 2 - Image 1')
colorbar



  