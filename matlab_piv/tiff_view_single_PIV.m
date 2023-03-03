%Basic m-file to display tiff images collected as individual tiff

% Last updated 10/15/21 by EAC

clear all
close all

s_path='/Users/eac/classes/CEE6370/labs/PIVLab2021/Wed06_2021/Images/';
prefix='wed06img_';
suffix='.tif';
iStart=1; %First image
iFiles=99; %number of files to loop over
Cmax=200; %maximum threshold value in pseudocolored image
Cmin=50;  %minimum threshold value in pseudocolored image
Nconvert=4; %accounts for 14-bits images being stored in top 16 bits

iPairSep=3; %number of images across pair (2 = adjacent images)

%Loop over images and display them
for i=iStart:iPairSep:iStart+iFiles-1
    sFile1=num2str(i,'%04d');
    fname=[s_path prefix sFile1 suffix];
    disp(fname)
    I1=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix
    sFile2=num2str(i+1,'%04d');
    fname=[s_path prefix sFile2 suffix];
    I2=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix
    min1=min(I1(:));
    max1=max(I1(:));
    min2=min(I2(:));
    max2=max(I2(:));
    mean1=mean(I1(:));
    mean2=mean(I2(:));
    median1=median(I1(:));
    median2=median(I2(:));
    disp(sprintf('Min1=%3d, Min2=%3d, Max1=%4d, Max2=%4d',min1,min2,max1,max2))
    disp(sprintf('Mean1=%6.2f, Mean2=%6.2f, Median1=%6.2f, Median2=%6.2f',mean1,mean2,median1,median2))
    
    figure(1)
    set(gcf,'position',[120 560 560 420])
    imagesc(I1,[median1+Cmin median1+Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title (num2str(i))
    
    figure(2)
    set(gcf,'position',[690 560 560 420])
    imagesc(I2,[median2+Cmin median2+Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title (num2str(i+1))
     
    figure(3)
    set(gcf,'position',[400 60 560 420])
    imagesc(I2-I1,[-Cmax Cmax])
    axis image
    colormap(gray(256))
    colorbar
    title ([sFile2 ' - ' sFile1])
    pause
end

