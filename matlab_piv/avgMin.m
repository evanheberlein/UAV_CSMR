function [Bkgrnd1,Bkgrnd2,iMean]=avgMin(iFileStart,iPairSep,iFileStop,sImagePath,sPrefix,sSuffix,Nconvert)
%[Bkgrnd1,Bkgrnd2,iMeans]=avgMin(iFileStart,iPairSep,iFileStop,sImagePath,sPrefix,sSuffix)
%calculate the global lowest value to occurr at each pixel location
%in each image pair and return this as imagesin Bkgrnd1 and Bkgrnd2
%The imgage means are returned in iMeans
%
%single tiff images

%Initialize the Bkgrnd files to a larger number than is possible to record
sFile1=num2str(iFileStart,'%04d'); %must change this if number of numeric 
%values in filename is differnt than 3
fname=[sImagePath sPrefix sFile1 sSuffix];
disp(fname)
I1=double(imread(fname)); 
Bkgrnd1=ones(size(I1))*2^16;
Bkgrnd2=Bkgrnd1;
i=0;
for iFile=iFileStart:iPairSep:iFileStop
   i=i+1;
% Read in Images into arrays I1 and I2
   sFile1=num2str(iFile,'%04d'); %must change this if number of numeric 
%values in filename is differnt than 3
   fname=[sImagePath sPrefix sFile1 sSuffix];
   I1=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix  
  
   sFile1=num2str(iFile+1,'%04d');  %must change this if number of numeric 
%values in filename is differnt than 3
   fname=[sImagePath sPrefix sFile1 sSuffix];
   I2=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix  
   sFrame=sprintf('%04d',iFile+iPairSep-1);
   disp(['Reading image ' fname])
 
   iMean(i,1)=median(mean(I1));
   iMean(i,2)=median(mean(I2));
   
   Bkgrnd1=min(Bkgrnd1,I1);
   Bkgrnd2=min(Bkgrnd2,I2);
end

