%Basic m-file to display tiff images collected as individual tiff

clear
close all

s_path='C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\MicaSense\0000SET\';

Nconvert=4; %accounts for 14-bits images being stored in top 16 bits

iPairSep=1; %number of images across pair (2 = adjacent images)

%Loop over images and display them
therm_imgs = dir(fullfile(s_path, '**\*6.tif'));
n_img = str2num(sprintf('%03d', length(therm_imgs)));
start = 141;
stop = 536;
col_df = zeros(80, stop-start);
col = 60;
j = 1;
for i=start:stop %Must be 3 digits, beginning w/ 001
    fname=[therm_imgs(i).folder '\' therm_imgs(i).name];
        disp(fname)
        I1=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix
        I1_nm = I1 - mean(I1);
        col_df(:,j) = I1_nm(21:100, col);
        j = j+1;
end
imagesc(col_df)