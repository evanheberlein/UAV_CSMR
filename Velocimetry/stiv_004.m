%Basic m-file to display tiff images collected as individual tiff

clear
close all

s_path='C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\MicaSense\0004SET\';

Nconvert=4; %accounts for 14-bits images being stored in top 16 bits

iPairSep=1; %number of images across pair (2 = adjacent images)

%Loop over images and display them
therm_imgs = dir(fullfile(s_path, '**\*6.tif'));
n_img = str2num(sprintf('%03d', length(therm_imgs)));
start = 50;
stop = 500;
ratio_4 = [0.264292635658915,0.536917098445596];
x_size = size(imread(therm_imgs(1)),2)
col_L = 
col_df = zeros(col_L, stop-start);
col_x = round(x_size*ratio_4(1));
j = 1;
for i=start:stop %Must be 3 digits, beginning w/ 001 - originally 050:550
    fname=[therm_imgs(i).folder '\' therm_imgs(i).name];
        disp(fname)
        I1=double(imread(fname)/Nconvert); %loads image into a double precision real variable matrix
        I1_nm = I1 - mean(I1);
        col_df(:,j) = I1_nm(21:100, col_x);
        j = j+1;
end
figure(1)
imagesc(col_df)

col_df_nm = col_df - mean(col_df)

sq_start = 230;
sq_stop = sq_start + col_L;
Rxy = ifft2(abs(fft2(col_df(:,sq_start:sq_stop))).^2);
Rxy = Rxy - Rxy(40,40);

figure(2)
imagesc(real(log(Rxy)))

% Create random noise randn