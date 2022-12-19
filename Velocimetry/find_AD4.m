
% Find Aquadopp in visible light images
s_path='C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\MicaSense\0004SET\';

% blue (1) green (2) & red edge (5) are closest to thermal sensor: https://support.micasense.com/hc/en-us/articles/360010025413-Altum-Integration-Guide#h.vtwsbws4yz1x

% load images
file = "002\IMG_0400_"

blue = strcat(s_path, file, '1', '.tif')
green = strcat(s_path, file, '2', '.tif')
rededge = strcat(s_path, file, '5', '.tif')

figure()
imshow(blue)
title("blue")
ADbl = [539 816]
bl_bot = [565 882]
bl_top = [592 761]
bl_hyp = sqrt((bl_bot(1)-bl_top(1))^2 + (bl_bot(2)-bl_top(2))^2)
bl_ang = acos(((bl_bot(2)-bl_top(2))/bl_hyp)*180/pi)
bl_gcp4 = [1167 973]
bl_gcp5 = [1697 453]
bl_gcp_hyp = sqrt((bl_gcp5(1)-bl_gcp4(1))^2 + (bl_gcp5(2)-bl_gcp4(2))^2)
bl_gcp_rat = bl_gcp_hyp / (472/39.37) % 472" c to c, 39.37 in/m

figure()
imshow(green)
% title("green")
ADgr = [552 842]
gr_bot = [578 904]
gr_top = [602 795]
gr_hyp = sqrt((gr_bot(1)-gr_top(1))^2  + (gr_bot(2)-gr_top(2))^2)
gr_ang = acos((gr_bot(2)-gr_top(2))/gr_hyp)*180/(2*pi)
gr_gcp4 = [1179 1000]
gr_gcp5 = [1708 481]
gr_gcp_hyp = sqrt((gr_gcp5(1)-gr_gcp4(1))^2 + (gr_gcp5(2)-gr_gcp4(2))^2)
gr_gcp_rat = gr_gcp_hyp / (472/39.37)
pix_per_m = (gr_gcp_rat + bl_gcp_rat)/2

ADloc = [(539+552)/2 (816+842)/2]
dims = size(imread(blue)) % [y x]
ADratio = [ADloc(1)/dims(2) ADloc(2)/dims(1)]
ADang = (bl_ang + gr_ang)/2

figure()
imshow(rededge)
title("red edge")
% Difficult to see

save('find_AD4.mat')