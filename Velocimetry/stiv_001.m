clear
close all

s_path='C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\MicaSense\0001SET\';
load('find_AD1.mat')
load('AD_process.mat')

Nconvert=4; %accounts for 14-bits images being stored in top 16 bits
iPairSep=1; %number of images across pair (2 = adjacent images)

%Loop over images and display them
therm_imgs = dir(fullfile(s_path, '**\*6.tif'));
n_img = str2num(sprintf('%03d', length(therm_imgs)));
start = 80; % first image of stationary hover
stop = 535; % last image of stationary hover 535
img1 = [therm_imgs(1).folder '\' therm_imgs(1).name];
% find column location over AquaDopp based on ratio size b/w IR and visible
IR_size = size(imread(img1)); % [y x]
ADloc_IR = [round(ADratio(1)*IR_size(2)) round(ADratio(2)*IR_size(1))];
ADsurfL = round(Lsurf1*pix_per_mv);

% Set up empty vectors for STIV
diam = min([ADsurfL IR_size(1)]);
col_df = zeros(diam, stop-start);
col_x = round(IR_size(2)*ADratio(1));
row_y = round(IR_size(1)*ADratio(2));

% Subtract minimum image of set
% rot_size = size(double(imrotate(imread(fname), ADang)))
% fname1=[therm_imgs(1).folder '\' therm_imgs(1).name];
% I1_min = zeros(rot_size(2), rot_size(1), stop-start); % 137x172 rotated dimension
% for i=start:stop %Must be 3 digits, beginning w/ 001 
%     fname=[therm_imgs(i).folder '\' therm_imgs(i).name];
%     j = i-start+1;
%         I1_minraw=double(imrotate(imread(fname), ADang));
%         I1_minraw(find(I1_minraw==0)) = NaN;
%         I1_min(:,:,j) = I1_minraw;
% end
% I1_min = min(I1_min, [], 3);

% Create STIV figure
j = 1;
for i=start:stop % Loop through imagery: must be 3 digits, beginning w/ 001 
    fname=[therm_imgs(i).folder '\' therm_imgs(i).name];
%         disp(fname)
        I1=double(imrotate(imread(fname), ADang)); %loads image into a double precision real variable matrix
        non0inds = find(I1(:,col_x)); % imrotate uses zero padding
        I1(find(I1==0)) = NaN; % convert zero-pad to nan so it doesn't affect mean
        I1 = I1 - mean(I1, 'all', 'omitnan'); % can -I1_min (min image) here
        col_df(:,j) = I1(non0inds, col_x);% can alternatively subtract column mean?
        j = j+1;
end
    
figure(1)
imagesc(col_df)
a=colorbar
xlabel('Time (May 23, 2022; hh:mm:ss)')
xticks([25:30:445])
xticklabels(datestr(d1time([d1_flight_t(3)+start+25:30:d1_flight_t(3)+start+30*14+25]),13))
yticklabels(flip([0:20:100]./pix_per_mir))
ylabel('y-distance (m)')
ylabel(a, 'Image-normalized pixel intensity (#)')
set(gca, 'FontSize', 28)

% Manual slope calculations
s11 = [215 120];
s12 = [275 75];
s1t = s12(1) - s11(1);
s1d = (s11(2) - s12(2))/pix_per_mir;
s1v = s1d/s1t

s21 = [155 67];
s22 = [215 17];
s2t = s22(1) - s21(1);
s2d = (s21(2) - s22(2))/pix_per_mir;
s2v = s2d/s2t

s31 = [15 111];
s32 = [76 49];
s3t = s32(1) - s31(1);
s3d = (s31(2) - s32(2))/pix_per_mir;
s3v = s3d/s3t

s41 = [312 83];
s42 = [372 24];
s4t = s42(1) - s41(1);
s4d = (s41(2) - s42(2))/pix_per_mir;
s4v = s4d/s4t

s51 = [205 120];
s52 = [265 71];
s5t = s52(1) - s51(1);
s5d = (s51(2) - s52(2))/pix_per_mir;
s5v = s5d/s5t

pause

% Automated STIV from Fujita et al. 2019 - work in progress
sq1 = 1;
Rxy1 = ifft2(abs(fft2(col_df(:,sq1:diam))).^2);
Rxy1 = Rxy1./Rxy1(diam/2, diam/2);
sq2= 177;
Rxy2 = ifft2(abs(fft2(col_df(:,sq2:sq2+diam-1))).^2);
Rxy2 = Rxy2./Rxy2(diam/2, diam/2);
sq3= 320;
Rxy3 = ifft2(abs(fft2(col_df(:,sq3:sq3+diam-1))).^2);
Rxy3 = Rxy3./Rxy3(diam/2, diam/2);

figure(2)
imagesc(real(log(Rxy1)))

Rxy_pol = zeros(diam);
xCenter = round(diam/2);
yCenter = round(diam/2);
radius = floor((diam-1)/2);

% circle tutorial https://matlab.fandom.com/wiki/FAQ#How_do_I_create_a_circle.3F
theta_c = linspace(0, 2*pi, round(4 * pi * radius));
xc = radius * cos(theta_c) + xCenter;
yc = radius * sin(theta_c) + yCenter;

circ_mask = poly2mask(xc, yc, diam, diam);
circ_nan = ones(length(circ_mask));
circ_nan(find(circ_mask==0)) = NaN;
theta = zeros(diam);
rho = zeros(diam);
Tx = zeros(diam);
Ty = zeros(diam);
M = 15;
% circ_mask(find(circ_mask==0)) = NaN;
for i = 1:diam
    for j = 1:diam
        if circ_mask(j,i) == 1
            Tx(j,i) = xCenter - i;
            Ty(j,i) = yCenter - j;
            theta(j,i) = atan(Ty(j,i)/Tx(j,i));
            rho(j,i) = M*log(sqrt(Tx(j,i)^2+Ty(j,i)^2));                
        end
    end
end     
Rxy1_pol = Rxy1.*circ_nan;
Rxy2_pol = Rxy2.*circ_nan;
Rxy3_pol = Rxy3.*circ_nan;

figure(3)
imagesc(Rxy1_pol)

rho = rho.*circ_nan;
theta = theta.*circ_nan;

figure(4)
imagesc(theta)

figure(5)
imagesc(rho)

% [Theta, Rho] = meshgrid
theta_deg = theta.*(180/pi);%+90; % +90 adds pi/2 - removed from eqn 11
pol_x = exp(rho).*cos(theta);
pol_y = exp(rho).*sin(theta);

figure(6)
pcolor(rho, theta_deg, Rxy1_pol)
colorbar
rho_round = round(rho, 1);

% in x/y grid, aggregate by polar coords
halfcirc = [-90:1:90];%[0:1:180];
theta1_int = zeros(length(halfcirc),1);
theta2_int = zeros(length(halfcirc),1);
theta3_int = zeros(length(halfcirc),1);
deg_start = halfcirc(1);
for m = 1:length(halfcirc)
    n = halfcirc(m);
    inds = find(round(theta_deg(:)) == n); % query cells from Rxy @ that index
    theta1_int(m) = sum(Rxy1_pol(inds)); % add up indexed cell values - numerical integral
    theta2_int(m) = sum(Rxy2_pol(inds));
    theta3_int(m) = sum(Rxy3_pol(inds));
end

max_rho = M*log(min(max(Tx(:)), max(Ty(:))));
mu1 = (1/max_rho).*theta1_int;
mu2 = (1/max_rho).*theta2_int;
mu3 = (1/max_rho).*theta3_int;
mu1_fit = fit([0:1:180]', mu1, "sin2");
mu2_fit = fit([0:1:180]', mu2, "sin2");
mu3_fit = fit([0:1:180]', mu3, "sin2");
fitx1 = mu1_fit([0:1:180]'); % to try sub-degree fit, change 1 to .1, use fitvec
fitx2 = mu2_fit([0:1:180]');
fitx3 = mu3_fit([0:1:180]');
[fitmin1, fitind1] = max(fitx1);
[fitmin2, fitind2] = max(fitx2);
[fitmin3, fitind3] = max(fitx3);
% fitvec = [0:.1:180]'
fitang1 = pi/2-(halfcirc(fitind1)*pi/180); % use fitvec here for sub-degree fit
fitang2 = pi/2-(halfcirc(fitind2)*pi/180);
fitang3 = pi/2-(halfcirc(fitind3)*pi/180);
Sx = 1/pix_per_mv;
St = 1; % Hz;
u1 = -Sx/St*tan(fitang1)
u2 = -Sx/St*tan(fitang2)
u3 = -Sx/St*tan(fitang3)
umean = mean([u1 u2 u3])

figure(7)
plot(mu1, 'o')
hold on
plot(mu1_fit)
hold off
xlabel('\theta (degrees)', 'FontSize', 24)
ylabel('\mu (-)', 'FontSize', 24)
xlim([0 180])
% title('STIV 1')

figure(8)
plot(mu2, 'o')
hold on
plot(mu2_fit)
hold off
xlabel('degrees')
ylabel('\mu')
title('STIV 2')

figure(9)
plot(mu3, 'o')
hold on
plot(mu3_fit)
hold off
xlabel('degrees')
ylabel('\mu')
title('STIV 3')