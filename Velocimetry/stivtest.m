clear
close all
x = (0:1:1000);
y = zeros(length(x));
y1 = zeros(length(x));
y2 = zeros(length(x));
for i = 1:length(x)
    y1(i,:) = sin(.15*(x-i/2));
    y2(i,:) = 0.9*sin(.05*(x-i));
    y(i,:) = y1(i,:)+y2(i,:);
end

% y = repmat(y1+y2, length(y1), 1);
figure(1)
imagesc(y)

Rxy = ifft2(abs(fft2(y(:,:))).^2);
figure(2)
imagesc(real(log(Rxy)))

Rxy_pol = zeros(length(x));
% 0c_grid = (-500:1:500);

xCenter = 501;
yCenter = 501;
radius = 500;
% Circumference for a circle of radius 350 should be 2*pi*r = 2199 pixels.
% To ensure that we have no gaps in the circle 
% we need to make sure we have at least as many coordinates in vectors x and y 
% as there are around the circumference of the circle.
% Make it double that just to make extra sure there are no gaps in the circle
% by going all 360 degrees (0 - 2*pi) with 4398 points.
theta_c = linspace(0, 2*pi, round(4 * pi * radius)); % Define angles
% Get x and y vectors for each point along the circumference.
xc = radius * cos(theta_c) + xCenter;
yc = radius * sin(theta_c) + yCenter;

% 
% % Write those (x,y) into the image with gray level 255.
% for k = 1 : length(x)
%     row = round(yc(k));
%     col = round(xc(k));
%     Rxy_pol(row, col) = Rxy(row, col);
% end
circ_mask = poly2mask(xc, yc, 1001, 1001);
circ_nan = ones(length(circ_mask));
circ_nan(find(circ_mask==0)) = NaN;
theta = zeros(length(x));
rho = zeros(length(x));
Tx = zeros(length(x));
Ty = zeros(length(x));
M = 15; % "coefficient of intensification"
% circ_mask(find(circ_mask==0)) = NaN;
for i = 1:length(x)
    for j = 1:length(y)
        if circ_mask(j,i) == 1
            Tx(j,i) = xCenter - i;
            Ty(j,i) = yCenter - j;
            theta(j,i) = atan(Ty(j,i)/Tx(j,i));
            rho(j,i) = M*log(sqrt(Tx(j,i)^2+Ty(j,i)^2));                
        end
    end
end     
Rxy_pol = Rxy.*circ_nan;
figure(3)
% plot(xc, yc);
% axis square;
% grid on;
% hold on
imagesc(Rxy_pol)

rho = rho.*circ_nan;
theta = theta.*circ_nan;

figure(4)
imagesc(theta)

figure(5)
imagesc(rho)

% [Theta, Rho] = meshgrid
theta_deg = theta.*(180/pi) + 90; % make all values positive
pol_x = exp(rho).*cos(theta);
pol_y = exp(rho).*sin(theta);

figure(6)
pcolor(rho, theta_deg, Rxy_pol)
colorbar
rho_round = round(rho, 1);

% in x/y grid, aggregate by polar coords
theta_int = zeros(181, 1);
deg_start = 1%-90;
deg_end = deg_start + 180
% for m = (deg_start:1:90)
%     n = m-deg_start+1
% %     p = sum(round(theta_deg(:)) == m) % number of cells matching theta
% %     theta_i = zeros(p);
% %     for q = 1:p; % find index of each cell @ that theta
%     inds = find(round(theta_deg(:)) == m); % query cells from Rxy @ that index
%     theta_int(m+91) = sum(Rxy_pol(inds));
% end

inds_vec = zeros(deg_end);
inds_counter = 0
for m = (deg_start:1:deg_end)
%     p = sum(round(theta_deg(:)) == m) % number of cells matching theta
%     theta_i = zeros(p);
%     for q = 1:p; % find index of each cell @ that theta
    inds = find(round(theta_deg(:)) == m); % query cells from Rxy @ that index
    inds_vec(m) = length(inds);
    theta_int(m) = sum(Rxy_pol(inds));
end

% max_rho = M*log(min(max(Tx), max(Ty)));
% mu_theta = 1/max_rho.*theta_int;

figure()
plot(inds_vec)

figure()
plot(theta_int)

Rxy_norm = Rxy(501,501);
theta_int_norm = theta_int./(Rxy_norm*mean(inds_vec));

figure()
plot(theta_int_norm)

figure()
imagesc(y1)

figure()
imagesc(y2)

% rho_vec = rho(501,1:1001);
% d = [0:1:1000].*(1001/360);