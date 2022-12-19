clear
close all
% AquaDopp HR data loading
% Load data
beam1d1 = importdata('0523_vel.v1');
beam2d1 = importdata('0523_vel.v2');
beam3d1 = importdata('0523_vel.v3');
beam1d2 = importdata('0524_vel.v1');
beam2d2 = importdata('0524_vel.v2');
beam3d2 = importdata('0524_vel.v3');
sen_d1 = importdata('0523_vel.sen');
sen_d2 = importdata('0524_vel.sen');
hr2_d1 = importdata('0523_vel.hdr');
hr2_d2 = importdata('0524_vel.hdr');

% Recreate time vectors for data streams
d1time = (datetime('2022-05-23 14:30:00') : seconds(1) : datetime('2022-05-23 19:00:41')).';
d1time.Format = 'yyyy-MM-dd HH:mm:ss';
d2time = (datetime('2022-05-24 07:00:00') : seconds(1) : datetime('2022-05-24 10:54:57')).';
d2time.Format = 'yyyy-MM-dd HH:mm:ss';
thermtime = (datetime('2022-05-23 07:00:00') : seconds(20) : datetime('2022-05-24 09:11:00')).';
thermtime.Format = 'yyyy-MM-dd HH:mm:ss';

% Import met data and create vectors
metdata = readtable('C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\Carp-May-metdata.csv');
mettime = table2array(metdata(:,1));
metdepth = table2array(metdata(:,2));

thermdata = load('C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\Thermistor\Carp_P4.mat');

phi = 60; %Angle of oblique head
blnk_dist = 0.096; %m
beamdist = load('AD_beamdist.txt'); % beamdist col 1 = index, col 2 = beam, col 3 = vertical
dbeam_vert = beamdist(2,3) - beamdist(1,3);
ylabs = [beamdist(1,3):dbeam_vert:beamdist(75,3)];
yticks = linspace(beamdist(1,3), beamdist(75,3), 75);

% Heatmap of bin correlations in time:
% c1d1 = importdata('0523_vel.c1')'; % transpose ' to get vertical bins in y
% c2d1 = importdata('0523_vel.c2')';
% c3d1 = importdata('0523_vel.c3')';
% c1d2 = importdata('0524_vel.c1')';
% c2d2 = importdata('0524_vel.c2')';
% c3d2 = importdata('0524_vel.c3')';

a1d1 = importdata('0523_vel.a1')'; % transpose ' to get vertical bins in y
a2d1 = importdata('0523_vel.a2')';
a3d1 = importdata('0523_vel.a3')';
a1d2 = importdata('0524_vel.a1')';
a2d2 = importdata('0524_vel.a2')';
a3d2 = importdata('0524_vel.a3')';

cdata_d1 = [2000 67; 3000 65; 4000 62; 5000 60;
    6000 56.5; 7000 54; 8000 51; 9000 48; 10000 47;
    11000 46; 12000 44; 13000 44; 14000 42; 15000 42]; % flipud bin #
x_d1 = cdata_d1(:,1);
y_d1 = cdata_d1(:,2);
curve_d1 = fit(x_d1(:), y_d1, "sin2");
curve_d1_t = string(d1time(x_d1));
d1_t_axis = [1:1800:16201];
d1_t_labs = ["14:30" "15:00" "15:30" "16:00" "16:30" "17:00" "17:30" "18:00" "18:30" "19:00"];
d1_flight_t = [find(d1time == datetime(2022,05,23,17,15,00))
    find(d1time == datetime(2022,05,23,17,25,00)) 
    find(d1time == datetime(2022,05,23,18,31,15))
    find(d1time == datetime(2022,05,23,18,41,15))];
% fo = fitoptions('Method','NonlinearLeastSquares');
% ft = fittype('sin2','options',fo);

figure()
plot(x_d1, y_d1, 'o')
hold on
plot(curve_d1)
hold off
xticks(d1_t_axis)
xticklabels(d1_t_labs)
xline(d1_flight_t)
legend('Location','southwest')

cdata_d2 = [2000 48; 3000 50; 4000 51; 5000 52;
    6000 54; 7000 55; 8000 56; 9000 58; 10000 59;
    11000 60; 12000 61; 13000 62]; % flipud bin # - 1 = 77, 75 = 3 
x_d2 = cdata_d2(:,1);
y_d2 = cdata_d2(:,2);
curve_d2 = fit(x_d2(:), y_d2, "sin2");
curve_d2_t = string(d2time(x_d2));
d2_t_axis = [1:1800:12601];
% t1_d2 = datetime(2022,05,24,07,00,00)
% t2_d2 = datetime(2022,05,24,10,30,00)
d2_t_labs = ["7:00" "7:30" "8:00" "8:30" "9:00" "9:30" "10:00" "10:30"];
d2_flight_t = [find(d2time == datetime(2022,05,24,08,00,15))
    find(d2time == datetime(2022,05,24,08,10,15))
    find(d2time == datetime(2022,05,24,08,50,00))
    find(d2time == datetime(2022,05,24,09,00,00))
    find(d2time == datetime(2022,05,24,10,06,45))
    find(d2time == datetime(2022,05,24,10,16,45))];

figure()
plot(x_d2, y_d2, 'o')
hold on
plot(curve_d2)
hold off
xticks(d2_t_axis)
xticklabels(d2_t_labs)
xline(d2_flight_t)
legend('Location','southeast')
% dateaxis('x', 15, datetime(2022,05,24,14,30,00)) % type 15 = hh:mm

% take max amplitude bin in vertical, find corresponding oblique bin (-1,
% -2 etc), convert to velocity using trig formula

% figure()
% imagesc(flipud(c1d1(3:77,:)))
% colorbar
% % axes normal/xy help axes
% set(gca, 'Ytick', yticks, 'YTickLabel', ylabs);

% figure()
% imagesc(flipud(c2d1(3:77,:)))
% colorbar

in_d1 = 1801
out_d1 = 14401
% d1_ticks = d1time([1801:1800:14401])
figure()
a = imagesc(flipud(a2d1(3:77,in_d1:out_d1)))
colormap(spring(256))
b = colorbar
xlabel('Time (May 23, 2022; hh:mm)')
ylabel('Distance from AquaDopp (m)')
ylabel(b, 'Amplitude reflected (#; max. 255)')
% xticks(a2d1(:,d1_t_axis(2:end-1)))
% xticklabels(d1_t_labs(2:end-1))
% xticks = d1_ticks
% xlabels = cellstr(datestr(d1time(xticks), datefmt))
% xline(d1_flight_t-1801)
set(gca, 'FontSize', 28, 'Xtick', [1:1800:out_d1-1800], 'XTickLabel', d1_t_labs(2:end-1),'YTick', [1:10:71], 'YTickLabel', flip(ylabs(5:10:75)))


figure()
imagesc(flipud(a3d1(3:77,:)))
colorbar

figure()
imagesc(flipud(a1d2(3:77,:)))
colorbar

figure()
imagesc(flipud(a2d2(3:77,:)))
colorbar
xline(d2_flight_t)

figure()
imagesc(flipud(a3d2(3:77,:)))
colorbar

figure()
imagesc(flipud(a1d2(3:77,:)))
colorbar

% figure()
% imagesc(beam2d1(:,3:77)')
% colorbar
% 
% figure()
% imagesc(flipud(beam1d1(:,3:77)'))
% colorbar

% figure()
% plot(median(beam1d1(:,65:75),2))
% 
% % Comparison of depth from P thermistor vs. met station (depth converted
% % but not calibrated?)
% figure()
% plot(mettime, metdepth)
% hold on
% plot(thermtime, thermdata.pressure)
% hold off

[amax_f0 binmax_f0] = max(a2d1(77-70:77-40,d1_flight_t(1):d1_flight_t(2))); % subtract flipud bin # from 77 to get actual bin w/ 1 @ bottom
[amax_f1 binmax_f1] = max(a2d1(77-70:77-40,d1_flight_t(3):d1_flight_t(4)));
[amax_f4 binmax_f4] = max(a2d2(77-65:77-45,d2_flight_t(1):d2_flight_t(2)));
[amax_f5 binmax_f5] = max(a2d2(77-65:77-45,d2_flight_t(3):d2_flight_t(4)));
[amax_f6 binmax_f6] = max(a2d2(77-65:77-45,d2_flight_t(5):d2_flight_t(6)));
binmax_f0 = binmax_f0 + (77-70);
binmax_f1 = binmax_f1 + (77-70);
binmax_f4 = binmax_f4 + (77-65);
binmax_f5 = binmax_f5 + (77-65);
binmax_f6 = binmax_f6 + (77-65);

% binmax_f0 + 2 % need to account for top 2 header rows? 75 bins
% binz = zeros(length(beamdist(:,1),1));
beamdist(:,4) = beamdist(:,3)./cos(phi*pi/180); % col4 = oblique beam dist to surface (based on vert beam 3)
for i = 1:length(beamdist(:,1)) % col5 = index of vert beam 3 bin corresponding to surface in oblique beams
    [mindiff beamdist(i,5)] = min(abs(beamdist(i,4)-beamdist(:,3))); %col3 original vertical bin
    x = beamdist(i,5);
    beamdist(i,6) = beamdist(x,3); % col6 = distance value of vert beam corresponding to index in col5
end
beamdist(:,7) = beamdist(:,6) - beamdist(:,4); %col7 = difference b/w vert&oblique bin centers (same to bin 36)

vel_d1 = (beam1d1' - beam3d1')./cos(phi);
vel_d2 = (beam1d2' - beam3d2')./cos(phi);

vel_f0 = zeros(4, length(binmax_f0));
for i = 1:length(vel_f0)
    j = i + d1_flight_t(1);
    obl_bin = beamdist(binmax_f0(i),5)+2; % +2 accounts for header rows in vel tables
    vel_f0(1,i) = vel_d1(obl_bin, j);
    vel_f0(2,i) = vel_d1(obl_bin-1, j);
    vel_f0(3,i) = vel_d1(obl_bin-2, j);
    vel_f0(4,i) = vel_d1(obl_bin-10, j);
end

figure()
plot(vel_f0(1,:));
hold on
plot(vel_f0(2,:));
hold on
plot(vel_f0(3,:));
hold on
plot(vel_f0(4,:));
hold off
legend()

figure()
histogram(vel_f0(1,:), 151)

% flight 001
vel_f1 = zeros(4, length(binmax_f1));
for i = 1:length(vel_f1)
    j = i + d1_flight_t(3);
    obl_bin = beamdist(binmax_f1(i),5)+2; % +2 accounts for header rows in vel tables
    vel_f1(1,i) = vel_d1(obl_bin, j);
    vel_f1(2,i) = vel_d1(obl_bin-1, j);
    vel_f1(3,i) = vel_d1(obl_bin-2, j);
    vel_f1(4,i) = vel_d1(obl_bin-10, j);
    if vel_f1(:,i) > -0.05
        vel_f1(:,i) = vel_f1(:,i-1);
    end
    if vel_f1(:,i) < -0.4
        vel_f1(:,i) = vel_f1(:,i-1);
    end
end
Lbin1 = beamdist(round(mean(binmax_f1)),5); %get avg obl_bin
Lbeam1 = beamdist(Lbin1, 3); % col 3 - vertical
Lsurf1 = 2*Lbeam1*sin(phi*pi/180); %m
u = mean(vel_f1(1,:))
u1 = mean(vel_f1(1,80:199))

figure()
plot(vel_f1(1,:));
hold on
plot(vel_f1(2,:));
hold on
plot(vel_f1(3,:));
hold on
plot(vel_f1(4,:));
hold off
legend()
title('Flight 001')

postervel = abs(movmean(vel_f1(1,80-60:535),[60 0]))
Pdatevec = d1time(d1_flight_t(3)+80:d1_flight_t(3)+535);
figure()
plot(Pdatevec, postervel(61:end), 'LineWidth', 4);
% plot(Pdatevec, vel_f1(1,80:535));
xlabel('Time (s)')
ylabel('ū-velocity (m/s)', 'interpreter', 'latex')
xlim([Pdatevec(1) Pdatevec(end)])
ylim([.1 .4])
set(gca, 'FontSize', 28)

pause

vel_f4 = zeros(4, length(binmax_f4));
obl_bin = zeros(1, length(binmax_f4));
for i = 1:length(vel_f4)
    j = i + d2_flight_t(1);
    obl_bin = beamdist(binmax_f4(i),5)+2; % +2 accounts for header rows in vel tables
    vel_f4(1,i) = vel_d2(obl_bin, j);
    vel_f4(2,i) = vel_d2(obl_bin-1, j);
    vel_f4(3,i) = vel_d2(obl_bin-2, j);
    vel_f4(4,i) = vel_d2(obl_bin-10, j);
end

figure()
plot(vel_f4(1,:));
hold on
plot(vel_f4(2,:));
hold on
plot(vel_f4(3,:));
hold on
plot(vel_f4(4,:));
hold off
legend()

vel_f5 = zeros(4, length(binmax_f5));
for i = 1:length(vel_f5)
    j = i + d2_flight_t(3);
    obl_bin = beamdist(binmax_f5(i),5)+2; % +2 accounts for header rows in vel tables
    vel_f5(1,i) = vel_d2(obl_bin, j);
    vel_f5(2,i) = vel_d2(obl_bin-1, j);
    vel_f5(3,i) = vel_d2(obl_bin-2, j);
    vel_f5(4,i) = vel_d2(obl_bin-10, j);
end

figure()
plot(vel_f5(1,:));
hold on
plot(vel_f5(2,:));
hold on
plot(vel_f5(3,:));
hold on
plot(vel_f5(4,:));
hold off
legend()

vel_f6 = zeros(4, length(binmax_f6));
for i = 1:length(vel_f6)
    j = i + d2_flight_t(5);
    obl_bin = beamdist(binmax_f6(i),5)+2; % +2 accounts for header rows in vel tables
    vel_f6(1,i) = vel_d2(obl_bin, j);
    vel_f6(2,i) = vel_d2(obl_bin-1, j);
    vel_f6(3,i) = vel_d2(obl_bin-2, j);
    vel_f6(4,i) = vel_d2(obl_bin-10, j);
end

figure()
plot(vel_f6(1,:));
hold on
plot(vel_f6(2,:));
hold on
plot(vel_f6(3,:));
hold on
plot(vel_f6(4,:));
hold off
legend()

save('AD_process.mat')

% To do:
% -convert time columns in .sen file to time vector
% -cross-ref flight times & data quality w/ correlations
% -how to determine actual depth?
% 60 degrees!!!
% Calculate v & look for non-zero mean