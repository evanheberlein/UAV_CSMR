clear
close all
% AquaDopp HR data loading
% Load data
beam1d1 = importdata('0523_vel.v1')';
beam1d1 = beam1d1(3:end,:);
beam2d1 = importdata('0523_vel.v2')';
beam2d1 = beam2d1(3:end,:);
beam3d1 = importdata('0523_vel.v3')';
beam3d1 = beam3d1(3:end,:);
beam1d2 = importdata('0524_vel.v1')';
beam1d2 = beam1d2(3:end,:);
beam2d2 = importdata('0524_vel.v2')';
beam2d2 = beam2d2(3:end,:);
beam3d2 = importdata('0524_vel.v3')';
beam3d2 = beam3d2(3:end,:);
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
metdata = readtable('/Users/evanheberlein/Library/CloudStorage/Box-Box/Cornell/Fall_2021/Carp_UAV/carp-may-metdata.csv');
mettime = table2array(metdata(:,1));
metdepth = table2array(metdata(:,2));

thermdata = load('/Users/evanheberlein/Library/CloudStorage/Box-Box/Cornell/Fall_2021/Carp_UAV/Thermistor/Carp_P4.mat');

phi = 60; %Angle of oblique head
blnk_dist = 0.096; %m
beamdist = load('AD_beamdist.txt'); % beamdist col 1 = index, col 2 = beam (25deg off vert), col 3 = vertical
dbeam_vert = beamdist(2,3) - beamdist(1,3);
ylabs = [beamdist(1,3):dbeam_vert:beamdist(75,3)];
yticks = linspace(beamdist(1,3), beamdist(75,3), 75);

% Heatmap of bin correlations in time:
c1d1 = importdata('0523_vel.c1')'; % transpose ' to get vertical bins in y
c1d1 = c1d1(3:end,:);
c2d1 = importdata('0523_vel.c2')';
c2d1 = c2d1(3:end,:);
c3d1 = importdata('0523_vel.c3')';
c3d1 = c3d1(3:end,:);
c1d2 = importdata('0524_vel.c1')';
c1d2 = c1d2(3:end,:);
c2d2 = importdata('0524_vel.c2')';
c2d2 = c2d2(3:end,:);
c3d2 = importdata('0524_vel.c3')';
c3d2 = c3d2(3:end,:);

% Same w/ amplitudes:
a1d1 = importdata('0523_vel.a1')'; % transpose ' to get vertical bins in y
a1d1 = a1d1(3:end,:);
a2d1 = importdata('0523_vel.a2')';
a2d1 = a2d1(3:end,:);
a3d1 = importdata('0523_vel.a3')';
a3d1 = a3d1(3:end,:);
a1d2 = importdata('0524_vel.a1')';
a1d2 = a1d2(3:end,:);
a2d2 = importdata('0524_vel.a2')';
a2d2 = a2d2(3:end,:);
a3d2 = importdata('0524_vel.a3')';
a3d2 = a3d2(3:end,:);

% Get indices of flight times 
d1_flight_t = [find(d1time == datetime(2022,05,23,17,15,00))
    find(d1time == datetime(2022,05,23,17,25,00)) 
    find(d1time == datetime(2022,05,23,18,31,15))
    find(d1time == datetime(2022,05,23,18,41,15))];

d2_flight_t = [find(d2time == datetime(2022,05,24,08,00,15))
    find(d2time == datetime(2022,05,24,08,10,15))
    find(d2time == datetime(2022,05,24,08,50,00))
    find(d2time == datetime(2022,05,24,09,00,00))
    find(d2time == datetime(2022,05,24,10,06,45))
    find(d2time == datetime(2022,05,24,10,16,45))];

% Peak amplitude response bin every 1000 seconds:
% cdata_d1 = [2000 67; 3000 65; 4000 62; 5000 60; % Took off header rows, now off by 2?
%     6000 56.5; 7000 54; 8000 51; 9000 48; 10000 47;
%     11000 46; 12000 44; 13000 44; 14000 42; 15000 42]; % flipud bin #
% x_d1 = cdata_d1(:,1);
% y_d1 = cdata_d1(:,2);
% curve_d1 = fit(x_d1(:), y_d1, "sin2");
% curve_d1_t = string(d1time(x_d1));
% d1_t_axis = [1:1800:16201];
% d1_t_labs = ["14:30" "15:00" "15:30" "16:00" "16:30" "17:00" "17:30" "18:00" "18:30" "19:00"];
% Get indices of start and stop times:

% Surface elevation + curve fit + flight times - day 1
% figure()
% plot(x_d1, y_d1, 'o')
% hold on
% plot(curve_d1)
% hold off
% xticks(d1_t_axis)
% xticklabels(d1_t_labs)
% xline(d1_flight_t)
% legend('Location','southwest')

% Peak amplitude response bin every 1000 seconds:
% cdata_d2 = [2000 48; 3000 50; 4000 51; 5000 52; % Took off header rows, now off by 2? 
%     6000 54; 7000 55; 8000 56; 9000 58; 10000 59;
%     11000 60; 12000 61; 13000 62]; % flipud bin # - 1 = 77, 75 = 3 
% x_d2 = cdata_d2(:,1);
% y_d2 = cdata_d2(:,2);
% curve_d2 = fit(x_d2(:), y_d2, "sin2");
% curve_d2_t = string(d2time(x_d2));
% d2_t_axis = [1:1800:12601];
% d2_t_labs = ["7:00" "7:30" "8:00" "8:30" "9:00" "9:30" "10:00" "10:30"];
% Get indices of start and stop times:


% Surface elevation + curve fit + flight times - day 2
% figure()
% plot(x_d2, y_d2, 'o')
% hold on
% plot(curve_d2)
% hold off
% xticks(d2_t_axis)
% xticklabels(d2_t_labs)
% xline(d2_flight_t)
% legend('Location','southeast')
% dateaxis('x', 15, datetime(2022,05,24,14,30,00)) % type 15 = hh:mm

% take max amplitude bin in vertical, find corresponding oblique bin (-1,
% -2 etc), convert to velocity using trig formula

in_d1 = 1801;
out_d1 = 14401;
% d1_ticks = d1time([1801:1800:14401])
in_d2 = 1831;
out_d2 = 13517;

% figure()
% a = imagesc(flipud(a2d1(:,in_d1:out_d1)));
% colormap(spring(256))
% b = colorbar;
% xlabel('Time (May 23, 2022; hh:mm)')
% ylabel('Distance from AquaDopp (m)')
% ylabel(b, 'Amplitude reflected (#; max. 255)')
% % xticks(a2d1(:,d1_t_axis(2:end-1)))
% % xticklabels(d1_t_labs(2:end-1))
% % xticks = d1_ticks
% % xlabels = cellstr(datestr(d1time(xticks), datefmt))
% % xline(d1_flight_t-1801)
% set(gca, 'FontSize', 28, 'Xtick', [1:1800:out_d1-1800], 'XTickLabel', d1_t_labs(2:end-1),'YTick', [1:10:71], 'YTickLabel', flip(ylabs(5:10:75)))

% Visualizing amplitude responses from each beam/day
% figure()
% imagesc(flipud(a2d1(:,:)))
% colorbar
% 
% figure()
% imagesc(flipud(a2d2(:,:)))
% colorbar
% 


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

% Find max amplitude value and index (bin) in beam 2:
% remove flipud, +2 for top 2 rows not shown in fig. (bin # from imagesc)
[amax_d1 binmax_d1] = max(a2d1(7:37,:)); % surface b/w bins 7&37 - removed header +2
[amax_d2 binmax_d2] = max(a2d2(12:32,:)); % surface b/w bins 12&32 - removed header +2
binmax_d1 = binmax_d1 + (7); % don't account for +2 header rows here
binmax_d2 = binmax_d2 + (12);

% binmax_f0 +2 % need to account for top 2 header rows? 75 bins
% binz = zeros(length(beamdist(:,1),1));
beamdist(:,4) = ((beamdist(:,3)-blnk_dist)./cos(deg2rad(phi)))+blnk_dist; % beamdist col 3 = vertical
% beamdist col4 = oblique beam dist to h of vert beam (col. 3);
% Length of oblique beam for a given vertical beam h
% col 4 calculates diagonal distance to corresponding h of each vert beam
% Sanity: (beamdist{4,3}-blnk_dist)/(beamdist{4,4}-blnk_dist) = cosd(60)
for i = 1:length(beamdist(:,1)) % col5 = index of vert beam 3 bin corresponding to surface in oblique beams
    [mindiff beamdist(i,5)] = min(abs(beamdist(i,4)-beamdist(:,3))); %col3 original vertical bin
    x = beamdist(i,5);
    beamdist(i,6) = beamdist(x,3); % col6 = distance value of vert beam corresponding to index in col5
end
beamdist(:,7) = beamdist(:,6) - beamdist(:,4); %col7 = difference b/w vert&oblique bin centers (same to bin 36)
% beamdist(:,8) = beamdist(:,3).*cos(deg2rad(phi)); %col8 = h to each obl beam bin
% beamdist.Properties.VariableNames = ["Bin number", "30deg beam dist (m)", "Vert beam 2 h (m)",...
%     "Obl dist to vert h (m)","Vert-corresponding obl bin", "Col 6", "Dif", "Oblique bin h (m)"];

figure()
imagesc(a2d1(:,:))
colorbar
set(gca,'YDir','normal')
xline(d1_flight_t)
title('Day 1 beam 2 amplitude response')

figure()
imagesc(c2d1(:,:))
b = colorbar;
ylabel(b, 'Correlation %')
set(gca,'YDir','normal')
xline(d1_flight_t,'r','LineWidth',3)
hold on
plot(binmax_d1-2,'r') % subtract 2 for header rows
hold off
title('Day 1 beam 2 correlation')

figure()
imagesc(c1d1(:,:)+c3d1(:,:))
colorbar
set(gca,'YDir','normal')
xline(d1_flight_t,'r','LineWidth',3)
hold on
plot(beamdist(binmax_d1, 5),'r')
hold off
title('Day 1 beam 1+3 correlation')

figure()
imagesc(a2d2(:,:))
colorbar
set(gca,'YDir','normal')
xline(d2_flight_t)
title('Day 2 beam 2 amplitude response')

figure()
imagesc(c2d2(:,:))
b = colorbar;
ylabel(b, 'Day 2 correlation %')
set(gca,'YDir','normal')
xline(d2_flight_t,'r','LineWidth',3) 
hold on
plot(binmax_d2-2,'r')% subtract 2 for header rows
hold off
title('Day 2 beam 2 correlation')

figure()
imagesc(c1d2(:,:)+c3d2(:,:))
colorbar
set(gca,'YDir','normal')
xline(d2_flight_t,'r','LineWidth',3)
hold on
plot(beamdist(binmax_d2, 5),'r')
hold off
title('Day 2 beam 1+3 correlation')

% horizontal velocities - from Veliz Carillo 2021
u_d1 = (beam1d1 - beam3d1)./cos(phi); 
u_d2 = (beam1d2 - beam3d2)./cos(phi);
% vertical velocities
v_d1 = -(beam1d1 + beam3d1)./(2*sin(phi)); 
v_d2 = -(beam1d2 + beam3d2)./(2*sin(phi));

vbar_beam_d1 = mean(beam2d1(:,in_d1:out_d1),'All');
vbar_trig_d1 = mean(v_d1(:,in_d1:out_d1),'All');
vbar_beam_d2 = mean(beam2d2(:,in_d2:out_d2),'All');
vbar_trig_d2 = mean(v_d2(:,in_d2:out_d2),'All');
vdiff_d1 = abs(vbar_beam_d1 - vbar_trig_d1);
vdiff_d2 = abs(vbar_beam_d2 - vbar_trig_d2);

% Velocity filtering using correlation 
% Day 1 filtering
u_filt_d1 = u_d1;
% u_filt_d1(beamdist(binmax_d1(:), 5),:) = NaN;
v_filt_d1 = v_d1;
beam2_filt_d1 = beam2d1;
v_diff_d1 = NaN(1, length(v_d1));
c_thr = 50;
v_thr = 0.05;
for i = 2:length(u_d1)
    for j = 1:height(u_d1)
        if c1d1(j,i) > c_thr & c3d1(j,i) > c_thr; %& c2d1(j,i) > c_thr & abs(v(j,i)-beam2d1(j,i)) > v_thr
            u_filt_d1(j,i) = u_d1(j,i);
            if c2d1(j,i) > c_thr % want to compare v's when all 3 beams are good
                v_filt_d1(j,i) = v_d1(j,i);
                beam2_filt_d1(j,i) = beam2d1(j,i);
            else
                v_filt_d1(j,i) = NaN;
                beam2_filt_d1(j,i) = NaN;
            end
        else
            u_filt_d1(j,i) = NaN;%u_filt_d1(j,i-1);
        end
        if abs(v_d1(j,i)) > 0.4
            v_filt_d1(j,i) = NaN;
            u_filt_d1(j,i) = NaN;
        end
        if j > binmax_d1(i)
            beam2_filt_d1(j:end,i) = NaN;
        end
        if j > beamdist(binmax_d1(i), 5)
            u_filt_d1(j:end,i) = NaN;
            v_filt_d1(j:end,i) = NaN;
        end
    end
    v_diff_d1(i) = mean(abs(beam2_filt_d1(1:binmax_d1(i), i) - v_filt_d1(beamdist(binmax_d1(i), 5), i)));
end

avg_v_diff_d1 = mean(v_diff_d1, 'all', 'omitnan');

imAlpha_d1=ones(size(u_filt_d1));
imAlpha_d1(isnan(u_filt_d1))=0;

figure()
imagesc(u_filt_d1, 'AlphaData',imAlpha_d1)
colorbar
set(gca,'YDir','normal')
xline(d1_flight_t,'r','LineWidth',3)
% hold on
% plot(beamdist(binmax_d1, 5),'r')
% hold off
title('Day 1 filtered velocities')

figure()
histogram(v_filt_d1, EdgeColor = "none")
hold on
histogram(beam2_filt_d1, EdgeColor = "none")
hold off
legend("Calculated", "Measured")
title("Day 1 v-vel distribution")

% Day 2 filtering
u_filt_d2 = u_d2;
v_filt_d2 = v_d2;
beam2_filt_d2 = beam2d2;
for i = 2:length(u_d2)
    for j = 1:height(u_d2)
        if c1d2(j,i) > c_thr & c3d2(j,i) > c_thr; % See additional filter params d1
            u_filt_d2(j,i) = u_d2(j,i);
            if c2d2(j,i) > c_thr % want to compare v's when all 3 beams are good
                v_filt_d2(j,i) = v_d2(j,i);
                beam2_filt_d2(j,i) = beam2d2(j,i);
            else
                v_filt_d2(j,i) = NaN;
                beam2_filt_d2(j,i) = NaN;
            end
        else
            u_filt_d2(j,i) = NaN;%u_filt_d2(j,i-1);
        end
        if abs(v_d2(j,i)) > 0.4
            v_filt_d2(j,i) = NaN;
            u_filt_d2(j,i) = NaN;
        end
        % if abs(beam2d1(j,i)) > 0.09
        %     v_filt_d2(j,i) = NaN;
        %     u_filt_d2(j,i) = NaN;
        % end
        if j > binmax_d2(i)
            beam2_filt_d2(j,i) = NaN;
        end
        if j > beamdist(binmax_d2(i), 5)
            u_filt_d2(j,i) = NaN;
            v_filt_d2(j,i) = NaN;
        end
    end
end

avg_v_diff_d2 = mean(abs(v_filt_d2 - beam2_filt_d2), 'all', 'omitnan');

imAlpha_d2=ones(size(u_filt_d2));
imAlpha_d2(isnan(u_filt_d2))=0;

figure()
imagesc(u_filt_d2, 'AlphaData',imAlpha_d2)
colorbar
set(gca,'YDir','normal')
xline(d2_flight_t,'r','LineWidth',3)
% hold on
% plot(beamdist(binmax_d2, 5),'r')
% hold off
title('Day 2 filtered velocities')

figure()
histogram(v_filt_d2, EdgeColor = "none")
hold on
histogram(beam2_filt_d2, EdgeColor = "none")
hold off
legend("Calculated", "Measured")
title("Day 2 v-vel distribution")


% curve fitting vertical velocity profile
mean_surf0 = round(mean(binmax_d1(d1_flight_t(1):d1_flight_t(2)))); % Vertical to surface
vel_bin_f0 = NaN(beamdist(mean_surf0, 5), 1); % Create 1D vector of NaNs for each bin - removed header +2
for j = 1:beamdist(mean_surf0, 5) % From first bin to surface in oblique bin - removed header +2
    %sum(u_filt_d1(j,d1_flight_t(1):d1_flight_t(2)), 'omitnan')
    if sum(u_filt_d1(j,d1_flight_t(1):d1_flight_t(2)), 'omitnan') == 0 % All NaNs sums to 0
        vel_bin_f0(j) = NaN; % removed header -2
    else
        vel_bin_f0(j) = mean(u_filt_d1(j,d1_flight_t(1):d1_flight_t(2)), 'omitnan'); % Omit NaNs calculates mean of remaining numbers, removed header -2
    end
end

figure()
plot(vel_bin_f0)
xlabel('bin number')
ylabel('velocity (m/s)')

mean_surf1 = round(mean(binmax_d1(:,d1_flight_t(3):d1_flight_t(4))));

% For each flight, extract every 3 bins down from surface 
% Day 1 flights: flight 000
vel_f0 = zeros(6, length(d1_flight_t(1):d1_flight_t(2)));
for i = 1:length(vel_f0)
    j = i + d1_flight_t(1) - 1; % Create second index based on time offset b/w AD & flight
    obl_bin = beamdist(binmax_d1(j),5); % removed header +2
    vel_f0(1,i) = u_d1(obl_bin, j);
    vel_f0(2,i) = u_d1(obl_bin-3, j);
    vel_f0(3,i) = u_d1(obl_bin-6, j);
    vel_f0(4,i) = u_d1(obl_bin-9, j);
    vel_f0(5,i) = u_d1(obl_bin-12, j);
    vel_f0(6,i) = u_d1(obl_bin-15, j);
    % c_throld filter attempt:
end
for i = 2:length(vel_f0)
    if vel_f0(:,i) < 0
        vel_f0(:,i) = vel_f0(:,i-1);
    end
end

figure()
for i = 1:size(vel_f0,1)
    hold on
    plot(vel_f0(i,:));
    hold on
    plot(mean(u_filt_d1(:, d1_flight_t(1):d1_flight_t(2))));
end
hold off
legend()
title('Flight 000')

% flight 001
vel_f1 = zeros(6, length(d1_flight_t(3):d1_flight_t(4)));
for i = 1:length(vel_f1)
    j = i + d1_flight_t(3) - 1;
    obl_bin = beamdist(binmax_d1(j),5); % removed header +2
%     vel_f1(1,i) = u_d1(obl_bin, j);
%     vel_f1(2,i) = u_d1(obl_bin-1, j);
%     vel_f1(3,i) = u_d1(obl_bin-2, j);
%     vel_f1(4,i) = u_d1(obl_bin-10, j);
    vel_f1(1,i) = u_d1(obl_bin, j);
    vel_f1(2,i) = u_d1(obl_bin-3, j);
    vel_f1(3,i) = u_d1(obl_bin-6, j);
    vel_f1(4,i) = u_d1(obl_bin-9, j);
    vel_f1(5,i) = u_d1(obl_bin-12, j);
    vel_f1(6,i) = u_d1(obl_bin-15, j);
    % c_throld filtering: done globally
    if vel_f1(:,i) > -0.05
        vel_f1(:,i) = vel_f1(:,i-1);
    end
    if vel_f1(:,i) < -0.4
        vel_f1(:,i) = vel_f1(:,i-1);
    end
end
Lbin1 = beamdist(mean_surf1,5); %get avg obl_bin
Lbeam1 = beamdist(Lbin1, 3); % col 3 - vertical
Lsurf1 = 2*Lbeam1*sin(phi*pi/180); %m - surface length b/w velocity measurement points
u = mean(vel_f1(1,:));
u1 = mean(vel_f1(1,80:199));

figure()
for i = 1:size(vel_f1,1)
    hold on
    plot(vel_f1(i,:));
end
hold on
plot(mean(u_filt_d1(:, d1_flight_t(3):d1_flight_t(4))))
hold off
legend()
title('Flight 001')

% Running average velocity figure for AGU poster:
postervel = abs(movmean(vel_f1(1,80-60:535),[60 0]));
Pdatevec = d1time(d1_flight_t(3)+80:d1_flight_t(3)+535);
figure()
plot(Pdatevec, postervel(61:end), 'LineWidth', 4);
% plot(Pdatevec, vel_f1(1,80:535));
xlabel('Time (s)')
ylabel('ū-velocity (m/s)', 'interpreter', 'latex')
xlim([Pdatevec(1) Pdatevec(end)])
ylim([.1 .4])
set(gca, 'FontSize', 28)

% Day 2 flights: flight 004
vel_f4 = zeros(6, length(d2_flight_t(1):d2_flight_t(2)));
for k = 2:length(vel_f4)
    i = k-1;
    j = i + d2_flight_t(1) - 1;
    obl_bin = beamdist(binmax_d2(j),5); % removed header +2
    vel_f4(1,i) = u_d2(obl_bin, j);
    vel_f4(2,i) = u_d2(obl_bin-3, j);
    vel_f4(3,i) = u_d2(obl_bin-6, j);
    vel_f4(4,i) = u_d2(obl_bin-9, j);
    vel_f4(5,i) = u_d2(obl_bin-12, j);
    vel_f4(6,i) = u_d2(obl_bin-15, j);
end
% Attempt at acceleration filter:
for i = 2:length(vel_f4)
    if abs(vel_f4(1,i)-vel_f4(1,i-1)) > .05
        vel_f4(:,i) = vel_f4(:,i-1); %max([vel_f4(1,i-1) vel_f4(1,i)])
    end
end

figure()
for i = 1:size(vel_f4,1)
    hold on
    plot(vel_f4(i,:));
end
hold on
plot(mean(u_filt_d2(:, d2_flight_t(1):d2_flight_t(2)), ...
    'omitnan'), 'LineWidth', 3);
hold off
legend()
title('Flight 004')

% Flight 005
vel_f5 = zeros(6, length(d2_flight_t(3):d2_flight_t(4)));
for i = 1:length(vel_f5)
    j = i + d2_flight_t(3) - 1;
    obl_bin = beamdist(binmax_d2(j),5); % removed header +2
    vel_f5(1,i) = u_d2(obl_bin, j);
    vel_f5(2,i) = u_d2(obl_bin-3, j);
    vel_f5(3,i) = u_d2(obl_bin-6, j);
    vel_f5(4,i) = u_d2(obl_bin-9, j);
    vel_f5(5,i) = u_d2(obl_bin-12, j);
    vel_f5(6,i) = u_d2(obl_bin-15, j);
end
% Attempt at acceleration filter:
for i = 2:length(vel_f5)
    if abs(vel_f5(1,i)-vel_f5(1,i-1)) > .05
        vel_f5(:,i) = vel_f5(:,i-1); 
    end
end

figure()
for i = 1:size(vel_f5,1)
    hold on
    plot(vel_f5(i,:));
end
hold off
legend()
title('Flight 005')

% Flight 6
vel_f6 = zeros(4, length(d2_flight_t(5):d2_flight_t(6)));
for i = 1:length(vel_f6)
    j = i + d2_flight_t(5) - 1;
    obl_bin = beamdist(binmax_d2(j),5); % removed header +2
    vel_f6(1,i) = u_d2(obl_bin, j);
    vel_f6(2,i) = u_d2(obl_bin-3, j);
    vel_f6(3,i) = u_d2(obl_bin-6, j);
    vel_f6(4,i) = u_d2(obl_bin-9, j);
    vel_f6(5,i) = u_d2(obl_bin-12, j);
    vel_f6(6,i) = u_d2(obl_bin-15, j);
end
% Attempt at acceleration filter:
for i = 2:length(vel_f6)
    if abs(vel_f6(1,i)-vel_f6(1,i-1)) > .05
        vel_f6(:,i) = vel_f6(:,i-1); 
    end
end

figure()
for i = 1:size(vel_f6,1)
    hold on
    plot(vel_f6(i,:));
end
hold off
legend()
title('Flight 006')

% Create spectra & histogram for each flight
flights = {[vel_f0]; [vel_f1]; [vel_f4]; [vel_f5]; [vel_f6]};
names = ['0';'1';'4';'5';'6'];
for i = 1:length(flights)
    vel = cell2mat(flights(i));
    N = length(vel);
    Raa = 1./N.* ifft(fft(vel).*conj(fft(vel)));
    Saa6 = abs(fft(Raa));
    % Spectra
    figure()
    loglog(Saa6')
    title(sprintf('Flight %s spectra', names(i)))
    % Histogram
    figure()
    histogram(vel(1,:), 151)
    title(sprintf('Flight %s histogram - top bin', names(i)))
end

save('AD_process.mat')


% To do:
% -convert time columns in .sen file to time vector
% -cross-ref flight times & data quality w/ correlations
% -how to determine actual depth?
% Calculate v & look for non-zero mean
% Remove header rows from all data/operations