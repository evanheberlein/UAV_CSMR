clear
close all

% Add AGW filter to path
addpath '/Users/evanheberlein/Library/CloudStorage/Box-Box/Cornell/Fall_2021/Exp_methods/Lab1/AGW-master'
addpath '/Users/evanheberlein/Library/CloudStorage/Box-Box/Cornell/Fall_2021/Exp_methods/Lab1/despiking_toolbox'
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

in_d1 = 1801;
out_d1 = 14401;
% d1_ticks = d1time([1801:1800:14401])
in_d2 = 1831;
out_d2 = 13517;

% Find max amplitude response bin in beam 2
[amax_d1 binmax_d1] = max(a2d1(7:37,:)); % surface b/w bins 7&37 - removed header +2
[amax_d2 binmax_d2] = max(a2d2(12:32,:)); % surface b/w bins 12&32 - removed header +2
binmax_d1 = binmax_d1 + (6); % add back starting index, don't account for +2 header rows here
binmax_d2 = binmax_d2 + (11);

% Beam geometry
phi = 60; %Angle of oblique head
blnk_dist = 0.096; %m
beamdist = load('AD_beamdist.txt'); % beamdist col 1 = index, col 2 = beam (25deg off vert), col 3 = vertical
beamdist = array2table(beamdist);
beamdist.beamdist2=[]; % Remove 30 degree beam distances - NA to this head

% beamdist(:,x) = ((beamdist(:,2)-blnk_dist)./cos(deg2rad(phi)))+blnk_dist;
% col. x gives distance along oblique to corresponding vertical bin
% Edited definition of col. x to account for same blnk_dist in obl.
% Bin sizes & blanking dist. assumed to be same 
% Sanity: (beamdist{4,2}-blnk_dist)/(beamdist{4,3}-blnk_dist) = cosd(60)
beamdist(:,3) = beamdist(:,2).*cos(deg2rad(phi)); % h of obl bins
for i = 1:length(beamdist{:,1}) % col5 = index of vert beam bin corresponding to oblique beam of same h
    [mindiff beamdist(i,4)] = min(abs(beamdist{i,2}-beamdist(:,3))); %col2 original vertical bin
    % To compare b/w vert/oblique of same h, vert 10 = oblique 23
    % Col 4 finds index of OBLIQUE bin corresponding to given VERTICAL bin
end

beamdist.Properties.VariableNames = ["Bin number","Vert beam 2 h (m)",...
    "Oblique beams 1&3 h (m)","Obl bin # for vert h"];

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

u_wat_d1 = u_d1;
v_wat_calc_d1 = v_d1;
v_wat_meas_d1 = beam2d1;
c_obl_d1 = (c1d1 + c3d1)./2;
a_obl_d1 = (a1d1 + a3d1)./(2*255).*100; % Normalize to maximum response

for i = 2:length(u_d1)
    obl_binmax = beamdist{binmax_d1(i), 4};
    for j = 1:height(u_d1)
        if j > obl_binmax % If vert. bin > surf, 
            u_wat_d1(j,i) = NaN;
            v_wat_calc_d1(j,i) = NaN;
            c_obl_d1(j,i) = NaN;
            a_obl_d1(j,i) = NaN;
        end
        if j > binmax_d1(i); % ID'd & applied to vertical bin - same
            v_wat_meas_d1(j,i) = NaN;
        end
    end
end

imAlpha_d1=ones(size(u_wat_d1));
imAlpha_d1(isnan(u_wat_d1))=0;

figure()
imagesc(u_wat_d1, 'AlphaData',imAlpha_d1)
colorbar
set(gca,'YDir','normal')
xline(d1_flight_t,'r','LineWidth',3)
hold on
plot(beamdist{binmax_d1, 4},'r')
hold off
title('Day 1 raw velocities')
yticks(5:5:height(u_wat_d1))
yticklabels(beamdist{:,3}(5:5:end))
xticks(2000:2000:16000)
xticklabels(datestr(d1time(2000:2000:16000),'HH:MM'))
ylabel("Depth (m)")
fontsize(16,'points')

u_wat_d2 = u_d2;
v_wat_calc_d2 = v_d2;
v_wat_meas_d2 = beam2d2;
c_obl_d2 = (c1d2 + c3d2)./2;
a_obl_d2 = (a1d2 + a3d2)./(2*255).*100; % Normalize to maximum response
for i = 2:length(u_d2)
    obl_binmax = beamdist{binmax_d2(i), 4};
    for j = 1:height(u_d2)
        if j > obl_binmax
            u_wat_d2(j,i) = NaN;
            v_wat_calc_d2(j,i) = NaN;
            c_obl_d2(j,i) = NaN;
            a_obl_d2(j,i) = NaN;
        end
        if j > binmax_d2(i);
            v_wat_meas_d2(j,i) = NaN;
        end
    end
end

imAlpha_d2=ones(size(u_wat_d2));
imAlpha_d2(isnan(u_wat_d2))=0;

figure()
imagesc(u_wat_d2, 'AlphaData',imAlpha_d2)
colorbar
set(gca,'YDir','normal')
xline(d2_flight_t,'r','LineWidth',3)
hold on
plot(beamdist{binmax_d2, 4},'r')
hold off
title('Day 2 raw velocities')
yticks(5:5:height(u_wat_d1))
yticklabels(beamdist{:,3}(5:5:end))
xticks(2000:2000:14000)
xticklabels(datestr(d2time(2000:2000:14000),'HH:MM'))
ylabel("Depth (m)")
fontsize(16,'points')

% Divide into flight times:
t_buffer = 60; % 1 min on either side
flight_us = NaN(height(beam1d1), ...
    d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1, 5);
flight_vs_calc = NaN(height(beam1d1), ...
    d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1, 5);
flight_vs_meas = NaN(height(beam1d1), ...
    d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1, 5);
flight_us_band = NaN(height(beam1d1)/15, ...
    d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1, 5);
flight_us_band_agw = NaN(height(beam1d1)/15, ...
    d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1, 5);
time_vecs = NaT(5, d1_flight_t(2)-d1_flight_t(1)+2*t_buffer+1);

flight_us(:,:,1) = u_wat_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
flight_us(:,:,2) = u_wat_d1(:,d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
flight_us(:,:,3) = u_wat_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
flight_us(:,:,4) = u_wat_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
flight_us(:,:,5) = u_wat_d2(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

flight_vs_calc(:,:,1) = v_wat_calc_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
flight_vs_calc(:,:,2) = v_wat_calc_d1(:,d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
flight_vs_calc(:,:,3) = v_wat_calc_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
flight_vs_calc(:,:,4) = v_wat_calc_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
flight_vs_calc(:,:,5) = v_wat_calc_d2(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

flight_vs_meas(:,:,1) = v_wat_meas_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
flight_vs_meas(:,:,2) = v_wat_meas_d1(:,d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
flight_vs_meas(:,:,3) = v_wat_meas_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
flight_vs_meas(:,:,4) = v_wat_meas_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
flight_vs_meas(:,:,5) = v_wat_meas_d2(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

flight_corrs(:,:,1) = c_obl_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
flight_corrs(:,:,2) = c_obl_d1(:,d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
flight_corrs(:,:,3) = c_obl_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
flight_corrs(:,:,4) = c_obl_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
flight_corrs(:,:,5) = c_obl_d2(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

flight_amps(:,:,1) = a_obl_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
flight_amps(:,:,2) = a_obl_d1(:,d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
flight_amps(:,:,3) = a_obl_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
flight_amps(:,:,4) = a_obl_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
flight_amps(:,:,5) = a_obl_d2(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

time_vecs(1,:) = d1time(d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer); %000
time_vecs(2,:) = d1time(d1_flight_t(3)-t_buffer:d1_flight_t(4)+t_buffer); %001
time_vecs(3,:) = d2time(d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer); %004
time_vecs(4,:) = d2time(d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer); %005
time_vecs(5,:) = d2time(d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer); %006

flight_names = ["Flight 0 (d1)", "Flight 1 (d1)", "Flight 4 (d2)", "Flight 5 (d2)", "Flight 6 (d2)"];
flight_us_agw = NaN(75,721,5);

for flight_ind = 1:5
    % Create histogram
    figure()
    u = histogram(flight_us(:,:,flight_ind));
    hold on
    vm = histogram(flight_vs_meas);
    hold on
    vc = histogram(flight_vs_calc(:,:,flight_ind));
    hold off
    u.EdgeColor = 'none';
    vc.EdgeColor = 'none';
    vm.EdgeColor = 'none';
    vc.FaceColor = 'green';
    legend("U","V (measured - beam 2)","V (calculated - beams 1 & 3)")
    title(append(flight_names(flight_ind), " Raw velocity distribution"))
    xlabel("Velocity (m/s)")
    ylabel("Count")
    fontsize(16, "points")
    grid on
    set(gcf,'position',[10,10,900,500])

     % AGW filtering
    for row = 1:75
        [row_filt agw_inds] = ...
            agw_filter(flight_us(row,:,flight_ind), (1:length(flight_us)), -5, 5);
        flight_us_agw(row, agw_inds, flight_ind) = row_filt;
    end

    % % 
    band_start = 1;
    leg_ind = 1;

    % Starting figure comparing raw & AGW-filtered data
    figure()
    for band_ind = 1:5
        % Plot raw data TS bands
        band_end = band_ind*15;
        flight_us_band(band_ind,:,flight_ind) = mean(...
            flight_us(band_start:band_end,:,flight_ind), 1, 'omitnan');
        scatter(time_vecs(flight_ind,:),...
            flight_us_band(band_ind,:,flight_ind), 200,'+')
        legendtxt{leg_ind} = ...
            [append(num2str(beamdist{band_start,3}), "-", ...
            num2str(beamdist{band_end,3}))];
        hold on

        leg_ind = leg_ind+1;

        flight_us_band_agw(band_ind,:,flight_ind) = mean(...
            flight_us_agw(band_start:band_end,:,flight_ind), 1, 'omitnan');
        plot(time_vecs(flight_ind,:),...
            flight_us_band_agw(band_ind,:,flight_ind),'LineWidth',3)
        legendtxt{leg_ind} = ...
            [append(num2str(beamdist{band_start,3}), "-", ...
            num2str(beamdist{band_end,3}), " AGW-filt")];
        hold on

        leg_ind = leg_ind+1;
        band_start = band_end+1;
    end
    hold off
    grid on
    xlabel("Time (HH:mm)")
    ylabel("Velocity (m/s)")
    leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
    title(leg, "Vertical profile band h from bottom (m)");
    set(gcf,'position',[10,10,1200,700])
    title(append(flight_names(flight_ind), " u-velocity time series"))
    fontsize(16, "points")

    % Create correlation histograms
    figure()
    histogram(flight_corrs(:,:,flight_ind));
    title(append(flight_names(flight_ind), " Raw oblique beam correlation distribution"))
    xlabel("Correlation % (/100)")
    ylabel("Count")
    fontsize(16, "points")
    grid on
    set(gcf,'position',[10,10,900,500])    
end


% Create vertical profiles of AGW-filtered data
figure()
for flight_ind = 1:5
    vert_prof = mean(flight_us_agw(:,:,flight_ind),2, 'omitnan');
    plot(vert_prof, beamdist{:,3}, 'LineWidth', 3);
    hold on
end
hold off
legend(flight_names, 'Location','southwest')
xlabel("u-velocity (m/s)")
ylabel("Vertical profile (m) [n = 75]")
title("Temporal-average vertical profiles for each flight (AGW-filtered only)")
fontsize(16, "points")
grid on
set(gcf,'position',[10,10,900,500])   

% Find average correlation %  for each flight
flight_corrs_mean = mean(flight_corrs,[1,2],'omitnan') ;
flight_amps_mean = mean(flight_amps,[1,2],'omitnan');

% Individual flights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight 0

% Try corr. filter first
flight_ind = 1;
flight0_us_filt = flight_us_agw(:,:,flight_ind);
corr_flight0 = c_obl_d1(:,d1_flight_t(1)-t_buffer:d1_flight_t(2)+t_buffer);
c_thresh = 35;

wrap_factor = 20.4;% [19.5:.1:20.5];

for row = 1:75
    wrapped_row = [flight0_us_filt(row,:)].*wrap_factor;
    unwrapped_row = unwrap(wrapped_row);
    flight0_us_filt(row,:) = unwrapped_row./wrap_factor;
end
 
for row = 1:height(flight0_us_filt)
    for i = 1:length(flight0_us_filt)
        if corr_flight0(row,i) < c_thresh;
            flight0_us_filt(row,i) = NaN;
        end
    end
end

band_start = 1;
leg_ind = 1;
figure()
for band_ind = 1:5
    band_end = band_ind*15;
    scatter(time_vecs(flight_ind,:),...
        flight_us_band_agw(band_ind,:,flight_ind), 200,'+')
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " AGW-filt")];
    hold on

    leg_ind = leg_ind+1;

    flight0_us_band(band_ind,:) = mean(...
        flight0_us_filt(band_start:band_end,:), 1, 'omitnan');

    % Manually correcting phase shift
    if band_ind == 4
        b4_phase_ind_1 = find(time_vecs(1,:) ==...
            "23-May-2022 17:19:47");
        b4_phase_ind_2 = find(time_vecs(1,:) ==...
            "23-May-2022 17:21:38");
        flight0_us_band(band_ind,b4_phase_ind_1:b4_phase_ind_2) =...
            flight0_us_band(4,b4_phase_ind_1:b4_phase_ind_2) +...
            (pi/wrap_factor); % Ok to shift by just pi instead of 2pi?
    end
     
    if band_ind == 5
        b5_phase_ind_1 = find(time_vecs(1,:) ==...
            "23-May-2022 17:21:03");
        b5_phase_ind_2 = find(time_vecs(1,:) ==...
            "23-May-2022 17:25:07");
        flight0_us_band(band_ind,b5_phase_ind_1:b5_phase_ind_2) =...
            flight0_us_band(5,b5_phase_ind_1:b5_phase_ind_2) +...
            (2*pi/wrap_factor);
        flight0_us_band(band_ind,b5_phase_ind_2:end) =...
            flight0_us_band(5,b5_phase_ind_2:end) +...
            (3*pi/wrap_factor); % Again, not a multiple of 2pi
    end

    % Acceleration filter of 0.07 (same as Flight 1)
    for i = 2:length(flight0_us_filt)
        if abs(flight0_us_band(band_ind, i) -...
                flight0_us_band(band_ind, i-1)) >= 0.07
            flight0_us_band(band_ind, i)  = NaN;
        end
    end

    plot(time_vecs(flight_ind,:),...
        flight0_us_band(band_ind,:),'LineWidth',3)
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " unwrap/acc")];
    hold on

    leg_ind = leg_ind+1;
    band_start = band_end+1;
end
hold off
grid on
xlabel("Time (HH:mm)")
ylabel("Velocity (m/s)")
leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
title(leg, "Vertical profile band h from bottom (m)");
set(gcf,'position',[10,10,1200,700])
title(append(flight_names(flight_ind), " u-velocity time series "))
fontsize(16, "points")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight 1
flight_ind = 2;

flight1_us_filt = NaN(size(flight_us_agw(:,:,flight_ind)));
% flight1_acc_tab = ones(size(flight_us_agw(:,:,2)));
wrap_factor = 20; % wrap only works in increments of pi:

% counter = 1;
for row = 1:75
    wrapped_row = [flight_us_agw(row,:,flight_ind)].*wrap_factor;
    unwrapped_row = unwrap(wrapped_row);
    flight1_us_filt(row,:) = unwrapped_row./wrap_factor;

    % Acceleration filter on individual rows: not so effective
    % for i = 2:length(flight1_us_filt)
    %     if abs(flight1_us_filt(row, i) - flight1_us_filt(row, i-1)) >= 0.05
    %         flight1_us_filt(row, i)  = NaN;
    %         % counter = counter+1
    %     end
    % end
end

% flight1_us_filt = flight1_us_filt.*flight1_acc_tab;

band_start = 1;
leg_ind = 1;
figure()
for band_ind = 1:5
    band_end = band_ind*15;
    scatter(time_vecs(flight_ind,:),...
        flight_us_band_agw(band_ind,:,flight_ind), 200,'+')
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " AGW-filt")];
    hold on

    leg_ind = leg_ind+1;

    flight1_us_band(band_ind,:) = mean(...
        flight1_us_filt(band_start:band_end,:), 1, 'omitnan');

    % Acceleration filter of 0.07m/s
    for i = 2:length(flight1_us_filt)
        if abs(flight1_us_band(band_ind, i) - flight1_us_band(band_ind, i-1)) >= 0.07
            flight1_us_band(band_ind, i)  = NaN;
            % counter = counter+1
        end
    end

    plot(time_vecs(flight_ind,:),...
        flight1_us_band(band_ind,:),'LineWidth',3)
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " unwrap/acc")];
    hold on

    leg_ind = leg_ind+1;
    band_start = band_end+1;
end
hold off
grid on
xlabel("Time (HH:mm)")
ylabel("Velocity (m/s)")
leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
title(leg, "Vertical profile band h from bottom (m)");
set(gcf,'position',[10,10,1200,700])
title(append(flight_names(flight_ind), " u-velocity time series"))
fontsize(16, "points")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight 4

flight_ind = 3;

flight4_us_filt = NaN(size(flight_us_agw(:,:,flight_ind)));

wrap_factor = 18.6; 
wrap_chunks = [1,180, 181,360, 361,540, 541,721, 0];
for row = 1:75
    wrapped_row = [flight_us_agw(row,:,flight_ind)].*wrap_factor;
    unwrapped_row = unwrap(wrapped_row);
    flight4_us_filt(row,:) = unwrapped_row./wrap_factor;

    %%% Does unwrapping in smaller chunks improve performance? No!
    % wrap_start = wrap_chunks(1);
    % for wrap_ts = 1:4
    %     wrap_end = wrap_chunks(wrap_ts *2);
    % 
    %     wrapped_row = [flight_us_agw(row,wrap_start:wrap_end,flight_ind)].*wrap_factor;
    %     unwrapped_row = unwrap(wrapped_row);
    %     flight4_us_filt(row,wrap_start:wrap_end) = unwrapped_row./wrap_factor;
    % 
    %     wrap_start = wrap_chunks(wrap_ts*2+1);
    % end
end

corr_flight4 = c_obl_d2(:,d2_flight_t(1)-t_buffer:d2_flight_t(2)+t_buffer);
c_thresh_1 = 78;
c_thresh_2 = 40;

for row = 1:30
    for i = 1:length(flight4_us_filt)
        if corr_flight4(row,i) < c_thresh_1;
            flight4_us_filt(row,i) = NaN;
        end
    end
end

% for row = 31:75
%     for i = 1:length(flight4_us_filt)
%         if corr_flight4(row,i) < c_thresh_2;
%             flight4_us_filt(row,i) = NaN;
%         end
%     end
% end


band_start = 1;
leg_ind = 1;
figure()
for band_ind = 1:5
    band_end = band_ind*15;
    scatter(time_vecs(flight_ind,:),...
        flight_us_band_agw(band_ind,:,flight_ind), 200,'+')
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " AGW-filt")];
    hold on

    leg_ind = leg_ind+1;

    flight4_us_band(band_ind,:) = mean(...
        flight4_us_filt(band_start:band_end,:), 1, 'omitnan');

    % Acceleration filter
    for i = 2:length(flight4_us_filt)
        if abs(flight4_us_band(band_ind, i) - flight4_us_band(band_ind, i-1)) >= 0.10
            flight4_us_band(band_ind, i)  = NaN;
            flight4_us_band(band_ind, i-1)  = NaN;
        end
    end

    % if band_ind == 2
    %     b2_phase_ind_1 = find(time_vecs(3,:) ==...
    %         "24-May-2022 08:07:59");
    %     b2_phase_ind_2 = find(time_vecs(3,:) ==...
    %         "24-May-2022 08:08:45");
    %     b2_phase_ind_3 = find(time_vecs(3,:) ==...
    %         "24-May-2022 08:09:52");
    %     flight4_us_band(band_ind,b2_phase_ind_1:b2_phase_ind_2) =...
    %         flight4_us_band(2,b2_phase_ind_1:b2_phase_ind_2) -...
    %         (pi/wrap_factor); 
    %     flight4_us_band(band_ind,b2_phase_ind_3:end) =...
    %         flight4_us_band(2,b2_phase_ind_3:end) +...
    %         (pi/wrap_factor); 
    % end

    if band_ind == 3
        b3_phase_ind_1 = find(time_vecs(3,:) ==...
            "24-May-2022 08:01:55");
        b3_phase_ind_2 = find(time_vecs(3,:) ==...
            "24-May-2022 08:04:42");
        b3_phase_ind_3 = find(time_vecs(3,:) ==...
            "24-May-2022 08:08:11");
        flight4_us_band(band_ind,b3_phase_ind_1:b3_phase_ind_2) =...
            flight4_us_band(band_ind,b3_phase_ind_1:b3_phase_ind_2) -...
            (pi/wrap_factor); 
        flight4_us_band(band_ind,b3_phase_ind_2+1:b3_phase_ind_3-1) =...
            flight4_us_band(band_ind,b3_phase_ind_2+1:b3_phase_ind_3-1) -...
            (3*pi/wrap_factor); 
        flight4_us_band(band_ind,b3_phase_ind_3:end) =...
            flight4_us_band(band_ind,b3_phase_ind_3:end) -...
            (5*pi/wrap_factor); 
    end

    if band_ind == 4
        b4_phase_ind_1 = find(time_vecs(3,:) ==...
            "24-May-2022 08:01:55");
        b4_phase_ind_2 = find(time_vecs(3,:) ==...
            "24-May-2022 08:04:42");
        b4_phase_ind_3 = find(time_vecs(3,:) ==...
            "24-May-2022 08:04:47");
        b4_phase_ind_4 = find(time_vecs(3,:) ==...
            "24-May-2022 08:08:01");
        b4_phase_ind_5 = find(time_vecs(3,:) ==...
            "24-May-2022 08:08:14");
        b4_phase_ind_6 = find(time_vecs(3,:) ==...
            "24-May-2022 08:08:58");
        % b4_phase_ind_7 = find(time_vecs(3,:) ==...
        %     "24-May-2022 08:09:03");
        flight4_us_band(band_ind,b4_phase_ind_1:b4_phase_ind_2) =...
            flight4_us_band(band_ind,b4_phase_ind_1:b4_phase_ind_2) -...
            (2*pi/wrap_factor); 
        flight4_us_band(band_ind,b4_phase_ind_2+1:b4_phase_ind_3) =...
            flight4_us_band(band_ind,b3_phase_ind_2+1:b4_phase_ind_3) -...
            (pi/wrap_factor); 
        flight4_us_band(band_ind,b4_phase_ind_3+1:b4_phase_ind_4) =...
            flight4_us_band(band_ind,b4_phase_ind_3+1:b4_phase_ind_4) -...
            (4*pi/wrap_factor); 
        flight4_us_band(band_ind,b4_phase_ind_4+1:b4_phase_ind_5) =...
            flight4_us_band(band_ind,b4_phase_ind_4+1:b4_phase_ind_5) -...
            (6*pi/wrap_factor); 
        flight4_us_band(band_ind,b4_phase_ind_5+1:b4_phase_ind_6) =...
            flight4_us_band(band_ind,b4_phase_ind_5+1:b4_phase_ind_6) -...
            (9*pi/wrap_factor); 
        flight4_us_band(band_ind,b4_phase_ind_6+1:end) =...
            flight4_us_band(band_ind,b4_phase_ind_6+1:end) -...
            (14*pi/wrap_factor);  
    end

    plot(time_vecs(flight_ind,:),...
        flight4_us_band(band_ind,:),'LineWidth',3)
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " unwrap/acc")];
    hold on

    leg_ind = leg_ind+1;
    band_start = band_end+1;
end
hold off
grid on
xlabel("Time (HH:mm)")
ylabel("Velocity (m/s)")
leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
title(leg, "Vertical profile band h from bottom (m)");
set(gcf,'position',[10,10,1200,700])
title(append(flight_names(flight_ind), " u-velocity time series "))
fontsize(16, "points")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight 5

flight_ind = 4;

flight5_us_filt = NaN(size(flight_us_agw(:,:,flight_ind)));

wrap_factor = 1; 
for row = 1:75
    wrapped_row = [flight_us_agw(row,:,flight_ind)].*wrap_factor;
    unwrapped_row = unwrap(wrapped_row);
    flight5_us_filt(row,:) = unwrapped_row./wrap_factor;
end

% Correlation filter:
corr_flight5 = c_obl_d2(:,d2_flight_t(3)-t_buffer:d2_flight_t(4)+t_buffer);
c_thresh_1 = 78;
c_thresh_2 = 60;

for row = 1:15
    for i = 1:length(flight5_us_filt)
        if corr_flight5(row,i) < c_thresh_1;
            flight5_us_filt(row,i) = NaN;
        end
    end
end

for row = 31:75
    for i = 1:length(flight4_us_filt)
        if corr_flight5(row,i) < c_thresh_2;
            flight5_us_filt(row,i) = NaN;
        end
    end
end


band_start = 1;
leg_ind = 1;
figure()
for band_ind = 1:5
    band_end = band_ind*15;
    scatter(time_vecs(flight_ind,:),...
        flight_us_band_agw(band_ind,:,flight_ind), 200,'+')
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " AGW-filt")];
    hold on

    leg_ind = leg_ind+1;

    flight5_us_band(band_ind,:) = mean(...
        flight5_us_filt(band_start:band_end,:), 1, 'omitnan');

    % Acceleration filter
    for i = 2:length(flight5_us_filt)
        if abs(flight5_us_band(band_ind, i) - flight5_us_band(band_ind, i-1)) >= 0.10
            flight5_us_band(band_ind, i)  = NaN;
            flight5_us_band(band_ind, i-1)  = NaN;
        end
    end

 
    plot(time_vecs(flight_ind,:),...
        flight5_us_band(band_ind,:),'LineWidth',3)
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " unwrap/acc")];
    hold on

    leg_ind = leg_ind+1;
    band_start = band_end+1;
end
hold off
grid on
xlabel("Time (HH:mm)")
ylabel("Velocity (m/s)")
leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
title(leg, "Vertical profile band h from bottom (m)");
set(gcf,'position',[10,10,1200,700])
title(append(flight_names(flight_ind), " u-velocity time series ", num2str(wrap_factor)))
fontsize(16, "points")

figure()
for i = 31:45
    plot(time_vecs(4,:),flight5_us_filt(i,:))
    hold on
end
hold off
title("Flight 5 bins 31:45 u-velocities")
xlabel("HH:MM")
ylabel("u (m/s)")
fontsize(16, "points")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight 6

flight_ind = 5;
flight6_us_filt = flight_us_agw(:,:,flight_ind);
corr_flight6 = c_obl_d1(:,d2_flight_t(5)-t_buffer:d2_flight_t(6)+t_buffer);
c_thresh = 60;

for row = 1:height(flight6_us_filt)
    for i = 1:length(flight6_us_filt)
        if corr_flight6(row,i) < c_thresh;
            flight6_us_filt(row,i) = NaN;
        end
    end
end

wrap_factor = 1;
for row = 1:75
    % wrapped_row = [flight_us_agw(row,:,flight_ind)].*wrap_factor;
    % unwrapped_row = unwrap(wrapped_row);
    % flight1_us_filt(row,:) = unwrapped_row./wrap_factor;

    % % Acceleration filter on individual rows: not so effective
    for i = 2:length(flight6_us_filt)
        if abs(flight6_us_filt(row, i) - flight6_us_filt(row, i-1)) >= 0.07
            flight6_us_filt(row, i)  = NaN;
            flight6_us_filt(row, i-1) = NaN;
            % counter = counter+1
        end
    end
end

band_start = 1;
leg_ind = 1;
figure()
for band_ind = 1:5
    band_end = band_ind*15;
    scatter(time_vecs(flight_ind,:),...
        flight_us_band_agw(band_ind,:,flight_ind), 200,'+')
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " AGW-filt")];
    hold on

    leg_ind = leg_ind+1;

    flight6_us_band(band_ind,:) = mean(...
        flight6_us_filt(band_start:band_end,:), 1, 'omitnan');

    % % Manually correcting phase shift
    % if band_ind == 4
    %     b4_phase_ind_1 = find(time_vecs(1,:) ==...
    %         "23-May-2022 17:19:47");
    %     b4_phase_ind_2 = find(time_vecs(1,:) ==...
    %         "23-May-2022 17:21:38");
    %     flight0_us_band(band_ind,b4_phase_ind_1:b4_phase_ind_2) =...
    %         flight0_us_band(4,b4_phase_ind_1:b4_phase_ind_2) +...
    %         (pi/wrap_factor); % Ok to shift by just pi instead of 2pi?
    % end
    % 
    % if band_ind == 5
    %     b5_phase_ind_1 = find(time_vecs(1,:) ==...
    %         "23-May-2022 17:21:03");
    %     b5_phase_ind_2 = find(time_vecs(1,:) ==...
    %         "23-May-2022 17:25:07");
    %     flight0_us_band(band_ind,b5_phase_ind_1:b5_phase_ind_2) =...
    %         flight0_us_band(5,b5_phase_ind_1:b5_phase_ind_2) +...
    %         (2*pi/wrap_factor);
    %     flight0_us_band(band_ind,b5_phase_ind_2:end) =...
    %         flight0_us_band(5,b5_phase_ind_2:end) +...
    %         (3*pi/wrap_factor); % Again, not a multiple of 2pi
    % end

    % % Acceleration filter of 0.07 (same as Flight 1)
    for i = 2:length(flight0_us_filt)
        if abs(flight6_us_band(band_ind, i) -...
                flight6_us_band(band_ind, i-1)) >= 0.07
            flight6_us_band(band_ind, i)  = NaN;
            flight6_us_band(band_ind, i-1)  = NaN;

        end
    end

    plot(time_vecs(flight_ind,:),...
        flight6_us_band(band_ind,:),'LineWidth',3)
    legendtxt{leg_ind} = ...
        [append(num2str(beamdist{band_start,3}), "-", ...
        num2str(beamdist{band_end,3}), " unwrap/acc ")];
    hold on

    leg_ind = leg_ind+1;
    band_start = band_end+1;
end
hold off
grid on
xlabel("Time (HH:mm)")
ylabel("Velocity (m/s)")
leg = legend(legendtxt,'Location','southoutside','Orientation','horizontal','NumColumns',5);
title(leg, "Vertical profile band h from bottom (m)");
set(gcf,'position',[10,10,1200,700])
title(append(flight_names(flight_ind), " u-velocity time series "))
fontsize(16, "points")

usurf_index = zeros(5,1);
usurf_index(1) = mean(flight0_us_band,'all','omitnan')/.85;
usurf_index(2) = mean(flight1_us_band,'all','omitnan')/.85;
usurf_index(3) = mean(flight4_us_band,'all','omitnan')/.85;
usurf_index(4) = mean(flight5_us_band,'all','omitnan')/.85;
usurf_index(5) = mean(flight6_us_band,'all','omitnan')/.85;

usurf_bin = zeros(5,1);
for flight_ind = 1:5
    binmeans = mean(flight_us_agw(:,:,flight_ind),2,'omitnan');
    B = ~isnan(binmeans);
    B_ind = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(mean(flight_us_agw,2), 2))
    usurf_bin(flight_ind) = mean(flight_us_agw(B_ind,:,flight_ind),2,'omitnan')
end

varnames = ["Flight", "Surface u-bar (index method 1/0.85) [m/s]", "Surface u-bar (top AGW bin) [m/s]", "Mean correlation %", "Mean amplitude % (normalized)"]
T = table(flight_names', usurf_index, usurf_bin, squeeze(flight_corrs_mean), squeeze(flight_amps_mean), 'VariableNames',varnames)