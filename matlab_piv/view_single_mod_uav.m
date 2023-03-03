%Basic m-file to display tiff images collected as individual tiff

clear
close all

s_path='C:\Users\evanh\Box\Cornell\Fall_2021\Carp_UAV\0036SET\';

Nconvert=4; %accounts for 14-bits images being stored in top 16 bits

iPairSep=1; %number of images across pair (2 = adjacent images)

% v = VideoWriter('carp_vec.avi');
% open(v);

%Loop over images and display them
therm_imgs = dir(fullfile(s_path, '**\*6.tif'));
for i=450:478 %330:478 %Must be 3 digits, from 001 to 992, range for mtg 330:510, stationary 330:478
    fname=[therm_imgs(i).folder '\' therm_imgs(i).name];
        disp(fname)
        I1=double(imrotate(imread(fname), 180)/Nconvert); %loads image into a double precision real variable matrix
%         min1=min(I1(:));
%         max1=max(I1(:));
%         mean1=mean(I1(:));
%         median1=median(I1(:));
%         disp(sprintf('Min1=%3d, Max1=%4d',min1,max1));
%         disp(sprintf('Mean1=%6.2f, Median1=%6.2f',mean1,median1));
    
    figure(1)
    set(gcf,'WindowState','fullscreen');
    imagesc(I1,[7140 7180]); %Set thermal resolution scale here
    axis image
    colormap(jet(256));
%      a = colorbar;
%      title (num2str(fname))
%     xlabel('x-distance (pixels)');
%     ylabel('y-distance (pixels)');
%     ylabel(a, 'Raw pixel intensity (#)');
%     set(gca, 'FontSize', 28);
    pause
%     hold on
%     quiver(IgXc,IgYc,Umed,-Vmed,'Color',[0 0 0],'LineWidth',3)
%     axis ij
%     axis image
    
%     pause(0.01);
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end
% close(v)