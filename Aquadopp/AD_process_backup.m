% Spectra??

% d1rms = sqrt(beam1d1.^2 + beam2d1.^2 + beam3d1.^2);
% d1rms(:,1:2) = [];
% d2rms = sqrt(beam1d2.^2 + beam2d2.^2 + beam3d2.^2);
% d2rms(:,1:2) = [];

% N=[length(beam1d1(:,1)) length(beam1d2(:,1))];
% bin_num = length(beam1d1(1,:));

% for i=1:N(1)
%     i
%     for bin=3:bin_num
%         Raa1(i,bin)=0;
%         Rcc1(i,bin)=0;
%         Rac1(i,bin)=0;
%         for t=1:N(1)
%             k=t+(i-1);
%             if k<=N(1)
%                 k=k;
%             else
%                 k=k-N(1);
%             end
%             Raa1(i,bin)=Raa1(i,bin)+beam1d1(t,bin)*beam1d1(k,bin);
%             Rcc1(i,bin)=Rcc1(i,bin)+beam3d1(t,bin)*beam3d1(k,bin);
%             Rac1(i,bin)=Rac1(i,bin)+beam1d1(t,bin)*beam3d1(k,bin);
%         end
%     end
% end
% Raa1=Raa1./N(1);
% Rcc1=Rcc1./N(1);
% Rac1=Rac1./N(1);