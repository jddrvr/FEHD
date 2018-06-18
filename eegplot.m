function handle = eegplot(data);
%EEGPLOT Plot optimized to display swatches of EEG data
%   HANDLE = EEGPLOT(DATA) accept DATA whose rows are channels of data
%   points and returns the handle to the plot.
%
%   If DATA has 16 channels, the channels are labelled according to the
%   standard EEG montage. If DATA has any other number of channels, they
%   are numbered in order.
%
%   The channel-to-electrode conversion table for the standard EEG
%   montage is:
%      1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
%      Fp1 F3  F7  C3  T3  T5  P3  O1  Fp2 F4  F8  C4  T4  T6  P4  O2
%
%   See also PLOT.
%
%   Michael A. Repucci
%   Last Modified: 9/13/2000 2:51 AM
%   Adapted from work by Matthew Davey

[m n] = size(data);
if m > n
   error('Input DATA matrix must be size=[m n] where m<n');
end

data = data-repmat(min(data,[],2),[1 n]);
maxs_minus_mins = [max(data,[],2)-min(data,[],2);0];
space = maxs_minus_mins(2:end);
for i = m-1:-1:1
   space(i) = space(i) + space(i+1);
end
plotspace = repmat(space,[1 n]);

handle = figure('numbertitle','off','name','EEG-Style Plot');
plot((plotspace+data)');
set(gca,'fontsize',16,'position',[0.08 0.14 0.89 0.83]);
set(handle,'color',[1 1 1]);

axis([1 n 0 max(data(1,:))+space(1)]);
ytickpts = space + .5*maxs_minus_mins(1:end-1);
set(gca,'ytick',flipud(ytickpts));

if m == 16,
   yticknames = [' O2';' P4';' T6';' T4';' C4';' F8';' F4';'Fp2'; ...
        ' O1';' P3';' T5';' T3';' C3';' F7';' F3';'Fp1'];   
else
   yticknames = m:-1:1;
end

set(gca,'yticklabel',yticknames);
xlabel('Time Points');
ylabel('Component');
set(gcf,'paperorientation','landscape')
set(gcf,'paperposition',[.5 .5 10 8])
zoom on; 
grid on;
