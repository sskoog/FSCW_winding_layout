% Run batch of Cros' method to find fundamental winding factor and optimal winding layout for
% pole/slot combinations in Fractional-Slot, non-overlapping
% Concentrated-Winding Permanent Magnet Machines

clear all
close all
clc

% Set up vector of machines to analyze
cros( 24,32,1,1 ); % Example machine 1
cros( 18,20,1,1 ); % Example machine 2


p = 2:2:80;
Q = 3:3:90;
% Make large batch of FSCW designs and evaluate winding factor
for i=1:length(p)
    for k=1:length(Q)
        kw(i,k) = cros( Q(k),p(i),0,0 );        
    end
end

%% Plot results nicely in 3D graphs
f = figure();
b = bar3(kw);
ax1 = b.Parent;
% Apply gradient colors based on height
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colorbar;

ax1.XTick = 1:2:length(Q);
ax1.XTickLabels = Q(1:2:end);
ax1.YTick = 1:2:length(p);
ax1.YTickLabels = p(1:2:end);
ax1.YDir = 'normal';
ylabel('p = Poles');
xlabel('Q_s = Slots');
title('Fundamental winding factor');

% Plot only results over 85% winding ratio
f = figure();
b = bar3(kw);
ax1 = b.Parent;
ax1.ZLim = [0.85 1];
ax1.XTick = 1:2:length(Q);
ax1.XTickLabels = Q(1:2:end);
ax1.YTick = 1:2:length(p);
ax1.YTickLabels = p(1:2:end);
ax1.YDir = 'normal';
% Apply gradient colors based on height
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
caxis([0.85 0.96]);
colorbar('Limits',[0.85 1],'Ticks',[0.85:0.02:1.0]);
ylabel('p = Poles');
xlabel('Q_s = Slots');
title('Fundamental winding factor >0.85');
% EOF
