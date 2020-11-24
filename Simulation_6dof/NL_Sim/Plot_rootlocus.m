%% Compare the poles of the Nominal vs the PID updated model

%% Clear workspace
clear
% close all
clc

%% Load the pole zero information
% load Simulink_VLM_noAICeta.mat
load Simulink_VLM_AICeta.mat

%% Plot Simulink models
% % numc = size(simulink_6dofVLM_noAICeta.models,5)+1;
numc = size(simulink_6dofVLM_AICeta.models,5)+1;
cb=jet(numc);
figure
hold on

% pzplot(simulink_6dofVLM_noAICeta.models,'o');
pzplot(simulink_6dofVLM_AICeta.models,'o');
% pzmap(models_simulink,'k');

h = findobj(gca,'type','line'); 
for jj=2:length(h)
        set(h(jj),'MarkerSize',6,'Color',cb(numc-floor(jj/2),:));
end
grid;
title('Flutter: 29 m/s, 5.825 Hz')

colormap(jet);   % I am using parula, you can use jet, but I discourage 
c=colorbar;
caxis([20,40])

xlim([-35,10])
ylim([0,150])
