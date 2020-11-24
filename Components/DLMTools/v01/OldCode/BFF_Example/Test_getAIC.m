function [AIC] = Test_getAIC(griddata)

%% Obtian grid information
PanelData = griddata.PanelData;
NodeData = griddata.NodeData;
PanelArea = griddata.PanelArea;
n_hat_w = griddata.n_hat_wingbody;
n_hat_wl = griddata.n_hat_wl;

% Determine number of aerodynamic panels
[NumPanels,~]=size(PanelData);

%% Define downwash location (aka collocation points)
% 3/4 chord and half span of the aero panel

P0 = zeros(NumPanels,3); 
% i-th row contains x,y,z coordinate of the collocation point of i-th Panel
P0(:,1) = (NodeData(PanelData(:,2),2)+NodeData(PanelData(:,3),2) + ...
    3*(NodeData(PanelData(:,4),2)+NodeData(PanelData(:,5),2)))/8; %x
P0(:,2) = (NodeData(PanelData(:,2),3)+NodeData(PanelData(:,3),3))/2; %y
P0(:,3) = (NodeData(PanelData(:,2),4)+NodeData(PanelData(:,3),4))/2; %z

%% define doublet locations (1/4 chord, 0 and full span of the
% aero panel)
% P1 is doublet point at one of the zero span
P1 = zeros(NumPanels,3);
P1(:,1) = NodeData(PanelData(:,2),2)+ ...
    (NodeData(PanelData(:,4),2) - NodeData(PanelData(:,2),2))/4; %xp1
P1(:,2) = NodeData(PanelData(:,2),3); %yp1
P1(:,3) = NodeData(PanelData(:,2),4); %zp1


% P3 is doublet point at end of span of the panel
P3 = zeros(NumPanels,3);
P3(:,1) = NodeData(PanelData(:,3),2)+ ...
    (NodeData(PanelData(:,5),2) - NodeData(PanelData(:,3),2))/4; %xp3
P3(:,2) = NodeData(PanelData(:,3),3); %yp3
P3(:,3) = NodeData(PanelData(:,3),4); %zp3


%% Obtain AIC matrix
Dv = getVLM(P0,P1,P3,PanelArea,n_hat_w,n_hat_wl);    % VLM (steady state effects)
AIC=inv(Dv);

return