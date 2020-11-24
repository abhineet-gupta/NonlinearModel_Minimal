function [AIC] = getAIC(PanelData,NodeData,k,PanelLength,PanelSpan,PanelAreas,n_hat_w,n_hat_wl)
%
% This function computes the AIC Matrix using ref. 1 as the primary
% reference for doublet lattice (DLM) implementation and ref 2 for vortex lattice   
% (VLM) implementation. The AIC matrix is a combination of DLM and VLM
%
% ref 1: Albano and Rodden - A Doublet-Lattic Method for Calculating 
%        Lift Distributions on Oscillating Surfaces in Subsonic Flows
% 
% ref 2: Katz & Plodkin, 'Low speed aerodynamics', second edition (for VLM
%        implementation)
%
% Inputs:
% 1. PanelData = 5 column array defining [panel number, southwest node #, northwest node #, southeast node #, northeast node #]
%             the node number is the same number as defined in the first column of NodeData
%
% 2. NodeData = 4 column array defining [node #, x coordinate, y coordinate, z coordinate]
%
% 3. Mach = Mach number
%
% 4. k = k1 term from ref 1.: omega/U
%    omega = frequency of oscillation
%    U = freestream velocity
% 
% 5. S: panel areas
%
% 6. n_hat_w: normal vector information for panels on the wing
% 7. n_hat_wl: normal vector information for panels on winglets
% 
% n_hat_w and n_hat_wl indicate the direction of the normal vector for all
% panels. n_hat_w contains cos(gamma) where 'gamma' is the dihedral angle
% of a panel. n_hat_wl contains sin(gamma) of the panels. This information
% is useful in the VLM code
%
% Outputs
% AIC = aerodynamic influence coefficient matrix
% rows of D are the downwash locations while columns define the effects of
% each aero panels doublet line/horseshoe vortex on that rows' downwash location
%
% Code written by
% Aditya Kotikalpudi
% Graduate Assistant
% University of Minnesota
%
%
% Code based on:
% Original code by Frank R. Chavez, Iowa State University (~2002)
% Modified version by:
% Brian P. Danowsky
% P. Chase Schulze
% (c) Systems Technology, Inc. 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of aero panels present
[N,~]=size(PanelData);

% define downwash location (3/4 chord and half span of the aero panel)
P0 = zeros(N,3);
P0(:,1) = (NodeData(PanelData(:,1),1)+NodeData(PanelData(:,2),1) + ...
    3*(NodeData(PanelData(:,3),1)+NodeData(PanelData(:,4),1)))/8; %xcp
P0(:,2) = (NodeData(PanelData(:,1),2)+NodeData(PanelData(:,2),2))/2; %ycp
P0(:,3) = (NodeData(PanelData(:,1),3)+NodeData(PanelData(:,2),3))/2; %zcp


% define doublet locations (1/4 chord and 0, half span and full span of the
% aero panel), kernel is computed at that 3 points and a parabolic function
% is fitted to approximate the kernel along the doublet line (ref 1,
% equation 7).

% P1 is doublet point at one of the zero span
P1 = zeros(N,3);
P1(:,1) = NodeData(PanelData(:,1),1)+ ...
    (NodeData(PanelData(:,3),1) - NodeData(PanelData(:,1),1))/4; %xp1
P1(:,2) = NodeData(PanelData(:,1),2); %yp1
P1(:,3) = NodeData(PanelData(:,1),3); %zp1


% P3 is doublet point at end of span of the panel
P3 = zeros(N,3);
P3(:,1) = NodeData(PanelData(:,2),1)+ ...
    (NodeData(PanelData(:,4),1) - NodeData(PanelData(:,2),1))/4; %xp3
P3(:,2) = NodeData(PanelData(:,2),2); %yp3
P3(:,3) = NodeData(PanelData(:,2),3); %zp3

% P2 is doublet point at the half-span location of the panel
P2 = zeros(N,3);
P2(:,1) = (P1(:,1)+P3(:,1))/2;
P2(:,2) = (P1(:,2)+P3(:,2))/2;
P2(:,3) = (P1(:,3)+P3(:,3))/2;

% define half span length and chord at centerline for each panel
% s = 0.5*((NodeData(PanelData(:,3),3)-NodeData(PanelData(:,2),3))' ...
%     + (NodeData(PanelData(:,3),4) - NodeData(PanelData(:,2),4))');
% c = ((NodeData(PanelData(:,4),2) - NodeData(PanelData(:,2),2) + ...
%     NodeData(PanelData(:,5),2) - NodeData(PanelData(:,3),2))/2)';

% get the downwash effect from DLM and VLM implementations for the defined geometry 
Dd = getDLM(P0,P1,P2,P3,PanelSpan/2,PanelLength,k,0);  % DLM (oscillatory effects)
Dv = getVLM(P0,P1,P3,PanelAreas,n_hat_w,n_hat_wl);    % VLM (steady state effects)

% rows of D are the downwash locations while columns define the effects of
% each aero panels doublet line/horseshoe vortex on that rows' downwash 
% location
D = Dv + Dd;
% D = Dv;

AIC=inv(D);

return