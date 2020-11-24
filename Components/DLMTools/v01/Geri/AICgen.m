%AICgen   Obtain AIC matrix using VLM or DLM
%   This script uses VLM/DLM model to obtain the aerodynamic influence
%   coefficient matrix. It further obtaines the generalized aerodynamic
%   matrices based on the rigid body, flexible and control surface modes. 
%   The script saves the relevant data in a mat file to be used by the
%   SIMULINK model

%% Clear workspace and setup
clear;
close all;
clc;

% Choose whether 'VLM' or 'DLM' model is used
AeroType = 'DLM';   % Select AeroType = 'VLM' or 'DLM'

%% Identify project path
addpath('../')
addpath('../GAM_generation')
addpath('./Gridding')
addpath('./GAM_matrices')

%% Define FEM model
load('Data_files/Geri_strucmodel.mat');

%% Obtain grid information
griddata = GridData_Geri(modes.rcg);
NumPanels = griddata.NumPanels;

%% Obtain aircraft data
aircraftdata = AircraftData(modes);
cref = aircraftdata.c;

%% Obtain differentiation and transformation matrices
[geom,modes] = GAM_matrices(griddata,cref,geom,modes);    

%% get AIC matrices for 8 red freq
k = [0 0.0469 0.0938 0.1875 0.375 0.75 1.5 3];  % freq. considered by DLR
kr = (2/cref)*k;
Mach = 0;

AIC = zeros(NumPanels,NumPanels,length(k));
for i = 1:length(k)
    AIC(:,:,i) = getAIC(griddata,kr(i),AeroType);
    fprintf('Cycle completed with kr = %4.3f\n',kr(i))
end

%% get AIC in modal coordinates 
% % Phibf = [modes.Phib modes.Phif];
[genAIC,A0,A1] = modal_aero(AIC,griddata.Phib_aero,modes.Phif,...
                    griddata.PhiCS_aero,geom.Tfa,griddata.S,...
                    geom.D1,geom.D2,k);

%% Rational Fraction Approximation
b_poles = [0.1100 0.2200];    % lag poles

if strcmp(AeroType,'DLM')
    [A_0,A_1,A_2,A_p,A_lag,B_lag,C_lag] = ...
                    rational_approx(genAIC,b_poles,k);
    aero.RFA_Qhhx_nln.A0 = A_0;
    aero.RFA_Qhhx_nln.A1 = A_1;
    aero.RFA_Qhhx_nln.A2 = A_2;
    aero.RFA_Qhhx_nln.Ap = A_p;
    aero.RFA_Qhhx_nln.Alag = A_lag;
    aero.RFA_Qhhx_nln.Blag = B_lag;
    aero.RFA_Qhhx_nln.Clag = C_lag;
    aero.RFA_Qhhx_nln.poles = b_poles;
    
elseif strcmp(AeroType,'VLM')
    NumBFModes = 12;
    NumCSModes = 8;
    
    aero.RFA_Qhhx_nln.A0 = A0;
    aero.RFA_Qhhx_nln.A1 = A1;
    aero.RFA_Qhhx_nln.A2 = zeros(NumBFModes,NumBFModes + NumCSModes);
    aero.RFA_Qhhx_nln.Alag = ...
                    blkdiag(-b_poles(1)*eye(NumBFModes + NumCSModes), ...
                            -b_poles(2)*eye(NumBFModes + NumCSModes));
    aero.RFA_Qhhx_nln.Blag = ...
                    [eye(NumBFModes + NumCSModes);...
                     eye(NumBFModes + NumCSModes)];
    aero.RFA_Qhhx_nln.Clag = zeros(NumBFModes,2*(NumBFModes + NumCSModes));
    aero.RFA_Qhhx_nln.poles = b_poles;
end
%% Output structure to be passed to simulation
miniMUTT_dyn.aero = aero;
miniMUTT_dyn.geom = geom;
miniMUTT_dyn.modes = modes;

%% Save data in relevant .mat file 
save ('Geri_Aero_VLM.mat','miniMUTT_dyn');