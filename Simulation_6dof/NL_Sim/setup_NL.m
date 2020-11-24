%% setup_NL  Setup the data required for nonlinear simulation
function linmodel = setup_NL(vi)
%% Clear workspace
% clear
% % close all
% bdclose all
% clc
%% Simulation Parameters
% Select Aircraft
addpath ../Configuration
addpath ../Libraries
% load AircraftData_Geri.mat
load AircraftData_Updated_Geri.mat

AircraftID = 'Geri';

% Select aerodynamic model type ('DLM needs to be updated')
AeroType = 'VLM';

if ~strcmp(AeroType,'VLM') && ~strcmp(AeroType,'DLM')
    error('AeroType should be "VLM" or "DLM"')
end

% Select if VLM AIC depends on eta's on not (AICdueEta = true or false)
AICdueEta = false;

%% Add library path and model name
% addpath('../Libraries')
mdl = 'NL_Simulation';

%% Configure aircraft
% addpath('../Configuration/')
[AC,Env] = ConfigureUAV(AircraftID,AircraftData);

% Index of all CS panels
Idx_PanelIndexAll  = AC.CS.CSPanelIndex_all;

%% Initialize the model
Vinf = vi;      % Aircraft velocity for initialization
StateVals = SetInitial(Vinf,AeroType);

%% Load the model
load_system(mdl)
mdlWks = get_param(mdl,'ModelWorkspace');

assignin(mdlWks,'AC',AC);
assignin(mdlWks,'Env',Env);
assignin(mdlWks,'StateVals',StateVals);

%% Comment out the unsteady aerodynamics blocks if AeroType = VLM
if strcmp(AeroType,'VLM')
    set_param(['NL_Simulation/Force and Moment /Aerodynamic Model/',...
        'Unsteady Aerodynamics'],'commented','on')
elseif strcmp(AeroType,'DLM')
    set_param(['NL_Simulation/Force and Moment /Aerodynamic Model/',...
        'Unsteady Aerodynamics'],'commented','off')
end


%% Comment out AIC (VLM) dependency on eta, if required
if AICdueEta
    set_param(['NL_Simulation/Force and Moment /',...
        'Aerodynamic Model/Steady Aerodynamics/Vortex Strength/',...
        'AIC Calculation/AIC due to Flex'],'commented','off')
else
    set_param(['NL_Simulation/Force and Moment /',...
        'Aerodynamic Model/Steady Aerodynamics/Vortex Strength/',...
        'AIC Calculation/AIC due to Flex'],'commented','on')
end

%% Trim the model
set_param(mdl,'LoadExternalInput','off')
set_param(mdl,'LoadInitialState','off')
% mdl = 'NL_Simulation';
% 
% % Create the operating point specification object.
% opspec = operspec(mdl);

Op_Trim = ObtainTrimPoint(StateVals,AeroType,AC); 

% %% Linearize the model
linmodel = LinearizeModel(mdl,Op_Trim);
end