function [AircraftData] = MakeAircraft(PreMultiplyingMat,PostMultiplyingMat)
%MakeAircraft  Obtain aerodynamic model and generate data for simulation
%
%   Loads FEM data, runs VLM and DLM code and saves the information needed
%   for simulation with 6DOF nonlinear flight dynamics mechanics, 6DOF 
%   linear FEM based structural model and 6DOF geometrically non-linear 
%   aerodynamics.

%% Setup
close all

% Select aircraft
AircraftID = 'Geri';

% Select if figuers are displayed
DisplayFig = 0;

% Select aerodynamic type
%   Aeroypte = 'VLM' for only steady aerodynamic model
%   Aeroypte = 'DLM' for steady and unsteady aerodynamic model
AeroType = 'VLM';

% Display aerodynamic model choice
% if strcmp(AeroType,'VLM')
%     fprintf('Steady aerodynamic model chosen\n')
% elseif strcmp(AeroType,'DLM')
%     fprintf('Unsteady aerodynamic model chosen\n')
%     addpath('DLMCode/')
% else
%     error('AeroTYpe should be "VLM" or "DLM"')
% end

%% Aircraft configuration
% Configure the aircraft
AircraftConfig = ConfigAircraft(AircraftID);

% Data used from aircraft configutaion 
    % Mean Aerodynamic Chord (reference length)
    MAC = AircraftConfig.MAC;
    % Structural Model details
    StrucModelFile = AircraftConfig.StrucModelFile; 
    % Script to obtain aerodynamic grid
    GriddingScript = AircraftConfig.GriddingScript; 
    % Script to obtain spline grid
    SplineScript = AircraftConfig.SplineScript;

%% Structural model
% Load the structural model
load(StrucModelFile);

% Data used from structural model
    % Coordinates of FEM nodcdes 
    %   (See FEM code for details about the coordinates system)
    Coords_FEM = strucgrid.coords;
    % Position vector of the structural node corresponing to the CG
    Rcg = strucmodes.rcg;                        
    % Number of sensors
    NumSensors = strucgrid.numsensors; 
    % Sensor coordinates    
    SensorLoc = strucgrid.coords_sensors;
    
%% Obtain aerodynamic grid
% Run gridding script
griddata = feval(GriddingScript);

    % Data used from the grid
    NumPanels = griddata.NumPanels;         % Number of aerodynamic panell
    Aerogrid = griddata.Aerogrid;           % Centroid of aerodynamic panels
    BodySpan = griddata.BodySpan;           % Span of the center body
    PanelSpan = griddata.PanelSpan;         % Span of each panel
    PanelLength = griddata.PanelLength;     % Flowwise length of each panel
    PanelData = griddata.PanelData;         % PanelData
    NodeData = griddata.NodeData;           % NodeData
    Nhat0 = griddata.Nhat0;                 % Normal vector of each panel
    CSHingeDC = griddata.CSHingeDC;         % CS hinge direction cosine
    CSPanelIndex = griddata.CSPanelIndex;   % Index of each CS panels
    CSPanelIndex_all = griddata.CSPanelIndex_all;   % Index of all CS panels

% Find nodes in aerodynamic grid corresponding to sensor locations
SensorNodes = zeros(1,NumSensors);
for i = 1:NumSensors
    [~,SensorNodes(i)] = min(vecnorm((Aerogrid-SensorLoc(:,i))));
end

% If correction matrices are not provided, assume them to be identity
switch nargin
    case 0
        PreMultiplyingMat = eye(NumPanels);
        PostMultiplyingMat = eye(NumPanels);
        disp('Identity matrices chosen as pre and post correction matrices')
    case 1
        error('Provided none or both the correction matrices, [] can be used')
    case 2
        if isempty(PreMultiplyingMat)
            PreMultiplyingMat = eye(NumPanels);
        elseif (size(PreMultiplyingMat,1) ~= NumPanels) || ...
                (size(PreMultiplyingMat,2) ~= NumPanels)
            error('PreCorrectionMat should be a square matrix with size = NumPanels-by-NumPanels')
        end
        if isempty(PostMultiplyingMat)
            PostMultiplyingMat = eye(NumPanels);
        elseif (size(PostMultiplyingMat,1) ~= NumPanels) || ...
                (size(PostMultiplyingMat,2) ~= NumPanels)
            error('PostCorrectionMat should be a square matrix with size = NumPanels-by-NumPanels')
        end
    otherwise
        error('Check pre and post correction matrices')
end


%% Obtain spline grid
[Coords_spline,Connections] = SplineScript(Coords_FEM,BodySpan); 

%% Transformation matrices
% Structural to spline 
Tsf = calc_Tsf(Coords_FEM,Coords_spline,Connections);   
% Spline to aerodynamic
Tas = calc_Tas(Coords_spline,Aerogrid);
Tas_node = calc_Tas(Coords_spline,NodeData);
% Structural to aerodynamic
Taf = Tas*Tsf;
Taf_node = Tas_node*Tsf;

%% Obtain flexible mode shapes
% Flexible mode shapes in aerodynamic grid
Phif_aero = Taf*strucmodes.Phif;     
Phif_aero_node = Taf_node*strucmodes.Phif;     

% Number of flexible modes
NumfModes = size(Phif_aero,2);  

% Display flexible modes
% NOTE: The flexible modes are scaled for better visualization. This
% scaling does not affect the calculated modes
if DisplayFig
    flexfigscale = 0;
    for i = 1:size(Phif_aero,2)
        figure
        for pan = 1:NumPanels
        x = [NodeData(1,PanelData(1,pan)), NodeData(1,PanelData(2,pan));
            NodeData(1,PanelData(3,pan)), NodeData(1,PanelData(4,pan))] ...
            + flexfigscale*[Phif_aero_node(6*PanelData(1,pan)-5,i), Phif_aero_node(6*PanelData(2,pan)-5,i);
            Phif_aero_node(6*PanelData(3,pan)-5,i), Phif_aero_node(6*PanelData(4,pan)-5,i)];
        y = [NodeData(2,PanelData(1,pan)), NodeData(2,PanelData(2,pan));
            NodeData(2,PanelData(3,pan)), NodeData(2,PanelData(4,pan))] ...
            + flexfigscale*[Phif_aero_node(6*PanelData(1,pan)-4,i), Phif_aero_node(6*PanelData(2,pan)-4,i);
            Phif_aero_node(6*PanelData(3,pan)-4,i), Phif_aero_node(6*PanelData(4,pan)-4,i)];
        z = [NodeData(3,PanelData(1,pan)), NodeData(3,PanelData(2,pan));
            NodeData(3,PanelData(3,pan)), NodeData(3,PanelData(4,pan))] ...
            + flexfigscale*[Phif_aero_node(6*PanelData(1,pan)-3,i), Phif_aero_node(6*PanelData(2,pan)-3,i);
            Phif_aero_node(6*PanelData(3,pan)-3,i), Phif_aero_node(6*PanelData(4,pan)-3,i)];
        surf(x,y,z)
        hold on
        end
        axis equal
        view(-60,5)
        title(['Mode Shape: Mode', num2str(i)])
%         xlabel('x'); ylabel('y'); zlabel('z');
        grid on;
    end
end

if DisplayFig
    flexfigscale = 4;
    for i = 1:size(Phif_aero,2)
        figure
        plot3(Aerogrid(1,:),Aerogrid(2,:),Aerogrid(3,:),'b*')
        hold on
        plot3(Aerogrid(1,j),Aerogrid(2,j),Aerogrid(3,j),'r+')
        plot3(Aerogrid(1,k),Aerogrid(2,k),Aerogrid(3,k),'r+')
%         plot3(Aerogrid(1,:)+flexfigscale*Phif_aero(1:6:end,i)',...
%                 Aerogrid(2,:)+flexfigscale*Phif_aero(2:6:end,i)', ...
%                 Aerogrid(3,:)+flexfigscale*Phif_aero(3:6:end,i)','r*')
        xlabel('x'); ylabel('y'); zlabel('z'); grid on;
        legend('Undeflected nodes','Deflected nodes')
        axis equal
        view(110,60)
        title(['Aerodynamic grid flexible mode: ', num2str(i)])
    end
end



%             
%% Obtain panel direction dependency on flexible deflection
Nhat_eta = getNhateta(Nhat0,Phif_aero,NumPanels,NumfModes);

%% Obtain steady aerodynamic model
VLMData = getSteadyAero(NumPanels,PanelData,NodeData,Nhat0,Nhat_eta,Rcg);
D0 = VLMData.D0;
% warning('Not all aerodynamic quantities are corrected as of now')
VLMData.Ainv0 = PreMultiplyingMat*VLMData.Ainv0*PostMultiplyingMat;

%% If needed, obtain the unsteady aerodynamic model
if strcmp(AeroType,'DLM')
    DLMData = getUnsteadyAero(NumPanels,PanelData,NodeData, ...
                                       PanelSpan,PanelLength,MAC,D0);
end

%% Save AircraftData for Simulation
AircraftData.Aero.VLMData = VLMData;
AircraftData.Aero.NumPanels = NumPanels;
AircraftData.Aero.Aerogrid = Aerogrid;
AircraftData.Aero.Nhat0 = Nhat0;
AircraftData.Aero.Nhat_eta = Nhat_eta;
AircraftData.Aero.MAC = MAC;
if strcmp(AeroType,'DLM')
    AircraftData.Aero.DLMData = DLMData;
end

AircraftData.CS.CSTau = 300;
AircraftData.CS.CSHingeDC = CSHingeDC;
AircraftData.CS.CSPanelIndex = CSPanelIndex;
AircraftData.CS.CSPanelIndex_all = CSPanelIndex_all;

AircraftData.SensorData.SensorNodes = SensorNodes;
AircraftData.SensorData.NumSensors = NumSensors;
AircraftData.SensorData.SensorAeroGrid = Aerogrid(:,SensorNodes);
AircraftData.SensorData.SensorCG = Aerogrid(:,SensorNodes)-Rcg;

AircraftData.MassProp.mass = strucmodes.Mass;
AircraftData.MassProp.inertia = strucmodes.Inertia;

AircraftData.ModalStruc.InvMass = 1./strucmodes.Mf;
AircraftData.ModalStruc.Damping = strucmodes.Bf;
AircraftData.ModalStruc.Stiffness = strucmodes.Kf;
AircraftData.ModalStruc.Phif_aero = Phif_aero;

save(['../../Simulation_6dof/Configuration/AircraftData_Updated_',AircraftID],'AircraftData')

end