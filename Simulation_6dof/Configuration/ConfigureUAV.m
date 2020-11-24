function [AC,Env] = ConfigureUAV(AircraftID,AircraftData)
%% ConfigureUAV Configues and pre-process data of aircraft for simulation
%
%   Input:
%       - AicraftID: Name of the aircraft
%
%   Output:
%       - AC: Aircraft specific data passed to the simulation
%       - Env: Environment specefic data passed to the simulation

%%
switch AircraftID
    case 'Geri'
%         load('AircraftData_Geri','AircraftData')
        AC =  AircraftData;
        
        %% Environment details
        Env.GroundAlt = 0;
        
        %% Postprocessing
        % Information about the aircraft is pre-processed to ensure a
        % faster simulation
        
        % Number of panels
        NumPanels = AC.Aero.NumPanels;
        
        % Index of panels belonging to each control surface
        AC.CS.CSPanelIndex1 = AC.CS.CSPanelIndex{1};
        AC.CS.CSPanelIndex2 = AC.CS.CSPanelIndex{2};
        AC.CS.CSPanelIndex3 = AC.CS.CSPanelIndex{3};
        AC.CS.CSPanelIndex4 = AC.CS.CSPanelIndex{4};
        AC.CS.CSPanelIndex5 = AC.CS.CSPanelIndex{5};
        AC.CS.CSPanelIndex6 = AC.CS.CSPanelIndex{6};
        AC.CS.CSPanelIndex7 = AC.CS.CSPanelIndex{7};
        AC.CS.CSPanelIndex8 = AC.CS.CSPanelIndex{8};
        
        % Direction cosines of the control surface hinges
        AC.CS.CSHingeDC1 = AC.CS.CSHingeDC(:,1);
        AC.CS.CSHingeDC2 = AC.CS.CSHingeDC(:,2);
        AC.CS.CSHingeDC3 = AC.CS.CSHingeDC(:,3);
        AC.CS.CSHingeDC4 = AC.CS.CSHingeDC(:,4);
        AC.CS.CSHingeDC5 = AC.CS.CSHingeDC(:,5);
        AC.CS.CSHingeDC6 = AC.CS.CSHingeDC(:,6);
        AC.CS.CSHingeDC7 = AC.CS.CSHingeDC(:,7);
        AC.CS.CSHingeDC8 = AC.CS.CSHingeDC(:,8);
        
        % Dependence of panel normal vectors on 
        Nhat0 = AC.Aero.Nhat0;          
        Nhat_eta = AC.Aero.Nhat_eta;
        NumfModes = size(Nhat_eta,3);

        Nhat_eta_all = zeros(3*NumPanels,NumfModes+1);
        Nhat_eta_all(:,1) = reshape(Nhat0,3*NumPanels,1);
        for i = 1:NumfModes
            Nhat_eta_all(:,i+1) = reshape(Nhat_eta(:,:,i), ...
                3*NumPanels,1);
        end
        AC.Aero.Nhat_eta_all = Nhat_eta_all;
        
        % Obtain translation part from flexible moode shapes
        idx_trans = sort([1:6:6*NumPanels, 2:6:6*NumPanels, ...
            3:6:6*NumPanels]);
        AC.ModalStruc.Phif_aero_trans = ...
            AC.ModalStruc.Phif_aero(idx_trans,:);
        
        % Obtain mode flexible information at sensor locations
        idx_sensor_Trans = [];
        idx_sensor_Rot = [];
        for i = 1:AC.SensorData.NumSensors
            idx_sensor_Trans = [idx_sensor_Trans, ...
                AC.SensorData.SensorNodes(i)*6-5, ...
                AC.SensorData.SensorNodes(i)*6-4, ...
                AC.SensorData.SensorNodes(i)*6-3];
            idx_sensor_Rot = [idx_sensor_Rot, ...
                AC.SensorData.SensorNodes(i)*6-2, ...
                AC.SensorData.SensorNodes(i)*6-1, ...
                AC.SensorData.SensorNodes(i)*6];
        end
        
        AC.ModalStruc.Phif_aero_sensor_Trans = ...
            AC.ModalStruc.Phif_aero(idx_sensor_Trans,:);
        
        AC.ModalStruc.Phif_aero_sensor_Rot = ...
            AC.ModalStruc.Phif_aero(idx_sensor_Rot,:);
        
        AC.SensorData.SensorOrientation0 = ...
            Nhat0(:,AC.SensorData.SensorNodes);
        
        % Obtain actuator model of the aircraft
        % The actuator is modeled by a second order The transfer function 
        % which can be written as: b/(s^2 + as +b) where:
        AC.Actuator.a = 840;
        AC.Actuator.b = 96710;
        
end
