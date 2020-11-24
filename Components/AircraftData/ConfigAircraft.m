function AircraftConfig = ConfigAircraft(AircraftID)
%ConfigAircraft  Loads properties of the aircraft
%
%   This script loads the information about the structural model data,
%   gridding and spline script location and dimensional properties of the
%   aircraft.
%
%   Input:
%       AircraftID: Name of the aircraft
%
%   Output: 
%       AircraftConfig: A structure containing the following fields
%           - StrucModelFile: .mat file containging the stuctural model
%             details
%           - GriddingScript: Function handle of the script used to grid
%             the aircraft
%           - SplineScript: Function handle of the script used to obtain
%             the spline grid for the aircraft

%% Obtain location of structural data, gridding script and 
switch AircraftID
    case 'Geri'
        addpath('ConfigurationFiles_Geri')
        
        StrucModelFile = 'Geri_strucmodel';     % Structural model details
        GriddingScript = @Geri_GridData;        % Gridding function
        SplineScript = @Geri_SplineGrid;        % Spline grid function
        
        MAC = 0.3937;   % Mean aerodynamic chord
        % Area = 1.085;   % Area (m^2)
        % Span = 3.048;   % Total wing span
    
    otherwise
        error('Aircraft configuration not found')
end

%% Output Data
AircraftConfig.StrucModelFile = StrucModelFile;
AircraftConfig.GriddingScript = GriddingScript;
AircraftConfig.SplineScript = SplineScript;
AircraftConfig.MAC = MAC;

end