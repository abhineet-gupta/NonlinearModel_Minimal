function aircraftdata = AircraftData(modes)
%AircraftData   Stores mass properties and dimensions of the aircraft
%   This script assembles and stores aircraft mass properties and
%   dimensions for further use. All data is stored in SI units
% 
%   Input:
%       - modes: A structure with following fields (might contain other
%           fields not used here) 
%           - Mass: Mass of the aircraft in kg
%           - Inertia: Inertia matrix of the aircraft in kg-m^2
%           - xcg: Distance of c.g. of the aircraft from the nose in m% 
% 
%     Output: 
%       - aircraftdata: A structure containing
%           - m_AC: Mass of the aircraft (Kg)
%           - I_B: Inertia matrix of the aircraft (Kg-m^2)
%           - cg: Distance of CG of the aircraft from nose (m)
%           - s: Planform area of the aircraft (m^2)
%           - c: Mean aerodynamic chord (m)
%           - b: Total wingspan (m)

%% Generate aircraft data
aircraftdata.m_AC = modes.Mass;           % Mass of the aircraft
aircraftdata.I_B = modes.Inertia;         % Inertia matrix
aircraftdata.cg = modes.rcg(1);                % CG location from nose
aircraftdata.s = 1.085;                     % Area (m^2)
aircraftdata.c = 0.3937;                    % Mean aerodynamic chord
aircraftdata.b = 3.048;                     % Total wing span

save('Data_files/AircraftData','aircraftdata'); 