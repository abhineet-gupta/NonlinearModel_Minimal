function aircraftdata = AircraftData(Components)
%% AircrafttData: Generates aircaftdata.mat
% Manual code written to generate aircrft data matrix

%% Define savelocation
saveloc = Components.path_data;
filename = Components.AircraftData_file;   
%% Load FEM file
load([Components.path_data,Components.FEM_file])

%% Generate aircraft data
aircraftdata.m_AC = modes.m_AC;
aircraftdata.I_B = modes.I_B;
aircraftdata.cg = modes.xCG;
aircraftdata.s = 1.08;
aircraftdata.c = 0.4;
aircraftdata.b = 3.04;

save([saveloc,filename],'aircraftdata'); 
end