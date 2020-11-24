function [geom] = GAM_matrices(Components)
%% GAM matrices: Generates matrices
% This script constructs the matrices required to carry out transformantion
% of the aerodynamic model from aerodynamic grid to the structural grid
% in modal coordinates

saveloc = Components.path_data;
filename = Components.GAM_file;

addpath([Components.path_main, '../GAM_generation'])

% Load FEM Data
load([Components.path_data,Components.FEM_file])

% Load GridData
load([Components.path_data,Components.Grid_file])

% Load Aircraft Data
load([Components.path_data,Components.AircraftData_file])
%% calculate transformation matrices
aerogrid = griddata.aerogrid';
Npanel = griddata.num_of_panels;
Plength = griddata.PanelLength';
c_ref = aircraftdata.c;

% Create spling grid between 
[offset_s,connections] = spline_grid(geom.strcgrid.offset); 

% Calculate differentiation matrix
% Ref: Kier T.M., Looye G.H, 
% "Unifying Manoeuver and Gust Loads Analysis Models"
[D1,D2] = getDiff(Npanel,Plength,c_ref);
geom.D1 = D1;
geom.D2 = D2;

geom.c_ref = c_ref;
geom.splinegrid.offset = offset_s;
geom.splinegrid.connec = connections;

% Calculate transformation from structural grid to spline grid
geom.Tsg = calc_Tsg(geom.strcgrid.offset, ...
           geom.splinegrid.offset,geom.splinegrid.connec);

% Calculate transformation from spline to aero grid
geom.Tks = calc_Tks(geom.splinegrid.offset,aerogrid');

% Calculate transformation from structural to aero grid
geom.Tkg = geom.Tks*geom.Tsg;

geom.aerogrid = aerogrid';
geom.phiCS = griddata.phiCS;
geom.S = griddata.S;

save([saveloc,filename],'geom');
end