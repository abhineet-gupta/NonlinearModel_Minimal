%ExtractFEM_Geri   Extracts information from FEM model
%
%   Extracts modal information from the FEM model 5.2 of mAEWing1-Geri.
%   The coordinate frame chosen for the final FEM results is defined such 
%   that:
%       A) Origin at the vehicle nose, 
%       B) +x pointing to rear or aircraft, and 
%       C) +y pointing starboard (to the right when looking forward)
%       D) +z pointing upward 
%   The processed data is stored in the a mat file with the following
%   variables
%
%   strucmodes: A structure with following fields
%       - Omegaf: A numflexmodes-by-1 vector containing frequencies of the
%           flexible modes in rad/sec
%       - Zetaf: A numflexmodes-by-1 vector containing damping ratio of the
%           flexible modes
%       - Phif: A (6*numnodes)-by-numflexmodes matrix where the columns
%           contain the flexible modes modes where each mode contains 
%           [surge,sway,heave,roll,pitch,yaw] degrees of freedom of each 
%           FEM node in the final FEM cooridante described above
%       - Mf: A numflexmodes-by-1 vector containing the modal masses of
%           the flexible modes (SI units are used)
%       - Bf: A numflexmodes-by-1 vector containing the modal damping of
%           the flexible modes (SI units are used)
%       - Kf: A numflexmodes-by-1 vector containing the modal stiffness of
%           the flexible modes (SI units are used)
%       - M: A (6*numnodes)-by-(6*numnodes) matrix containing the
%           mass matrix from the FEM with 
%           [surge,sway,heave,roll,pitch,yaw] as degrees of freedom 
%           (SI units are used)
%       - K: A (6*numnodes)-by-(6*numnodes) matrix containing the
%           mass matrix from the FEM with 
%           [surge,sway,heave,roll,pitch,yaw] as degrees of freedom (SI 
%           units are used)
%       - Mass: Mass of the aircraft in kg
%       - Inertia: Inertia matrix of the aircraft in kg-m^2
%       - rcg: A 3-by-1 vector containing [x;y;z] coordinates of the c.g. 
%               of the aircraft in the coordinate frame described above
%
%   strucgrid: A structure with the following fields
%       - numnodes: Number of FEM nodes
%       - coords: A 3-by-numnodes matrix where the i-th column contains 
%           [x;y;z] coordinate of the i-th node
%           The coordinate frame is such that:
%               A) Origin at the vehicle nose, 
%               B) +x pointing to rear or aircraft, and 
%               C) +y pointing starboard (to the right when looking
%                   forward)
%               D) +z pointing upward 
%       - numsensors: Number of sensors
%       - coords_sensors: A 3-by-numsensors matrix where the i-th column 
%           contains [x;y;z] coordinate of the i-th sensor
%           The coordinate frame is such that:
%               A) Origin at the vehicle nose, 
%               B) +x pointing to rear or aircraft, and 
%               C) +y pointing starboard (to the right when looking
%                   forward)
%               D) +z pointing upward 

%% Clear workspace
clear;
close all;
clc;

%% Setup
% Units conversions
units.in2m = 0.0254;          % Convert from inch to m, m/inch
units.m2in = 1/units.in2m;          % Convert from m to inch, inch/m
units.lbf2N = 4.44822;        % Convert from lbf to N, N/lbf
units.N2lbf = 1/units.lbf2N;        % Convert from N to lbf, lbf/N
units.lb2kg = 0.453592;       % Convert from lbs to kg, Kg/lb
units.kg2lb = 1/units.lb2kg;        % Convert from Kg to lbs, lb/Kg
units.slugs2kg = 14.5939;     % Convert from slugs to Kg, Kg/slugs
units.kg2slugs = 1/units.slugs2kg;  %  Convert from Kg to slugs, slugs/Kg

% Option to toggle figure displays
displayfigs = 1;                 

%% Options for model extraction
% Number of flexible modes 
%   (some modes could be discarded later, if found undesirable, eg. 
%   duplicates etc)
numflex_all = 8;                   

% Number of rigid modes 
%   (some modes could be discarded later, if found undesirable, eg. 
%   duplicates etc)
numrigid_all = 6;                       

% Total number of modes = Number of rigid + flexible modes
nummodes_all = numflex_all + numrigid_all;      

%% Aircraft mass properties
% This data can be obtained from 
%   1) The NASTRAN FEM model document on PAAW git
%   2) Mass properties used in the FD-Flex model. These are very close to 
%   the mass properties from the FEM model (obtained from Dave)
%   3) The mass properties test conducted at the UAV lab at UMN. The test 
%   result are summarized in:
%   'mAEWing1 - Mass Property Test Geri.xlsx' (included in the folder)

% Mass properties (from FD-Flex model)
Mass = 6.2513;                                  % Mass, kg     
Inertia = diag([2.6404, 0.4637, 2.7193]);       % Inertia, kg-m^2 

% CG Location
% This is the CG location in the final FEM coordinate frame with:
%   A) Origin at the vehicle nose, 
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward 
rcg = [0.6033;0;0];  % CG location, m

%% Load FEM information
% The description provided by Virginia Tech is as follows:
% Each node of the beam has 6 degrees of freedom (UX,UY,UZ,RX,RY,RZ). The
% number of effective nodes that compose the beam/rod elements is 
% numnodes=71.
% NASTRAN allows us to use massless points to interpolate the mode shapes
% for aerodynamic model in the aeroelastic analysis. These massless points 
% are connected to the beam/rod nodes through rigid massless bars. These 
% massless nodes are not taken into account in the global stiffness and 
% mass matrices assembling because the rigid massless bars have neither
% mass nor stiffness. Also, the mode shapes for those massless nodes can 
% be directly calculated based on that for the connected beam/rod nodes 
% because of the rigid connections. Hence, the number of the total degrees 
% of freedom is 71 \times 6 =426.

% Generate node labels:
%   Nodelabel: a numnodes-by-1 vector containing the label of the effective
%   nodes
run mAEWing1_nodes_label.m
clear node_label_mAEWing1
numnodes = numel(nodelabel);

% Load all FEM node coordinates: 
%   FEM.point_coordinates is an Nfem-by-4 matrix with columns containing
%   [Labels, x position, y position, z position]
%   where the coordinates are in *units of inches* and in a frame with:
%   A) Origin on vehicle centerline but ** -r0(1)(meter, defined later) 
%   aft of the nose**
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward
%   These coordinates include effective FEM beam nodes as well as extra 
%   nodes connected to the FEM beam nodes via massless elements. The FEM 
%   points do not contain the winglets or nodes for control surfaces.
%   Note that FEM is a structure with several fields but only 
%   points_coordinates field is used
load mAEWing1_Geri_WS4_FEMv5p2_modeinfo_v2.mat
FEMlabel_all = FEM.points_coordinates(:,1)';    % Labels
FEMcoord_all = FEM.points_coordinates(:,2:4)';  % [x,y,z] position, inch

% Location of Nose in temporary FEM coordinate frame in meters
%   A) Origin on vehicle centerline but ** -r0(1)(meter) 
%   aft of the nose**
%   B) +x pointing to rear or aircraft,
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward
r0 = [-0.3950; 0; 0];

%% Remove uneffective nodes
% Construct 3-by-numnodes coordinate matrix with columns equal to 
%   [x position; y position; z position]
%   wh%   The origin is shifted to the CG to ease the calculation of the
%   rigid body mode shapes.ere the coordinates are *units of m* and in a frame with:
%   A) Origin at the Nose, 
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward
%   This step removes nodes that are connected by massless elements.

% Index of effective nodes
idx_effective_nodes = ismember(FEMlabel_all,nodelabel);
% Coordinates of effective nodes in temporary FEM coordinate in meters
%   A) Origin on vehicle centerline but ** -r0(1)(meter) 
%   aft of the nose**
%   B) +x pointing to rear or aircraft,
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward
coords = units.in2m*FEMcoord_all(:,idx_effective_nodes);
% Shifting origin to the vehicle nose to obtain coordiantes in final
% coordinates:
%   A) Origin at the vehicle nose, 
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward 
coords = coords - r0;       

% Plot effective nodes
if displayfigs
    figure(1)
    plot3(coords(1,:),coords(2,:),coords(3,:),'*')
    xlabel('x'); ylabel('y'); zlabel('z'); grid on;
    view(90,90)
    title('Effective FEM nodes and CG location')
    
    % Add Location of CG
    hold on;
    plot3(rcg(1),rcg(2),rcg(3),'ro','MarkerSize',10,'LineWidth',4);
    legend('Effective FEM nodes','CG location','location','best')
    hold off;
    
    % Text # labels correspond to their order
    labels = cellstr(num2str(nodelabel));  
    text(coords(1,:), coords(2,:), coords(3,:), ...
        labels,'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right')
end

%% Obtain sensor node locations

% List of nodes corresponding to sensors
%   The order is [CBfor,CBaft,Lfor,Laft,Rfor,Raft,CG]
nodelabels_sensor = [9001:9006 18888];
% Number of sensors
numsensors = length(nodelabels_sensor);
% Find coordinates of sensor nodes
coords_sensor = zeros(3,numsensors);

for i = 1:numsensors
    % Index of sensor nodes in the full FEM nodes list () 
    idx_sensor = ismember(FEM.points_coordinates(:,1), ...
        nodelabels_sensor(i));
    
    % Sensor coordinates in final FEM frame where the coordinates 
    % are *units of m* and in a frame with:
    %   A) Origin at the Nose,
    %   B) +x pointing to rear or aircraft,
    %   C) +y pointing starboard (to the right when looking forward)
    %   D) +z pointing upward
    coords_sensor(:,i) = units.in2m*FEMcoord_all(:,idx_sensor) - r0;
    
    % Display sensors on the 
    if displayfigs
        % Add sensor locations in plot
        hold on;
        plot3(coords_sensor(1,i),coords_sensor(2,i),coords_sensor(3,i), ...
            'co','MarkerSize',8,'LineWidth',2);
        legend('Effective FEM nodes','CG location','Sensor Nodes', ...
            'location','best')
        hold off;
    end
end

%% Extract flexible model from FEM
% Load K and M matrices
%   These are obtained from VT and stored in FEM.K and FEM.M.
%   Both matrices are (6*numnodes)-by-(6*numnodes) because each node has 
%   six degrees of freedom. The matrices are expressed assuming nodal 
%   degrees of freedom are in units of (inches, lbf, rads). Units are
%   converted to SI for further processing.
load FEM_KM.mat

% Extract stiffness(K_temp) and mass(M_temp) matrices with nodal degrees of
% freedom [UX,UY,UZ,RX,RY,Rz] where
%   Ux = displacement towards tail (+X axis)
%   Uy = displacement towards right wing (starboard) (+Y axis)
%   UZ = displacement upward (+ Z axis)
%   RX = roll rotation about (rotaion aobut +X-axis)
%   RY = pitch rotation (rotation about +Y-axis)
%   Rz = yaw rotation (rotaition about +Z-axis)
K_temp=FEM.K;
M_temp=FEM.M;

% Convert mass/stiffness matrices (M_temp,K_temp) units from (in,lbf,rad) 
% to SI (m,N,rads)
% Displacement degrees of freedom are in (in,in,in,rad,rad,rad)
S_disp = repmat([units.m2in;units.m2in;units.m2in;1;1;1],[numnodes,1]);
% Force degrees of freedom are in (lbf,lbf,lbf,lbf-in,lbf-in,lbf-in)
S_force = repmat([units.N2lbf;units.N2lbf;units.N2lbf;
    units.N2lbf*units.m2in;units.N2lbf*units.m2in; ...
    units.N2lbf*units.m2in],[numnodes,1]);
% Scale K and M matrices to convert to SI units
K = lrscale(K_temp,1./S_force,S_disp);
M = lrscale(M_temp,1./S_force,S_disp);

%% Obtain model in modal coordinates
% Compute eigenvals, eigenvectors, and modal frequencies for all modes
[V_all,D_all]=eigs(K,M,nummodes_all,'sm');
wf=real(sqrt(diag(D_all)));     % All frquencies in rad/sec

% Flexible model
flexidx_all = (numrigid_all+1):nummodes_all;
Phif_all = V_all(:,flexidx_all);            % Flexible mode shapes 
Omegaf_all = wf(flexidx_all);               % Flexible modal freqs (rad/s)
Zetaf_all = 2/100*ones(numflex_all,1);      % Assumed damping

% Make surge, sway and yaw displacements zero
%   For Geri, the first six flexible modes contains negligible amount of
%   surge, sway and yaw deflection. These dofs are zero-ed out for
%   numerical accuracy.
%   Note: This step should be reconsidered for a new aircraft
Phif_all([1:6:end,2:6:end,6:6:end],:) = 0;

% Plot the flexible mode shapes
if displayfigs
    for i = 1:numflex_all  
    % Plot undeflected nodes
    figure
    plot3(coords(1,:),coords(2,:),coords(3,:),'b*')
    hold on
    
    % Plot deflected nodes
    plot3(coords(1,:)+Phif_all(1:6:end,i)',coords(2,:)+Phif_all(2:6:end,i)', ...
            coords(3,:)+Phif_all(3:6:end,i)','r*')
    xlabel('x'); ylabel('y'); zlabel('z'); grid on;
    view(110,60)
    title(['Flexible mode: ', num2str(i), ', Frequency: ', ...
        num2str(Omegaf_all(i)/2/pi),' Hz'])
    legend('Undeflected nodes','Deflected nodes')
    end
end

%% Discard unwanted modes
idx_modef_keep = [1,2,3,4,6];
% Note that 5th mode is not something identified during the FEM update and
% thus, it is ignored.

% Obtain modal properties of the flexible modes of interest.
numflex = length(idx_modef_keep);       % Numer of flexible modes kept
Phif = Phif_all(:,idx_modef_keep);      % Mode shapes
Omegaf = Omegaf_all(idx_modef_keep);    % Modal frequencies
Zetaf = Zetaf_all(idx_modef_keep);      % Modal damping

% Plot the final modes, if required
if displayfigs
    for i = 1:numflex     
    % Plot undeflected nodes
    figure
    plot3(coords(1,:),coords(2,:),coords(3,:),'b*')
    hold on
    
    % Plot deflected nodes
    plot3(coords(1,:)+Phif(1:6:end,i)',coords(2,:)+Phif(2:6:end,i)', ...
            coords(3,:)+Phif(3:6:end,i)','r*')
    xlabel('x'); ylabel('y'); zlabel('z'); grid on;
    view(110,60)
    title(['Flexible mode: ', num2str(i), ', Frequency: ', ...
        num2str(Omegaf(i)/2/pi),' Hz'])
    legend('Undeflected nodes','Deflected nodes')
    end
end

%% Modal mass, stiffness and damping matrices
Mf = diag(Phif'*M*Phif);
Kf = diag(Phif'*K*Phif);
Bf = 2*Zetaf.*Omegaf.*Mf;

%% Save data
strucmodes.Omegaf = Omegaf;      % Flexible frequencies (rad/sec)
strucmodes.Zetaf = Zetaf;        % Assumed modal damping
strucmodes.Phif = Phif;          % Flexible mode shapes, 6 dof
strucmodes.Mf = Mf;              % Modal mass matrix
strucmodes.Bf = Bf;              % Modal damping matrix
strucmodes.Kf = Kf;              % Modal stiffness matrix
strucmodes.M = M;                % Mass matrix
strucmodes.K = K;                % Stiffness matrix
strucmodes.Mass = Mass;          % Aircraft mass
strucmodes.Inertia = Inertia;    % Aircraft inertia matrix
strucmodes.rcg = rcg;            % CG location (nose as origin)

% Number of nodes and their coordinates
strucgrid.numnodes = numnodes;
strucgrid.coords = coords;

% Number of sensors and their coordinates
strucgrid.numsensors = numsensors;
strucgrid.coords_sensors = coords_sensor;

save Geri_strucmodel strucgrid strucmodes