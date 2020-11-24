function griddata = GridData_Geri(rcg)
%GridData_Geri   Grids the aircraft
%   Function to generate grid for the Geri aircraft and return information 
%   required by the getAIC function. The aircraft geometry is given by the
%   function GeriDrawing. The grid is defined in a coordinate frame with-
%       A) Origin at the vehicle nose, 
%       B) +x pointing to rear or aircraft, and 
%       C) +y pointing starboard (to the right when looking forward)
%       D) +z pointing upward 
%
%   All the quantities are defined in SI
%
%   Input:
%       - rcg: A 3-by-1 vector cotaining [x;y;z] coordinates of the c.g of
%           the aircarft in the coordinate frame defined above
%   Outputs
%       - griddata: A structure containing the following fields
%           - NodeData: A 3-by-NumNodes matrix where i-th column contains
%               [x;y;z] coordinates of the i-th node
%           - PanelData: A 4-by-NumPanels matrix where i-th column contains
%               the NW, NE, SW, SE node indexes of the i-th panel
%           - NumNodes: Total number of nodes in the aerodynamic grid
%           - NumPanels: Total number of panels in the aerodynamic grid
%           - Aeorgrid: An 3-by-NumPanels matrix where i-th column contains
%               the [x;y;z] coordinates of the centroid of the i-th panel
%           - PanelLength: A 1-by-NumPanels vector where i-th element is 
%               the average length of the panel in x-direction
%           - PanelAera: A 1-by-NumPanel vector where i-th element is the 
%               area of the i-th panel
%           - PanelSpan: A 1-by-NumPanels vector where i-th element is the
%               span (distance between two edges parallel to the x-axis)
%               for the i-th panel
%           - S: A (6*NumPanels)-by-NumPanels matrix where where blocks of 
%               six rows give the [Ax; Ay; Az; Mx; My; Mz] data.
%               Ax is the panel area projected perpendicular to the x-axis
%               Ay and Az are defined similarly.
%               Mx is th moment of the area about the x-axis. My and Mz are
%               defined similarly. 
%           - nhatPF: A 1-by-NumPanels vector where i-th element is 1 if 
%               the panel belongs to the planform and 0 if it belongs to 
%               the winglets
%           - nhatWL: A 1-by-NumPanels vector where i-th element is 1 if 
%               the panel belongs to the winglets and 0 if it belongs to 
%               the planform
%           - BodySpan: Span of the centrebody defined by the outer edge of
%               the body flaps
%           - PhiCS_aero: A (6*NumPanels)-by-NumCS matrix where 
%               PhiCS_aero(6*i-5:6*i,k) contains 
%               [x,y,z,theta_x,theta_y,theta_z] degrees of freedom for the 
%               i-the panel for 1 radian deflection mode of k-th control 
%               surface.
%           - Phib_aero: A (6*NumPanels)-by-NumRigidModes(=6) matrix where 
%               the i-th column of Phib_aero contains 
%               [x,y,z,theta_x,theta_y,theta_z] degrees of freedom for the
%               the aerodynamic panels for the 
%               [surge,sway,heave,roll,pitch,yaw] rigid body modes

%% Setup
displayfigs = true;     % Option to toggle figure display

%% Grid parameters
% Note: xpfdiv defines the divisions in the +x direction of the planform
% but not including area associated with control surfaces.
yCSDiv = 5;    % Number of spanwise (+y) divisions of control surfaces
xCSDiv = 3;    % Number of streamwise (+x) divisions of control surfaces
xPFDiv = 9;    % Number of streamwise (+x) divisions of planform 
zWLDiv = 14;   % Number of hightwise (+z) divisions of winglet

%% Load aircraft geometry
[acgeom,LE,TE,lWLLE,lWLTE,rWLLE,rWLTE,BodySpan,CSVertices] = Drawing_Geri;
% acgeom contains dimensions of the full aircraft e.g rootchort, wingspan
% LE, TE contains x,y,z coordinates of points on the leading edge and
% trailing edge. The LE and TE can be obtained by linear interpolation of
% the points
% lWLLE contains points on the left-winglet-leading-edge
% lWLTE, rWLLE and rWLTE similarly define the edges of the winglets

%% Planform grid
% Grid is defined by:
%   1) Parallel lines along the flow (+x) direction as specified by the 
%      points in the vector yPlanform.
%   2) Lines (not necessarily parallel) as specified by points in the 
%      chordwise direction.  The surface flaps can be gridded at
%      a different density than the remainder of the planform.

% y-grid: Includes centerline, wingtips and a specified number of
% streamwise divisions of control surfaces.
yFlapEdge = acgeom.AllFlapYEdges;
NFlap = size(yFlapEdge,1);
yPFGrid = zeros(NFlap,yCSDiv+1);
for i=1:NFlap
    yPFGrid(i,:) = linspace(yFlapEdge(i,1),yFlapEdge(i,2),yCSDiv+1);
end
b = acgeom.WingSpan;
yPFGrid = unique([0; -b/2; b/2; yPFGrid(:)]);

% x-grid: This involves chordwise gridding that splits the planform into
% two regions. Specifically, define CS as a single line from the left to
% right wingtip that parallels the TE and contains all control surfaces.
% The chordwise gridding is defined with one gridding from leading edge 
% to the CS line and another gridding from CS to the trailing edge.
%
% xPFGrid(i,:) provides the x positions for the grid locations at
% y-position yPFGrid(i).
NyPF = numel(yPFGrid);
NxPF = xCSDiv+xPFDiv+1;
xPFGrid = zeros(NyPF,NxPF);

BodySpan = acgeom.BodySpan;
cbf = acgeom.BodyFlapChordSize;
cwf = acgeom.WingFlapChordSize;
for i=1:NyPF
    % Linearly interpolate to find location on LE and TE
    yi = yPFGrid(i);
    xLE = interp1( LE(2,:), LE(1,:), yi);
    xTE = interp1( TE(2,:), TE(1,:), yi);
    
    % Define line CS that contains control surfaces
    if abs(yi)<= BodySpan/2
        % Inside Body / Fuselage
        xCS = xTE - cbf;
    else
        % On wing
        xCS = xTE - cwf;
    end
    
    % Define grids from LE to CS and from CS to TE
    x1 = linspace(xLE,xCS,xPFDiv+1);
    x2 = linspace(xCS,xTE,xCSDiv+1);
    
    % Combine and store x grid
    xPFGrid(i,:) = [x1 x2(2:end)];
end

%% Winglet grid
% Grid is defined by:
%   1) Parallel lines along the flow (+x) direction as specified by the 
%      points in the vector zWLGrid, applicable to both winglets
%   2) Lines (not necessarily parallel) as specified by points in the 
%      chordwise direction. The number of divisions in x-direction on the
%      winglets is same as the total number of divisions in x-direction of
%      the wing planform (including control surfaces). Although it is
%      possible that the gridpoints on the winglets do not align with that
%      on the wing planform (including control surfaces)

% z-grid: Includes top,bottom edges of the winglet at specified number of
% divisions.
zWLBot = lWLLE(3,1);
zWLTop = lWLLE(3,end);
zWLGrid = linspace(zWLBot,zWLTop,zWLDiv+1);

% x-grid: Chordwise gridding for the winglets at each z-grid location
%
% xPFGrid(i,:) provides the x positions for the grid locations at
% z-position zW:Grid(i).
NzWL = numel(zWLGrid);
NxWL = xCSDiv+xPFDiv+1;
xWLGrid = zeros( NzWL, NxWL );

for i=1:NzWL
    % Linearly interpolate to find corresponding location on LE and TE
    zi = zWLGrid(i);
    xLE = interp1( lWLLE(3,:), lWLLE(1,:), zi);
    xTE = interp1( lWLTE(3,:), lWLTE(1,:), zi);
    
    % Define grids from LE to CS and from CS to TE
    xWLGrid(i,:) = linspace(xLE,xTE,xPFDiv+xCSDiv+1);
end

%% Generate NodeData
% NodeData is 3-by-NumNodes matrix where i-th column contains [x;y;z] 
% coordinates of the i-th node

% For planform (including control surfaces)
NumNodesPF = numel(xPFGrid);
NodeDataPF = zeros(3,NumNodesPF);  

for ii = 1:NyPF         % For grid line parallel to x-axis
    for jj = 1:NxPF     % For each node on the grid line
        nodei = (ii-1)*NxPF + jj;
        NodeDataPF(:,nodei) = [xPFGrid(ii,jj); yPFGrid(ii); 0];
    end
end

% For winglets
NumNodesWL = numel(xWLGrid);
NodeDataWL1 = zeros(3,NumNodesWL);      % Left winglet

for ii = 1:NzWL         % For each grid line parallel to x-axis
    for jj = 1:NxWL     % For each node on the grid line
        nodei = (ii-1)*NxWL + jj;
        NodeDataWL1(:,nodei) = [xWLGrid(ii,jj); -b/2; zWLGrid(ii)];
    end
end

% Right winglet
NodeDataWL2 = NodeDataWL1;
NodeDataWL2(2,:) = b/2;

% Combine planform and winglet nodedata
NodeData = [NodeDataPF, NodeDataWL1, NodeDataWL2];
NumNodes = size(NodeData,2);        % Number of nodes

%% Generate PanelData
% For planform (including control surfaces)

% Collect node-index of node at (-x,-y) corner (Northwest in top view with 
% nose point north) for each panel
NWnodesPF = 1:NxPF*(NyPF-1);        % Exclude points on TE 
NWnodesPF((1:NyPF-1)*NxPF) = [];    % Delete points on right wingtip

% PanelDataPF contains nodes index of 4 nodes 
% Note: Viewing from top, if LE is toward north, left is towards west. The
% panel include [NW, NE, SW, SE nodes].
% Note that this order is assumed in DLM / VLM codes further on
PanelDataPF = [NWnodesPF', NWnodesPF'+NxPF, NWnodesPF'+1, ...
                NWnodesPF'+NxPF+1]';
NumPanelsPF = size(NWnodesPF,2);            

% For Left Winglet
% Collect node-index of node at (-x,-z) corner (NorthEest in side view from
% with right nose pointing north) for each panel
NEnodesWL = (1:NxWL*(NzWL-1));
NEnodesWL((1:NzWL-1)*NxWL) = [];
NEnodesWL = NEnodesWL + NumNodesPF;
NumPanelsWL = size(NEnodesWL,2);

% PanelDataWL1 contains nodes index of 4 nodes of each panel on left
% winglet
PanelDataWL1 = [NEnodesWL', NEnodesWL'+NxWL, NEnodesWL'+1, ...
                NEnodesWL'+NxWL+1]';

% For right winglet            
NEnodesWL = NEnodesWL + NumNodesWL;
PanelDataWL2 = [NEnodesWL', NEnodesWL'+NxWL, NEnodesWL'+1, ...
                NEnodesWL'+NxWL+1]'; 

% Combine to create PanelData for full aircraft
PanelData = [PanelDataPF, PanelDataWL1, PanelDataWL2]; 
NumPanels = size(PanelData,2);      % Number of panels

%% Booleans to identify nodes on planform or winglets
nhatPF = [ones(1,NumPanelsPF), zeros(1,2*NumPanelsWL)];
nhatWL = -1*~nhatPF;

%% Panel lengths, areas and centroid points
Aerogrid = zeros(3,NumPanels);      % x,y,z coord of centre of each panel
PanelLength = zeros(1,NumPanels);   % Average xlength of each panel
PanelArea = zeros(1,NumPanels);     % Area of each panel
PanelSpan = zeros(1,NumPanels);     % Span of each panel

for ii = 1:NumPanels
    x1 = NodeData(1,PanelData(1,ii));   % X coord of each node 
    x2 = NodeData(1,PanelData(2,ii));
    x3 = NodeData(1,PanelData(3,ii));
    x4 = NodeData(1,PanelData(4,ii));    

    y1 = NodeData(2,PanelData(1,ii));   % y coord of each node
    y2 = NodeData(2,PanelData(2,ii));
    y3 = NodeData(2,PanelData(3,ii));
    y4 = NodeData(2,PanelData(4,ii));

    z1 = NodeData(3,PanelData(1,ii));   % z coord of each node
    z2 = NodeData(3,PanelData(2,ii));
    z3 = NodeData(3,PanelData(3,ii));
    z4 = NodeData(3,PanelData(4,ii));
    
    PanelLength(ii) = ((x3+x4)/2) - ((x1+x2)/2); % Avg xlengh of a panel
    PanelSpan(ii) = sqrt((y2 - y1)^2 + (z2 - z1)^2);
    
    %Area of each panel
    PanelArea(ii) = ((x3-x1)+(x4-x2))*(((y2-y1)/2) + ((z2-z1)/2));
    
    % Centroid of each panel
    Aerogrid(1,ii) = 0.25*(x1+x2+x3+x4); 
    Aerogrid(2,ii) = 0.25*(y1+y2+y3+y4); 
    Aerogrid(3,ii) = 0.25*(z1+z2+z3+z4);
end

%% Construct rigid-body mode shapes
% There are six rigid body mode shapes corresponding to translation and
% rotation along the (x,y,z) axes using the following axes:
%   A) Origin at the **vehicle CG**, 
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward

% Note that the final modes are expressed in the same coordinate system as
% before.
%   A) Origin at the vehicle nose, 
%   B) +x pointing to rear or aircraft, and 
%   C) +y pointing starboard (to the right when looking forward)
%   D) +z pointing upward 

% Aerodynamic grid with origin at the CG
Aerogrid_cg = Aerogrid-rcg;

% Translational mode shapes along the (x,y,z) directions.
TransMode = zeros(6*NumPanels,3);
TransMode(1:6:(6*NumPanels),1) = 1;     % Translation along +x direction
TransMode(2:6:(6*NumPanels),2) = 1;     % Translation along +y direction
TransMode(3:6:(6*NumPanels),3) = 1;     % Translation along -z direction

% Rotational mode shape along the x direction.
% The matrix for a rotation by TH (rads) around the x-axis is: 
% Rx = [1 0 0; 0 cos(TH) -sin(TH); 0 sin(TH) cos(TH)];
% This causes each nodal point to shift from [x;y;z] to Rx*[x;y;z].
% It also causes the nodal rotational degree of freedom about x to be TH.  
% The mode shape is constructed with an infinitesimal rotation, TH<<1.
% This gives the linearized rotation matrix:
%  Rx = eye(3) + [0 0 0; 0 0 -TH; 0 TH 0]
% The next translational is Rx*[x;y;z] - [x;y;z]. This is given by
%     TH*Tx*[x;y;z] where Tx = [0 0 0; 0 0 -1; 0 1 0]
% The (small) factor of TH is not included in the mode shape.
Tx = [0 0 0; 0 0 -1; 0 1 0];
Trans = (Tx*Aerogrid_cg)';                          % Translation
Rot = [-ones(NumPanels,1) zeros(NumPanels,2)];      % Rotation
xRotMode = [Trans Rot]';
xRotMode = xRotMode(:);

% Rotational mode shape along the y direction.
% The matrix for a rotation by TH (rads) around the y-axis is: 
%   Rx = [cos(TH) 0 sin(TH);0 1 0; -sin(TH) 0 cos(TH)];
% Following the same discussion above, yields the net translation:
%   TH*Ty*[x;y;z] where Ty = [0 0 1; 0 0 0; -1 0 0]
Ty = [0 0 1; 0 0 0; -1 0 0];
Trans = (Ty*Aerogrid_cg)';                   % Translation
Rot = [zeros(NumPanels,1) ones(NumPanels,1) zeros(NumPanels,1)];    % Rotation
yRotMode = [Trans Rot]';
yRotMode = yRotMode(:);

% Rotational mode shape along the z direction.
% The matrix for a rotation by TH (rads) around the z-axis is: 
%   Rz = [cos(TH) -sin(TH) 0; sin(TH) cos(TH) 0; 0 0 1];
% Following the same discussion above, yields the net translation:
%   TH*Tz*[x;y;z] where Tz = [0 1 0; -1 0 0; 0 0 0]
Tz = [0 1 0; -1 0 0; 0 0 0];
Trans = (Tz*Aerogrid_cg)';                   % Translation
Rot = [zeros(NumPanels,2) -ones(NumPanels,1)];    % Rotation
zRotMode = [Trans Rot]';
zRotMode = zRotMode(:);

% Stack translational and rotational mode shapes
Phib = [TransMode xRotMode yRotMode zRotMode];

% Plot rigid body modes
if displayfigs
    for i = 1:size(Phib,2)
        figure
        plot3(Aerogrid(1,:),Aerogrid(2,:),Aerogrid(3,:),'b*')
        hold on
        plot3(Aerogrid(1,:)+Phib(1:6:end,i)',...
                Aerogrid(2,:)+Phib(2:6:end,i)', ...
                Aerogrid(3,:)+Phib(3:6:end,i)','r*')
        xlabel('x'); ylabel('y'); zlabel('z'); grid on;
        legend('Undeflected nodes','Deflected nodes')
        axis equal
        view(110,60)
        title(['Rigid mode: ', num2str(i)])
    end
end

%% Obtain control surface modes
% PanelIndexCS (xCSDiv+xPFDiv)-by-(NyPF-1) matrix. This is dimension of the
% panels of the planform.
% if PanelIndexCS(i,j) == k then the i,j-th panel belongs to the control
% surface indexed k. The order of the control surfaces is
% [L1,L2,L3,L4,R1,R2,R3,R4]

PanelIndexCS = zeros(xCSDiv+xPFDiv,NyPF-1);
for i = 1:NyPF-1
    yi = 0.5*(yPFGrid(i) + yPFGrid(i+1));
    idx = find( yi>=yFlapEdge(:,1) & yi<=yFlapEdge(:,2) );
    if ~isempty(idx)
        PanelIndexCS(xPFDiv+1:end,i) = idx;
    end
end

% PhiCS(:,k) contain the control surface modes for 1 radian deflection of
% k-th control surface.
PhiCS = zeros(NumPanels*6,NFlap);
for CSidx = 1:NFlap % For each CS panel
    CSPanel = find(PanelIndexCS == CSidx); % Current panel number
    NumCSPanel = length(CSPanel);
    HingePoint1 = CSVertices(:,2,CSidx);    % Find ends of the hingeline
    HingePoint2 = CSVertices(:,3,CSidx);
    HingeAngle = atan2d(HingePoint2(2)-HingePoint1(2), ...  % Hinge angle
                        HingePoint2(1)-HingePoint1(1));
    % Perpendicular distance of the panels from hinge line
    CSPanelArmLen = vecnorm(cross(repmat(HingePoint2-HingePoint1, ...
                                                    1,NumCSPanel), ...
                    HingePoint2-Aerogrid(:,CSPanel)),2,1) / ...
                    norm(HingePoint2-HingePoint1);
    
    % Translation degrees of freedome along x,y,z direction
    PhiCS(CSPanel*6-5,CSidx) = zeros(NumCSPanel,1);
    PhiCS(CSPanel*6-4,CSidx) = zeros(NumCSPanel,1);
    PhiCS(CSPanel*6-3,CSidx) = -CSPanelArmLen';
    % Rotation degrees of freedom about x,y,z direction
    PhiCS(CSPanel*6-2,CSidx) = repmat(cosd(HingeAngle),1,NumCSPanel);
    PhiCS(CSPanel*6-1,CSidx) = repmat(sind(HingeAngle),1,NumCSPanel);
    PhiCS(CSPanel*6,CSidx) = zeros(NumCSPanel,1);
end

% Plot rigid body modes
if displayfigs
    for i = 1:size(PhiCS,2)
        figure
        plot3(Aerogrid(1,:),Aerogrid(2,:),Aerogrid(3,:),'b*')
        hold on
        plot3(Aerogrid(1,:)+PhiCS(1:6:end,i)',...
                Aerogrid(2,:)+PhiCS(2:6:end,i)', ...
                Aerogrid(3,:)+PhiCS(3:6:end,i)','r*')
        xlabel('x'); ylabel('y'); zlabel('z'); grid on;
        legend('Undeflected nodes','Deflected nodes')
        axis equal
        view(100,30)
        title(['Control surface mode: ', num2str(i)])
    end
end
%% Create area matrix (Skj)
% S is a 6*NumPanels-by-NumPanels.... (Check why this needs to be large and sparse
% after reviewing the DLM/VLM code)
% where blocks of six rows give the following data:
%    [Ax; Ay; Az; Mx; My; Mz]
% Ax is the panel area projected perpendicular to the x-axis. Ay and
% Az are defined similarly.  Mx is th moment of the area about the
% x-axis 
% S matrix contains areas and moment of areas about(quarter chord point)
% in rows, (3,9,15, ... ) a and (5,11,17,...) etc
%
% This seems to use local coordinates for S. If yes, then is there
% any reason to have area projections? In other words, define the
% axes such that the panel normal is along +/- z (so that only Az
% is non-zero for all panels including those on planform, winglets...)
%
% Double check the moment calculation for non-rectangular panels.
% The perpendicular distance from the centroid to the quarter chord
% line is not quite (?) PanelLength/4.
% S matrix contains areas and moment of areas about(quarter chord point)
% in rows, (3,9,15, ... ) a and (5,11,17,...) etc
S = zeros(6*NumPanels,NumPanels);
S(6*(1:NumPanels) - 3,:) = diag(PanelArea);
S(6*(1:NumPanels) - 1,:) = diag(PanelArea.*PanelLength*0.25);

%% Save final output
griddata.NodeData = NodeData;
griddata.PanelData = PanelData;
griddata.NumNodes = NumNodes;
griddata.NumPanels = NumPanels;
griddata.Aerogrid = Aerogrid;
griddata.PanelLength = PanelLength;
griddata.PanelArea = PanelArea;
griddata.PanelSpan = PanelSpan;
griddata.S = S;
griddata.nhatPF = nhatPF;
griddata.nhatWL = nhatWL;
griddata.BodySpan = BodySpan;
griddata.PhiCS_aero = PhiCS;
griddata.Phib_aero = Phib;
% save('Gridding/GridData_Geri.mat','griddata'); 