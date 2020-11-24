function griddata = Geri_GridData
%GridData_Geri   Grids the aircraft
%
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
%           - PanelSpan: A 1-by-NumPanels vector where i-th element is the
%               span (distance between two edges parallel to the x-axis)
%               for the i-th panel
%           - nhat0: A 3-by-NumPanels vector where i-th column contains the
%               direction cosinse of the panel normal vector
%           - BodySpan: Span of the centrebody defined by the outer edge of
%               the body flaps
%           - NumFlaps: Number of control surface flaps
%           - PanelIndexCS: A (xCSDiv+xPFDiv)-by-(NyPF-1) matrix. This is 
%               the dimension of the panels of the planform. If 
%               PanelIndexCS(i,j) == k then the i,j-th panel belongs to the
%               control surface indexed k. The order of the control 
%               surfaces is [L1,L2,L3,L4,R1,R2,R3,R4]
%           - CSVertices: 3-by-4-by-NFlaps array where CSVertices(:,:,i) 
%               describes the [x,y,z] coordinates of the 4 vertices of the 
%               ith control surface in the order [L1,L2,L3,L4,R1,R2,R3,R4].
%               The vertices are ordered as [inner TE, inner hinge-line, 
%               outer hinge-line, outer TE]
%           - CSHingeDC: A 3-by-NumFlaps matrix where ith column contains
%               the direction cosine of the ith control surface hinge 
%               lines
%           - CSPanelIndex: A NumCS-by-1 cell vector contain the indices of
%               panels belonging to each control surfaces in order
%               [L1,L2,L3,L4,L5,L6,L7,L8];
%           - CSPanelIndex_all: A NumCSPanel vector contain the indices of
%               panels belonging to all control surfaces in order
%               [L1,L2,L3,L4,L5,L6,L7,L8];
%           - PFPanelIndex_all: A NumPFPanel vector contain the indices of
%               panels belonging to the planform but not to any control 
%               surfaces in order
%           - WLPanelIndex_all: A NumWLPanel vector contain the indices of
%               panels belonging to the winglets

%% Grid parameters
% Define number of divisions of the aircraft surface
yCSDiv = 5;    % Number of spanwise (+y) divisions of control surfaces
xCSDiv = 3;    % Number of streamwise (+x) divisions of control surfaces
xPFDiv = 9;    % Number of streamwise (+x) divisions of planform 
zWLDiv = 14;   % Number of hightwise (+z) divisions of winglet
% Note: xpfdiv defines the divisions in the +x direction of the planform
% but not including area associated with control surfaces.

%% Load aircraft geometry
[acgeom,LE,TE,lWLLE,lWLTE,rWLLE,rWLTE,CSVertices] = Geri_Outline;
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
NumFlap = size(yFlapEdge,1);
yPFGrid = zeros(NumFlap,yCSDiv+1);
for i=1:NumFlap
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
Nhat0 = [repmat([0;0;1],1,NumPanelsPF), ...
         repmat([0;1;0],1,2*NumPanelsWL)];

%% Panel lengths, areas and centroid points
Aerogrid = zeros(3,NumPanels);      % x,y,z coord of centre of each panel
PanelLength = zeros(1,NumPanels);   % Average xlength of each panel
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
    
    % Centroid of each panel
    Aerogrid(1,ii) = 0.25*(x1+x2+x3+x4); 
    Aerogrid(2,ii) = 0.25*(y1+y2+y3+y4); 
    Aerogrid(3,ii) = 0.25*(z1+z2+z3+z4);
end

%% Obtain Control surface panels and hinge data
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

% CS hingeline direction cosines
CSHingeDC = zeros(3,NumFlap);

% Cell of incides of control surface panels
CSPanelIndex = cell(NumFlap,1);

% Index of panels belonging to control surfaces
CSPanelIndex_all = [];
for i = 1:NumFlap
    HingePoint1 = CSVertices(:,2,i);    % Find ends of the hingeline
    HingePoint2 = CSVertices(:,3,i);
    CSHingeDC(:,i) = (HingePoint2-HingePoint1)/ ...
                        (norm(HingePoint2-HingePoint1));
    CSPanelIndex{i} = find(PanelIndexCS==i);
    CSPanelIndex_all = [CSPanelIndex_all find(PanelIndexCS==i)'];
end

% Index of planform panels not belonging to any control surface
PFPanelIndex_all = setdiff(1:NumPanelsPF,CSPanelIndex_all);
WLPanelIndex_all = NumPanelsPF+1:NumPanels;

%% Save final output
griddata.NodeData = NodeData;
griddata.PanelData = PanelData;
griddata.NumNodes = NumNodes;
griddata.NumPanels = NumPanels;
griddata.Aerogrid = Aerogrid;
griddata.PanelLength = PanelLength;
griddata.PanelSpan = PanelSpan;
griddata.Nhat0 = Nhat0;
griddata.BodySpan = BodySpan;
griddata.NumFlap = NumFlap;
griddata.PanelIndexCS = PanelIndexCS;
griddata.CSVertices = CSVertices;
griddata.CSHingeDC = CSHingeDC;
griddata.CSPanelIndex = CSPanelIndex;
griddata.CSPanelIndex_all = CSPanelIndex_all;
griddata.PFPanelIndex_all = PFPanelIndex_all;
griddata.WLPanelIndex_all = WLPanelIndex_all;

end