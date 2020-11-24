function griddata = GridData(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate grid for the BFF and return information required by
% the getAIC function. The aircraft geometry is given by function
% Geri_drawing() function is used to obtain halfwing edge data
%
% Input:
%   -NStr - no. of spanwise divisions on the wing (strips)
%   -N: no. of chordwise divisions of airframe exculding control surfaces
%   -C: no. of chordwise divisions of control surfaces
%   -Nz: no. of divisions along z axis for winglets
%
% Output:
%    -griddata:
%         - NodeData: NumNodes by 4 matrix where each row contains
%           [node_index xcoord ycoord zcoord] of each node
%         - PanelData: NumPanels-by-5 matrix where each row contains
%           [panel_index node_index1 node_index2 node_index3 node_index4]
%         - num_of_nodes
%         - num_of_panels
%         - aerogrid: 3-by-NumPanels matrix where each column contains 
%           [xcoord ycoord zcoord] of centriod of each panel
%         - PanelLength: 1-by-NumPanel vector where each element contains
%           the average x-length of each panel
%         - PanelAera: 1-by-NumPanel vector containing area of each panel
%         - S: NumPanel*6-by_NumPanel matrix containing aeras and moment
%           of areas about quarter chord point in rows, (3,9,15, ... ) 
%           and (5,11,17,...) etc
%         - phiCS: NumPanel*6-by-NumCS matrix containing 
%           [surge,sway,heave,roll,pitch,yaw] dofs of each panel for each
%           CS deflection
%         - PSpan: 1-by-NumPanel vector, span (or vertical) length of each 
%           panel
%         - n_hat_wingbody: NumPanel-by-1 vector containing 1 for panels
%           belonging to wingbody, others are 0
%         - n_hat_wl: NumPanel-by-1 vector containing 1 for panels beloning
%           to winglet, others are 0
%         - CSPanel: NumPanel-by-NumCS matrix where (i,j) element is 1 if 
%           it ith panel belongs to the jth control surface.
%          
% 
% Code developed by Aditya Kotikalpudi
% Date: 8/12/2014
% University of Minnesota
%
% Modified by: 
% Abhineet Gupta
% Graduate assistant
% University of Minnesota
% Date: 11/22/2017
% 
%
% Original Grid parameters were:
% NStr = 70;
% N = 9;
% Nz = 15;
% C = 3;

%% Obtain gridding parameters
switch nargin
    case 1  % If gridding data is not provided
        
        %%% CS divisions %%%
        % Spanwise division of each CS
        csdiv = 5;
        % Spanvise div of panels close to TE (including CS)
        NStr = [1,csdiv,csdiv,csdiv,1,csdiv,1,1,csdiv,1,csdiv,csdiv ...
            ,csdiv,1]; 
        % Chordwise divs of panels close to TE (including CS)
        C = round(3);  
        
        %%% Wing divisions %%%
        N = round(9);           % Chordwise divisions (except CS)
        
        %%% Winglet divs %%%
        Nz = round(15);             % Vertical divs
        WL_chordwise = round(N+C);   % Chordwise divs
        
        Components = varargin{1};
    case 5 % If gridding data is provided
        error('Pass single argument')
%         NStr = varargin{1};
%         N = varargin{2}; 
%         csdiv = varargin{3};
%         NStr_cs = [2,csdiv,csdiv,csdiv,1,csdiv,1,1,csdiv,1,csdiv,csdiv ...
%             ,csdiv,2]; 
%         C = varargin{4};
%         Nz = varargin{5};
%         WL_chordwise = varargin{6};
%         Components = varargin{7};        
    otherwise
        error('Check number of arguments to GridData.m')
end

saveloc = Components.path_data;   
filename = Components.Grid_file;

%% Obtain aircraft edge coordinates
% Half wing coordinates
[x,y,xwu,zwu,xwl,zwl,y_sec] = Geri_drawing;
% x: X-coords of LE and TE of body+wing
% y: Y-coords of LE or TE of body+wing
% xwu: x-coords of LE and TE of upper winglet
% zwu: z-coords of LE or TE of upper winglet
% xwl: x-coords of LE and TE of lower winglet
% zwl: z-coords of LE or TE of lower winglet

% Obtaining full wing coordinates from half wing
x_full_wingbody = [fliplr(x) x(:,2:size(x,2))]; 
y_full_wingbody = [-fliplr(y) y(2:length(y))];
z_full_wingbody = zeros(size(y_full_wingbody));
x_full_winglet = [fliplr(xwl) xwu(:,2:end)];
z_full_winglet = [fliplr(zwl) zwu(2:end)];
y_full_winglet_L = ones(size(z_full_winglet))*y_full_wingbody(1);
y_full_winglet_R = ones(size(z_full_winglet))*y_full_wingbody(end);

% Panels close to trailing edge (including CS) are divided into spanwise
% sections to better resolve the grid
% Obtain section inromation CS and panels close to CS
Ycoord_sec_full = [-fliplr(y_sec) y_sec(2:end)];

L_z0 = length(zwl); % Number of points on the lower half of the winglets

% Plot the leading edge and trailing edge of the wing and the winglets
figure
plot3(x_full_wingbody,y_full_wingbody,z_full_wingbody,'k');
title('Aerodynamic Grid')
axis([-1.75 1.75 -1.75 1.75 -1.75 1.75]);
hold on;
plot3(x_full_winglet,y_full_winglet_L,z_full_winglet,'k');
plot3(x_full_winglet,y_full_winglet_R,z_full_winglet,'k');
axis equal

%% define spanwise division (strips) based on Aspect Ratio-AR  

% Find the indices in edgecoordinates closest to gridline
idx = [];
for i = 1:(length(Ycoord_sec_full)-1)    % For each section
    ystart = Ycoord_sec_full(i);      % ycord of left edge
    yend = Ycoord_sec_full(i+1);      % ycord of right edge
    
    % ycoord of chordwise gridlines
    ygrid = linspace(ystart,yend,NStr(i)+1); 
    
    % Delete the rightmost gridline for all sections except for the last
    % one as they will be included as first edge of next section
    if i~=(length(Ycoord_sec_full) -1)
        ygrid(end) = [];
    end
    
    % Same as for wingbody
    for j = 1:length(ygrid)
        % Find indices of edgecoordinates withing a tolerence of gridline
        idxF = find(abs(y_full_wingbody - ygrid(j))<1e-5); 
        if length(idxF)>1
            % Pick the one just right of the median of all selected points
            idx_temp = ceil((idxF(1)+idxF(end))/2); 
        else
            idx_temp = idxF;
        end
        idx = [idx idx_temp];
    end

end

% No. of spanwise strips in TE region (including CS)
Nstrip = length(idx); 

%% Wingbody (without CS) chordwise division
% xgrid_wingbody_nocs contains xcoord of all the gridpoints in a x-by-y 
% grid
% N+1 since N panels in airframe cause N+1 nodes
xgrid_wingbody_noCS = zeros(N+1,Nstrip);   

% xcoords of beginning and end of strip lines
% (includes CS as of now, CS will be removed)
xgrid_begend_wingbody_noCS = x_full_wingbody(:,idx);

% ycoords of strip lines
ygrid_wingbody = y_full_wingbody(idx);

% length of strip lines (Includes CS as of now)
stripLength_wingbody = xgrid_begend_wingbody_noCS(2,:) ...
                       - xgrid_begend_wingbody_noCS(1,:);   

% xcoords of trailing point of strip lines
xcoord_TE_wingbody = xgrid_begend_wingbody_noCS(2,:); 

% xcoords of CS root points on the strip line
xctrlroot_wingbody = zeros(1,Nstrip); 


% CS root coordinate can be found by substracting CS chord from TE coords
% Body flaps
xctrlroot_wingbody((stripLength_wingbody>0.2992)) = ...
            xcoord_TE_wingbody((stripLength_wingbody>0.2992)) - 0.0946;  
% Wing flaps
xctrlroot_wingbody((stripLength_wingbody<0.2992)) = ...
            xcoord_TE_wingbody((stripLength_wingbody<0.2992)) - .0651;

% x lengths of the wingbody gridlines without CS
striplen_wingbody_chord = xctrlroot_wingbody - xgrid_begend_wingbody_noCS(1,:); 

% x-length of the pannels in wingbody without CS
panellen_wingbody_chord = striplen_wingbody_chord/N;   

% Divide each stripline of wingbody
for i = 1:Nstrip
    x1 = xgrid_begend_wingbody_noCS(1,i);   
    panLen = panellen_wingbody_chord(i);  
    
    % xcoord of each gridpoint
    xcol = [];
    for j = 1:N+1
        xcol = [xcol; x1 + (j-1)*panLen];
    end
    xgrid_wingbody_noCS(:,i) = xcol;
end
   
%% Chordwise division of panels close to TE including CS
xgrid_cs = zeros(C+1, Nstrip);   

% xcoords of each beginning and end of strip line
% (includes wingbody as of now)
xgrid_begend_CS = x_full_wingbody(:,idx);

% ycoordss of strip lines
ygrid = y_full_wingbody(idx);

% length of each strip line (Includes wingbody as of now)
stripleength_cs_chord = xgrid_begend_CS(2,:) - xgrid_begend_CS(1,:);   

% xcoord of trailing point of strip lines
xcoord_TE_CS = xgrid_begend_CS(2,:); 

%xcoord of CS root points on the strip line
xctrlroot_cs = zeros(1,Nstrip);

% CS root coordinate can be found by substracting CS chord from TE coords
% Body flaps
xctrlroot_cs((stripleength_cs_chord>0.2992)) = ...
    xcoord_TE_CS((stripleength_cs_chord>0.2992)) - 0.0946;
% Wing flaps
xctrlroot_cs((stripleength_cs_chord<0.2992)) = ...
    xcoord_TE_CS((stripleength_cs_chord<0.2992)) - 0.0651;

% x length of CS gridlines 
striplen_CS_chord = xgrid_begend_CS(2,:) - xctrlroot_cs;

% x length of CS panels
paanellen_CS_chord = striplen_CS_chord/C;

% Divide each stripline of CS panels
for i = 1:Nstrip
    x1 = xctrlroot_cs(i);
    CpanLen = paanellen_CS_chord(i);
    xcol = [];
    for j = 1:C+1
        xcol = [xcol; x1 + (j-1)*CpanLen];
    end
    xgrid_cs(:,i) = xcol;
end

   
%% plot the panels

% Plot wing strips
for ii = 1:Nstrip
x1 = xgrid_wingbody_noCS(:,ii); y1_l = ygrid_wingbody(ii)*ones(N+1,1);
plot(x1,y1_l,'k'); hold on;
end

% Plot wing chordwise divisions
for ll = 1:Nstrip-1
    for nn = 1:N+1
        x1 = [xgrid_wingbody_noCS(nn,ll);xgrid_wingbody_noCS(nn,ll+1)];
        y1_l = [ygrid_wingbody(ll);ygrid_wingbody(ll+1)];
        plot(x1,y1_l,'k'); hold on;
    end
end

% Plot CS strips
for ii = 1:Nstrip
x1 = xgrid_cs(:,ii); y1_l = ygrid(ii)*ones(C+1,1);
plot(x1,y1_l,'k'); hold on;
end

% Plot winglet chordwise divisions
for ll = 1:Nstrip-1
    for nn = 1:C+1
        x1 = [xgrid_cs(nn,ll);xgrid_cs(nn,ll+1)];
        y1_l = [ygrid(ll);ygrid(ll+1)];
        plot(x1,y1_l,'k'); hold on;
    end
end

%% Griding winglets
if~mod(Nz,2)
    Nz = Nz+1;  % Nz is always odd so that there is a line on the z=0
end

% zcoors of gridlines on upper and lower winglets
zcoord_grid_l = linspace(z_full_winglet(1), ...
                z_full_winglet(L_z0),0.5*(Nz-1)+1);
zcoord_grid_u = linspace(z_full_winglet(L_z0),...
    z_full_winglet(end),0.5*(Nz-1)+1);

zcoord_grid = [zcoord_grid_l zcoord_grid_u(2:end)]; % zcord of dividing approx panel lines

% Same as before
idz = []; 
for i = 1:length(zcoord_grid)
    idzF = find(abs(z_full_winglet - zcoord_grid(i))<1e-5);
    if length(idzF)>1
        idz_temp = ceil((idzF(1)+idzF(end))/2);
    else
        idz_temp = idzF;
    end
    idz = [idz idz_temp];
end

% xcord of LE and TE points of z strip lines on WL
xcoord_wl_LE_TE = x_full_winglet(:,idz);  

% xlength of panels on the z strip lines on WL
xwl_panellen = (xcoord_wl_LE_TE(2,:) - xcoord_wl_LE_TE(1,:))/WL_chordwise; 

% xcoord of gridpoint on the winglet
xcoord_wl_grid = (repmat([0:WL_chordwise]',1,length(xwl_panellen)).* ...
            repmat(xwl_panellen,WL_chordwise+1,1)) ... 
            + repmat(xcoord_wl_LE_TE(1,:),WL_chordwise+1,1);


% Plot winglet grid
for ll = 1:Nz-1
    for nn = 1:WL_chordwise+1
        x1 = [xcoord_wl_grid(nn,ll);xcoord_wl_grid(nn,ll+1)];
        y1_l = y_full_wingbody(1)*ones(2,1);
        y1_r = y_full_wingbody(end)*ones(2,1);
        z1 = [zcoord_grid(ll);zcoord_grid(ll+1)];
        plot3(x1,y1_l,z1,'r');
        plot3(x1,y1_r,z1,'r');        
    end
    x2 = xcoord_wl_LE_TE(:,ll);
    y2_l = y_full_wingbody(1)*ones(2,1);
    y2_r = y_full_wingbody(end)*ones(2,1);
    z2 = zcoord_grid(ll)*ones(2,1);
    plot3(x2,y2_l,z2,'r');
    plot3(x2,y2_r,z2,'r');
    hold on;
end

% Plot the topmost chord or wl
x2 = xcoord_wl_LE_TE(:,Nz);
y2_l = y_full_wingbody(1)*ones(2,1);
y2_r = y_full_wingbody(end)*ones(2,1);
z2 = zcoord_grid(Nz)*ones(2,1);
plot3(x2,y2_l,z2,'r');
plot3(x2,y2_r,z2,'r');
hold on;


%% Generate PanelData and NodeData
%%% Wingbody %%%
% Number of nodes on wingbody
num_of_nodes = size(xgrid_wingbody_noCS,1)*size(xgrid_wingbody_noCS,2);  

nodeIndex = 1:num_of_nodes; % Node index

% Nodedata contains [nodeindex,xcord,ycord, zcoord]
NodeData = zeros(num_of_nodes,4);  

% Calculating NodeData for wingbody 
% Note that node number goes [1:(N+1),(N+2):2(N+1)...] etc.
% In other words, they go from leftmost chord (LE to TE) and move right (LE
% to TE)
for ii = 1:Nstrip      % looping over every strip line
    for jj = 1:N+1      % looping over each node on a strip
        nodeNo = floor(nodeIndex(((ii-1)*(N+1)) + jj));
        NodeData(nodeNo,:) = [nodeNo xgrid_wingbody_noCS(jj,ii) ygrid_wingbody(ii) 0];
    end
end

%%% CS region %%%
% Number of nodes on cs
num_of_nodes_cs = size(xgrid_cs,1)*size(xgrid_cs,2);  
nodeIndex_cs = 1:num_of_nodes_cs; % Node index

% Nodedata contains [nodeindex,xcord,ycord, zcoord]
NodeData_cs = zeros(num_of_nodes_cs,4);  

% Calculating NodeData for wings
% Note that node number goes [1:(C+1),(C+2):2(C+1)...] etc.
% In other words, they go from leftmost chord (LE to TE) and move right (LE
% to TE)
for ii = 1:Nstrip       % looping over every strip line
    for jj = 1:C+1      % looping over each node on a strip
        nodeNo = floor(nodeIndex_cs(((ii-1)*(C+1)) + jj));
        NodeData_cs(nodeNo,:) = [nodeNo xgrid_cs(jj,ii) ...
            ygrid(ii) 0];
    end
end
NodeData_cs(:,1) = NodeData_cs(:,1) + num_of_nodes;

% NodeData for winglets ()
% Note that node number goes [1:(WL_chordwise+1), ...
% (WL_chordwise+2):2(WL_chordwise+1)...] etc.
% In other words, they go from lowermost chord (LE to TE) and move up (LE
% to TE)
num_of_nodes_wl1 = size(xcoord_wl_grid,1)*size(xcoord_wl_grid,2); 
nodeIndex_wl1 = 1:num_of_nodes_wl1; % node index
NodeData_wl1 = zeros(num_of_nodes_wl1,4);   % left winglet
NodeData_wl2 = zeros(num_of_nodes_wl1,4);   % right winglet

for ii = 1:Nz       % looping over every strip line
    for jj = 1:WL_chordwise+1   % looping over each node on a strip
        nodeNo = floor(nodeIndex_wl1(((ii-1)*(WL_chordwise+1)) + jj));
        NodeData_wl1(nodeNo,:) = [nodeNo xcoord_wl_grid(jj,ii) ...
                                  y_full_wingbody(1) zcoord_grid(ii)];
        NodeData_wl2(nodeNo,:) = [nodeNo xcoord_wl_grid(jj,ii) ...
                                  y_full_wingbody(end) zcoord_grid(ii)];

    end
end

NodeData_wl1(:,1) = NodeData_wl1(:,1) + num_of_nodes + num_of_nodes_cs;
NodeData_wl2(:,1) = NodeData_wl2(:,1) + num_of_nodes ...
                + num_of_nodes_cs + num_of_nodes_wl1;

Nnodes = NodeData_wl2(end,1);
NodeData = [NodeData;NodeData_cs;NodeData_wl1;NodeData_wl2];


%% PanelData
%%%  For Wings %%%
num_of_panels = (Nstrip-1)*N;  
panelIndex = 1:num_of_panels;

% Find panel data from the node data
SWnodes = nodeIndex;        % Same as number of nodes in the wings
SWnodes(end - N: end) = []; % Delete one of nodes along the each chord
Isw = [1:Nstrip-1]*(N+1); 
SWnodes(Isw) = []; % Delete one line of nodes along the span

%Now you're left with SWnodes containing 1 node coresponding to each
% panel (left leading point(double check if any doubt))

% Panel [number, nodes index of 4 nodes]
% Note: Viewing from top, if LE is toward north, left is towards west. The
% panel include [panel number, NW, NE, SW, SE nodes].
PanelData = [panelIndex' SWnodes' SWnodes'+N+1 SWnodes'+1 SWnodes'+N+2];

%%% For CS region %%%
% Same thing for winglets
num_of_panels_cs = (Nstrip-1)*C; 
panelIndex_cs = [1:num_of_panels_cs] + num_of_panels;

%%% For winglets %%%
num_of_panels_wl1 = (Nz-1)*(WL_chordwise);   
panelIndex_wl1 = [1:num_of_panels_wl1] + num_of_panels + num_of_panels_cs;
panelIndex_wl2 = panelIndex_wl1 + num_of_panels_wl1;

SWnodes_cs = NodeData_cs(:,1);
SWnodes_wl1 = NodeData_wl1(:,1);
SWnodes_wl2 = NodeData_wl2(:,1);

SWnodes_cs(end - C: end) = [];
Isw = [1:Nstrip-1]*(C+1);
SWnodes_cs(Isw) = [];

SWnodes_wl1(end - (WL_chordwise): end) = [];
Isw = [1:Nz-1]*(WL_chordwise+1);
SWnodes_wl1(Isw) = [];

SWnodes_wl2(end - (WL_chordwise): end) = [];
Isw = [1:Nz-1]*(WL_chordwise+1);
SWnodes_wl2(Isw) = [];


PanelData_cs = [panelIndex_cs' SWnodes_cs SWnodes_cs+C+1 ...
                SWnodes_cs+1 SWnodes_cs+C+2];
PanelData_wl1 = [panelIndex_wl1' SWnodes_wl1 SWnodes_wl1+WL_chordwise+1 ...
                SWnodes_wl1+1 SWnodes_wl1+WL_chordwise+2];
PanelData_wl2 = [panelIndex_wl2' SWnodes_wl2 SWnodes_wl2+WL_chordwise+1 ...
                SWnodes_wl2+1 SWnodes_wl2+WL_chordwise+2];
PanelData = [PanelData; PanelData_cs; PanelData_wl1; PanelData_wl2];  

%% Get normal vectors for the panels

% For wings, 1 for pannels containing 
n_hat_wingbodyCS = [ones(num_of_panels+num_of_panels_cs,1);
            zeros(2*num_of_panels_wl1,1)];
n_hat_wl = [zeros(num_of_panels+num_of_panels_cs,1);
            -ones(2*num_of_panels_wl1,1)];

%% Panel lengths, areas and center points
% Total number of panels
Npanel = num_of_panels + num_of_panels_cs +(2*num_of_panels_wl1);

aerogrid = zeros(3,Npanel); % x,y,z coord of centre of each panel
Plength = zeros(1,Npanel); % Average xlength of each panel
PAreas = zeros(1,Npanel); % Area of each panel

for ii = 1:Npanel
    x1 = NodeData(PanelData(ii,2),2);       % X coord of each node 
    x2 = NodeData(PanelData(ii,3),2);
    x3 = NodeData(PanelData(ii,4),2);
    x4 = NodeData(PanelData(ii,5),2);    

    y1_l = NodeData(PanelData(ii,2),3);     % y coord of each node
    y2_l = NodeData(PanelData(ii,3),3);
    y3 = NodeData(PanelData(ii,4),3);
    y4 = NodeData(PanelData(ii,5),3);

    z1 = NodeData(PanelData(ii,2),4);       % z coord of each node
    z2 = NodeData(PanelData(ii,3),4);
    z3 = NodeData(PanelData(ii,4),4);
    z4 = NodeData(PanelData(ii,5),4);
    
    Plength(ii) = ((x3+x4)/2) - ((x1+x2)/2); % Avg xlengh of a panel
    
    % As panels are either horizontal or verticle, therefore span of a
    % panel can be calculated as follows.
    PSpan(ii) = ((y2_l - y1_l) + (z2 - z1));
    
    %Area of each panel
    PAreas(ii) = ((x3-x1)+(x4-x2))*(((y2_l-y1_l)/2) + ((z2-z1)/2));
    
    
    % aerogrid contians the coordinate of centroid of each panel
    xr1=x1+(x3-x1)/2;   % 1/2chord x-location at the root of the panel
    xr2=x2+(x4-x2)/2;   % 1/2chord x-location at the tip of the panel
    
    % 1/2chord x-location at the half-span location of the panel
    aerogrid(1,ii) = 0.5*(xr1+xr2); 
    aerogrid(2,ii) = (y1_l+y2_l)/2;
    aerogrid(3,ii) = (z1+z2)/2;
end

% plot3(aerogrid(1,:),aerogrid(2,:),aerogrid(3,:),'bx');

%% create area matrix (Skj)
% S matrix contains aeras and moment of areas about(quarter chord point)
% in rows, (3,9,15, ... ) a and (5,11,17,...) etc

S = zeros(6*Npanel,Npanel);
rows1 = 6*[1:Npanel] - 3;
rows2 = 6*[1:Npanel] - 1;
S(rows1,:) = diag(PAreas(1,:));
S(rows2,:) = diag(PAreas.*Plength*0.25);

S = sparse(S);

%% create control surface definitions and mode shapes
% For 8 control surfaces and 6dof for each panel
phiCS = zeros(Npanel*6,8);  

% CS panel index
csIndex = ((num_of_panels + 1):(num_of_panels+num_of_panels_cs))';

% Binary variable, contains 1 in i,j position if ith panel is part of
% jth control surface
CSind = zeros(Npanel,8);

for ii = 1:length(csIndex) % For each CS panel
    Pno = csIndex(ii); % Current panel number
    ymid = aerogrid(2,Pno); % ycord of current CS panel midpoint
    
    % CSpanel markes the row number of the panel in CS
    CSpanel = mod(Pno-num_of_panels,C);
    
    % Ctarting with 1 at CS root
    if CSpanel == 0
       CSpanel = C;
    end
    
    % Calculates the mode shapes associated with CS deflections for
    % body flaps
    if (abs(ymid)<0.3405 && abs(ymid)>=0.0497)
        if ymid<0 % check if on the left half or right half
            % Calculates mode shape for 1 rad deflection, 
            phiCS((Pno*6 - 5):Pno*6,1) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3914; 0.9202;0];
            % Assign CSind as control surface 1
            CSind(Pno,1) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,5) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3914; 0.9202;0];
            CSind(Pno,5) = 1;
        end
    end
    
    % Same for other CS
    if (abs(ymid)<0.7437 && abs(ymid)>=0.3934)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,2) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3745; 0.9272;0];
            CSind(Pno,2) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,6) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3745; 0.9272;0];
            CSind(Pno,6) = 1;
        end
    end
            

    if (abs(ymid)<1.0941 && abs(ymid)>=0.7437)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,3) = ...
            [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3745; 0.9272;0];
            CSind(Pno,3) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,7) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3745; 0.9272;0];
            CSind(Pno,7) = 1;
        end
    end
    
    
    if (abs(ymid)<1.4444 && abs(ymid)>=1.0941)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,4) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3745; 0.9272;0];
            CSind(Pno,4) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,8) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3745; 0.9272;0];
            CSind(Pno,8) = 1;
        end
    end
    
end

phiCS = sparse(phiCS);

%% Save data
griddata.NodeData = NodeData;
griddata.PanelData = PanelData;
griddata.num_of_nodes = Nnodes;
griddata.num_of_panels = Npanel;
griddata.aerogrid = aerogrid;
griddata.PanelLength = Plength;
griddata.PanelArea = PAreas;
griddata.S = S;
griddata.phiCS = phiCS;
griddata.PSpan = PSpan;
griddata.n_hat_wingbody = n_hat_wingbodyCS;
griddata.n_hat_wl = n_hat_wl;
griddata.CSPanel = CSind;

save([saveloc,filename],'griddata'); 
end
