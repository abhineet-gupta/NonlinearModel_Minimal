function griddata = GridData_original(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate grid for the BFF and return information required by
% the getAIC function. The aircraft geometry is given by function
% BFFdrawing() .Note that BFFdrawing gives only half the wing.
% Input Variables: 
% NStr - no. of spanwise divisions on the wing (strips)
% N: no. of chordwise divisions of airframe exculding control surfaces
% C: no. of chordwise divisions of control surfaces
% Nz: no. of divisions along z axis for winglets

% Code developed by Aditya Kotikalpudi
% Date: 8/12/2014
% University of Minnesota
%
% Modified by: 
% Abhineet Gupta
% Graduate assistant
% University of Minnesota
% Date: 11/22/2017

switch nargin
    case 1
        scale=1;
        NStr = floor(scale*70);
        N = floor(scale*9); 
        C = floor(scale*3);
        Nz = floor(scale*15);
        Components = varargin{1};
    case 5
        NStr = varargin{1};
        N = varargin{2}; 
        C = varargin{3};
        Nz = varargin{4};
        Components = varargin{5};        
    otherwise
        error('Check number of arguments to GridData.m')
end

saveloc = Components.path_data;
filename = Components.Grid_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,xwu,zwu,xwl,zwl] = BFFdrawing_original;
% half wing coordinates, x has 2 rows; top row gives leading edge 
% coordinates, 2nd row gives trailing edge, for a given y
% x,y contains xcoord,ycoord of both LE and TE of the half wing
% xwu contains xcoord of winglet upper (both LE and TE) of one winglet
% Rest is easy to understand

% Obtaining full wing coordinates from half wing
x_full_wing = [fliplr(x) x(:,2:size(x,2))]; 
y_full_wing = [-fliplr(y) y(2:length(y))];
z_full_wing = zeros(size(y_full_wing));
x_full_winglet = [fliplr(xwl) xwu(:,2:end)];
z_full_winglet = [-fliplr(zwl) zwu(2:end)];
y_full_winglet_L = ones(size(z_full_winglet))*y_full_wing(1);
y_full_winglet_R = ones(size(z_full_winglet))*y_full_wing(end);

L_z0 = length(zwl); % Number of points on the lower half of the winglets

% Plot the leading edge and trailing edge of the wing and the winglets
figure(1)
 plot3(x_full_wing,y_full_wing,z_full_wing,'k');
 axis([-1.75 1.75 -1.75 1.75 -1.75 1.75]);
 hold on;
 plot3(x_full_winglet,y_full_winglet_L,z_full_winglet,'k');
 plot3(x_full_winglet,y_full_winglet_R,z_full_winglet,'k');

%% define spanwise division (strips) based on Aspect Ratio-AR

ystart = y_full_wing(1);    % ycord of left wing tip
yend = y_full_wing(end);    % ycord of right wing top
ypanel = linspace(ystart,yend,NStr); % ycoord of dividing lines
idx = [];
for i = 1:length(ypanel)
    % To find the index of gridline close to dividing line
    idxF = find(abs(y_full_wing - ypanel(i))<1e-5); 
    if length(idxF)>1
        % Pick the one just right of the center of all the selected points
        idx_temp = ceil((idxF(1)+idxF(end))/2); 
    else
        idx_temp = idxF;
    end
    idx = [idx idx_temp];   % index of the gridlines which will become strip
end

Nstrip = length(idx); % no. of strips

%% define chordwise divisions

% xplnel contains xcoord of all the gridpoints in a x-by-y grid
xpanel = zeros(N+1+C, Nstrip);   % N+1+C since N panels in airframe cause N+1 nodes, C extra nodes for control surfaces
xstrip = x_full_wing(:,idx);        % x coordinates of each beginning and end of strip line
ypanel = y_full_wing(idx);          % y coordinates of strip lines
StripLength = xstrip(2,:) - xstrip(1,:);   % length of each strip line

% xcoord of trailing point of strip lines
xtrail = xstrip(2,:); 

xctrl = zeros(1,Nstrip); % variable for xcoord of CS root points on the strip line
xctrl((StripLength>0.31)) = xtrail((StripLength>0.31)) - 0.102;   % Wing flaps top x coordinate
xctrl((StripLength<0.31)) = xtrail((StripLength<0.31)) - 0.0744;   % Body flaps top x coordinate

FrameLength = xctrl - xstrip(1,:);  % x-length of the airframe w/o CS
ctrlLength = xstrip(2,:) - xctrl;   % x-length of the CS
framePanelLength = FrameLength/N;   % x-length of the pannels in airframe w/o CS
ctrlPanelLength = ctrlLength/C;     % x-length of pannels in the CS

for i = 1:Nstrip
    x1 = xstrip(1,i);   % xcoord of leading points of each stripline
    
    FpanLen = framePanelLength(i);  % frame-planel length for stripline
    CpanLen = ctrlPanelLength(i);   % CS panel length for stripline
    xcol = [];          % Temporary variable used to create 'xpanel'
    for j = 1:N+1
        xcol = [xcol; x1 + (j-1)*FpanLen];
    end
    
    x2 = xcol(N+1,:);
    
    for j = 1:C
        xcol = [xcol; x2 + j*CpanLen];
    end
    
    % xpanel is a Num_chord_division - by - Nstrip matrix
    % It contains xcoord of both frame and CS combined
    xpanel(:,i) = xcol;
end
   

   
%% plot the panels
%plot(xpanel,ypanel,'ro')
% plot strips

for ii = 1:Nstrip
x1 = xpanel(:,ii); y1_l = ypanel(ii)*ones(N+C+1,1);
plot(x1,y1_l,'k'); hold on;
end

% plot chordwise divisions
for ll = 1:Nstrip-1
    for nn = 1:N+C+1
        x1 = [xpanel(nn,ll);xpanel(nn,ll+1)];
        y1_l = [ypanel(ll);ypanel(ll+1)];
        plot(x1,y1_l,'k'); hold on;
    end
end

%% griding winglets
if~mod(Nz,2)
    Nz = Nz+1;  % Nz is always odd so that there is a line on the z=0
end

zpanell = linspace(z_full_winglet(1),z_full_winglet(L_z0),0.5*(Nz-1)+1);
zpanelu = linspace(z_full_winglet(L_z0),z_full_winglet(end),0.5*(Nz-1)+1);
zpanel = [zpanell zpanelu(2:end)]; % zcord of dividing approx panel lines

idz = []; % variable for index of z-strip lines on winglet
for i = 1:length(zpanel)
    idzF = find(abs(z_full_winglet - zpanel(i))<1e-5);
    if length(idzF)>1
        idz_temp = ceil((idzF(1)+idzF(end))/2);
    else
        idz_temp = idzF;
    end
    idz = [idz idz_temp];
end

% xcord of LE and TE points of z strip lines on WL
x_winglet = x_full_winglet(:,idz);  

% xlength of panels on the z strip lines on WL
xw_panel_delta = (x_winglet(2,:) - x_winglet(1,:))/(N+C); 

% xw_panel is xcoord of gridpoint on the WL
xw_panel = (repmat([0:N+C]',1,length(xw_panel_delta)).* ...
            repmat(xw_panel_delta,N+C+1,1)) ... 
            + repmat(x_winglet(1,:),N+C+1,1);


% plot left winglet grid
for ll = 1:Nz-1
    for nn = 1:N+C+1
        x1 = [xw_panel(nn,ll);xw_panel(nn,ll+1)];
        y1_l = ystart*ones(2,1);
        y1_r = yend*ones(2,1);
        z1 = [zpanel(ll);zpanel(ll+1)];
       % plot3(x1(1),y1_l(1),z1(1),'bo');
       % plot3(x1(2),y1_l(2),z1(2),'bo');
        plot3(x1,y1_l,z1,'r');
       % plot3(x1(1),y1_r(1),z1(1),'bo');
       % plot3(x1(2),y1_r(2),z1(2),'bo');
        plot3(x1,y1_r,z1,'r');        
    end
    x2 = x_winglet(:,ll);
    y2_l = ystart*ones(2,1);
    y2_r = yend*ones(2,1);
    z2 = zpanel(ll)*ones(2,1);
    plot3(x2,y2_l,z2,'r');
    plot3(x2,y2_r,z2,'r');
    hold on;
end
x2 = x_winglet(:,Nz);
y2_l = ystart*ones(2,1);
y2_r = yend*ones(2,1);
z2 = zpanel(Nz)*ones(2,1);
plot3(x2,y2_l,z2,'r');
plot3(x2,y2_r,z2,'r');
hold on;


%% Generate PanelData and NodeData

% NodeData wings
num_of_nodes = size(xpanel,1)*size(xpanel,2);  % number of nodes on wings
nodeIndex = 1:num_of_nodes; % node index

% Nodedata contains [nodeindex,xcord,ycord, zcoord]
NodeData = zeros(num_of_nodes,4);  

% Calculating NodeData for wings
for ii = 1:Nstrip       % looping over every strip line
    for jj = 1:N+C+1      % looping over each node on a strip
        nodeNo = floor(nodeIndex(((ii-1)*(N+C+1)) + jj));
        NodeData(nodeNo,:) = [nodeNo xpanel(jj,ii) ypanel(ii) 0];
    end
end

% NodeData for winglets
num_of_nodes2 = size(xw_panel,1)*size(xw_panel,2);  % number of nodes
nodeIndex2 = 1:num_of_nodes2; % node index
NodeData2 = zeros(num_of_nodes2,4);   % left winglet
NodeData3 = zeros(num_of_nodes2,4);   % right winglet

for ii = 1:Nz       % looping over every strip line
    for jj = 1:N+C+1      % looping over each node on a strip
        nodeNo = floor(nodeIndex2(((ii-1)*(N+C+1)) + jj));
        NodeData2(nodeNo,:) = [nodeNo xw_panel(jj,ii) ystart zpanel(ii)];
        NodeData3(nodeNo,:) = [nodeNo xw_panel(jj,ii) yend zpanel(ii)];

    end
end

NodeData2(:,1) = NodeData2(:,1) + num_of_nodes;
NodeData3(:,1) = NodeData3(:,1) + num_of_nodes + num_of_nodes2;
Nnodes = NodeData3(end,1);
NodeData = [NodeData;NodeData2;NodeData3];


%% PanelData

%  For Wings

% Includes control surface panels, excludes winglets
num_of_panels = (Nstrip-1)*(N+C);  

panelIndex = 1:num_of_panels;

% Finding panel data from the node data
SWnodes = nodeIndex;        % Same as number of nodes in the wings
SWnodes(end - (N+C): end) = []; % Delete one of nodes along the each chord
Isw = [1:Nstrip-1]*(N+C+1); 
SWnodes(Isw) = []; % Delete one line of nodes along the span
%Now you're left with SWnodes containing 1 node coresponding to each
% panel (left leading point(double check if any doubt))

% Panel [number, nodes index of 4 nodes]
PanelData = [panelIndex' SWnodes' SWnodes'+N+C+1 SWnodes'+1 SWnodes'+N+C+2];


%% winglets
% Same thing for winglets
num_of_panels2 = (Nz-1)*(N+C);   
panelIndex2 = [1:num_of_panels2] + num_of_panels;
panelIndex3 = panelIndex2 + num_of_panels2;

SWnodes2 = NodeData2(:,1);
SWnodes3 = NodeData3(:,1);


SWnodes2(end - (N+C): end) = [];
Isw = [1:Nz-1]*(N+C+1);
SWnodes2(Isw) = [];

SWnodes3(end - (N+C): end) = [];
Isw = [1:Nz-1]*(N+C+1);
SWnodes3(Isw) = [];

PanelData2 = [panelIndex2' SWnodes2 SWnodes2+N+C+1 SWnodes2+1 SWnodes2+N+C+2];
PanelData3 = [panelIndex3' SWnodes3 SWnodes3+N+C+1 SWnodes3+1 SWnodes3+N+C+2];
PanelData = [PanelData; PanelData2; PanelData3];  

%% get normal vectors for the panels

% For wings, 1 for pannels containing 
n_hat_w = [ones(num_of_panels,1); zeros(2*num_of_panels2,1)];
n_hat_wl = [zeros(num_of_panels,1); -ones(2*num_of_panels2,1)];

%% panel lengths, areas and center points
Npanel = num_of_panels + (2*num_of_panels2); % Total number of panels

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
    
    %Area of each panel
    PAreas(ii) = ((x3-x1)+(x4-x2))*(((y2_l-y1_l)/2) + ((z2-z1)/2));
    
    % aerogrid contians the coordinate of centroid of each panel
    aerogrid(2,ii) = (y1_l+y2_l)/2;
    aerogrid(3,ii) = (z1+z2)/2;
    
    % As panels are either horizontal or verticle, therefore span of a
    % panel can be calculated as follows.
    PSpan(ii) = ((y2_l - y1_l) + (z2 - z1));
    
    
    xr1=x1+(x3-x1)/2;   % 1/2chord x-location at the root of the panel
    xr2=x2+(x4-x2)/2;   % 1/2chord x-location at the tip of the panel
    
    % 1/2chord x-location at the half-span location of the panel
    aerogrid(1,ii) = 0.5*(xr1+xr2); 
end

%plot3(aerogrid(1,:),aerogrid(2,:),aerogrid(3,:),'bx');

%% create area matrix (Skj)
S = zeros(6*Npanel,Npanel);
rows1 = 6*[1:Npanel] - 3;
rows2 = 6*[1:Npanel] - 1;
S(rows1,:) = diag(PAreas(1,:));
S(rows2,:) = diag(PAreas.*Plength*0.25);

% S matrix contains aeras and moment of areas about(quarter chord point)
% in rows, (3,9,15, ... ) a and (5,11,17,...) etc
S = sparse(S);
%% create control surface definitions and mode shapes

% For 8 control surfaces and 6dof for each panel
phiCS = zeros(Npanel*6,8);  
PIndex = PanelData(:,1); % Index of each panel

csIndex1 = find(mod(PIndex,N+C)>N); 
csIndex2 = find(mod(PIndex,N+C)==0);

% Contains panel index of CS, note that this includes the last panels on
% winglets but they are ignored later on.
csIndex = sort([csIndex1' csIndex2']);


% Binary variable, contains 1 in i,j position if ith panel is part of
% jth control surface
CSind = zeros(Npanel,8);

for ii = 1:length(csIndex) % For each CS panel
    Pno = csIndex(ii); % Current panel number
    ymid = aerogrid(2,Pno); % ycord of current CS panel midpoint
    
    % CSpanel markes the row number of the panel in CS
    CSpanel = mod(Pno,N+C);
    
    % Ctarting with 1 at CS root
    if CSpanel == 0
       CSpanel = C;
    else
       CSpanel = CSpanel - N; 
    end
    
    % Calculates the mode shapes associated with CS deflections for
    % body flaps
    if (abs(ymid)<0.36 && abs(ymid)>=0.0254)
        if ymid<0 % check if on the left half or right half
            % Calculates mode shape for 1 rad deflection, 
            phiCS((Pno*6 - 5):Pno*6,1) = ...
                [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3651; 0.931;0];
            % Assign CSind as control surface 1
            CSind(Pno,1) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,5) = [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3651; 0.931;0];
            CSind(Pno,5) = 1;
        end
    end
    
    % Same for other CS
    if (abs(ymid)<0.74 && abs(ymid)>=0.38)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,2) = [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3773; 0.9261;0];
            CSind(Pno,2) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,6) = [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3773; 0.9261;0];
            CSind(Pno,6) = 1;
        end
    end
            

    if (abs(ymid)<1.10 && abs(ymid)>=0.74)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,3) = [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3773; 0.9261;0];
            CSind(Pno,3) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,7) = [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3773; 0.9261;0];
            CSind(Pno,7) = 1;
        end
    end
    
    
    if (abs(ymid)<1.46 && abs(ymid)>=1.10)
        if ymid<0
            phiCS((Pno*6 - 5):Pno*6,4) = [0;0;-(CSpanel-0.5)*Plength(Pno); -0.3773; 0.9261;0];
            CSind(Pno,4) = 1;
        else
            phiCS((Pno*6 - 5):Pno*6,8) = [0;0;-(CSpanel-0.5)*Plength(Pno); 0.3773; 0.9261;0];
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
griddata.n_hat_w = n_hat_w;
griddata.n_hat_wl = n_hat_wl;
griddata.CSPanel = CSind;

save([saveloc,filename],'griddata');
