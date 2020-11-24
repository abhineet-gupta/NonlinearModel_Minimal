function [AIC] = getAIC(griddata,kr)
%getAIC   Obtain the Aerodynamic Influence Coefficient (AIC)
%   Calculates the Aerodynamic Influcence Coefficient matrix using Vortex 
%   lattice method(VLM) or doublet lattice method(DLM)
%
%   Inputs:
%       - griddata: A structure containing the following fields
%           - NodeData: A 3-by-NumNodes matrix where i-th column contains
%               [x;y;z] coordinates of the i-th node
%           - PanelData: A 4-by-NumPanels matrix where i-th column contains
%               the NW, NE, SW, SE node indexes of the i-th panel
%           - NumNodes: Total number of nodes
%           - NumPanels: Total number of panels
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
%                the body flaps
%           - PhiCS: A (6*NumPanels)-by-NumCS matrix where 
%               PhiCS(6*i-5:6*i,k) contains [x,y,z,theta_x,theta_y,theta_z]
%               degrees of freedom for the i-the panel for 1 radian
%               deflection mode of k-th control surface.
%       - kr: A 1-by-NumOmega vector containing the reduced frequencies
%       - AeroType: Variable to define the aerodynamic model type, should 
%           be either 'VLM' or 'DLM'

%% Obtian grid information
PanelData = griddata.PanelData;
NodeData = griddata.NodeData;
PanelSpan = griddata.PanelSpan;
PanelLength = griddata.PanelLength;

% Determine number of aerodynamic panels
NumPanels =size(PanelData,2);

%% Define downwash location (aka collocation points)
% The collocation points are defined at the 3/4 chord and half span of 
% the aero panels.

% P0 is NumPanels-by-3.  The i-th row contains the (x,y,z) coordinate of 
% the collocation point of i-th panel.
P0 = zeros(3,NumPanels); 

% Grab indices to panel nodes.
% Note: Viewing from top, if LE is toward north, left is towards west. 
% The panel data is [panel number, NW, NE, SW, SE nodes].
NWidx = PanelData(1,:);
NEidx = PanelData(2,:);
SWidx = PanelData(3,:);
SEidx = PanelData(4,:);

% Collocation points: x coordinate
% Computed as the average of the 3/4 chord point each side of the panel
x1 = 0.25*NodeData(1,NWidx) + 0.75*NodeData(1,SWidx);
x2 = 0.25*NodeData(1,NEidx) + 0.75*NodeData(1,SEidx);
P0(1,:) = (x1 + x2)/2;

% Collocation points: y coordinate
P0(2,:) = (NodeData(2,NWidx) + NodeData(2,NEidx))/2; 

% Collocation points: z coordinate
P0(3,:) = (NodeData(3,NWidx) + NodeData(3,NEidx))/2; 

%% Define doublet / vortex locations 
% These are defined at the 1/4 chord and at the edges (0 and full span)
% of each aerodynamic panel.

% P1 contains points at the 1/4 chord and zero span of each panel
% The i-th row contains the (x,y,z) coordinate of the point.
P1 = zeros(3,NumPanels);

P1(1,:) = NodeData(1,NWidx) + (NodeData(1,SWidx)-NodeData(1,NWidx))/4; %xp1
P1(2,:) = NodeData(2,NWidx); %yp1
P1(3,:) = NodeData(3,NWidx); %zp1

% P3 contains points at the 1/4 chord and full span
% The i-th row contains the (x,y,z) coordinate of the point.
P3 = zeros(3,NumPanels);
P3(1,:) = NodeData(1,NEidx)+ (NodeData(1,SEidx)-NodeData(1,NEidx))/4; %xp3
P3(2,:) = NodeData(2,NEidx); %yp3
P3(3,:) = NodeData(3,NEidx); %zp3

% P2 contains the point at 1/4 chord and half-span location of the panel
P2 = zeros(3,NumPanels);
P2(1,:) = (P1(1,:)+P3(1,:))/2;
P2(2,:) = (P1(2,:)+P3(2,:))/2;
P2(3,:) = (P1(3,:)+P3(3,:))/2;

%% Get relevant inputs for calculating Kernel function 
%(Page 1, reference 1,beginning with eq. 2)

% root doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x01= bsxfun(@minus,P0(1,:),P1(1,:)');
y01= bsxfun(@minus,P0(:,2),P1(:,2)');
z01= bsxfun(@minus,P0(:,3),P1(:,3)');

% semi-span doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x02= bsxfun(@minus,P0(:,1),P2(:,1)');
y02= bsxfun(@minus,P0(:,2),P2(:,2)');
z02= bsxfun(@minus,P0(:,3),P2(:,3)');

% tip doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x03= bsxfun(@minus,P0(:,1),P3(:,1)');
y03= bsxfun(@minus,P0(:,2),P3(:,2)');
z03= bsxfun(@minus,P0(:,3),P3(:,3)');

%% cos and sin angles of all doublet line dihedrals

% CosGamma = y/sqrt(z^2 + y^2)
cosGamma = (P3(:,2) - P1(:,2))./ ...
            (sqrt((P3(:,3)-P1(:,3)).^2 + (P3(:,2)-P1(:,2)).^2));
            
% SinGamma = z/sqrt(z^2 + y^2)
sinGamma = (P3(:,3) - P1(:,3))./(sqrt((P3(:,3)-P1(:,3)).^2 + ...
            (P3(:,2)-P1(:,2)).^2));

%% Obtain AIC matrix
Dd = getDLM(P0',P1',P2',P3',PanelSpan'/2,PanelLength',kr,0);  % DLM (oscillatory effects)
AIC = inv(Dd);
end