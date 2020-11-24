function DLMData = getUnsteadyAero(NumPanels,PanelData,NodeData, ...
                                    PanelSpan,PanelLength,MAC,D0)
%% getDLMAIC Obtains AIC matrix using Doublet Lattice Method 
%
%   Input:
%       - NumPanels: Total number of panels in the aerodynamic grid
%       - PanelData: A 4-by-NumPanels matrix where i-th column contains
%           the NW, NE, SW, SE node indexes of the i-th panel
%       - NodeData: A 3-by-NumNodes matrix where i-th column contains
%           [x;y;z] coordinates of the i-th node
%       - PanelSpan: A 1-by-NumPanels vector where i-th element is the
%           span (distance between two edges parallel to the x-axis)
%       - PanelLength: A 1-by-NumPanels vector where i-th element is 
%           the average length of the panel in x-direction
%       - MAC: Mean Aerodynamic chord of the aircraft used as reference
%           length
%       - D0: A (3*NumPanels)-by-NumPanels matrix which when multiplied
%               to the vortex strength vector (Gamma) results in 3DOF wash
%               at the collocation point
%
%   Output:
%       - DLMData: A structure containing the rational function
%           approximation of the AIC in the following fields
%           - A0,A1,Ap: NumPanels-by-NumPanels matrices such that 
%               AIC = A0 + A1s + Ap*(s)/(s+p)
%           - Alag,Blag,Clag: State space representation of Ap/(s+p)

%% Define frequency grid for DLM calculation
% The frequency domain dependancy of AIC and rational function
% approximation needs to be described properly.
% Cp(\omega){complex} = AIC(omega,V)*w/V(omega){complex}
% Now, note that 

% Frequencies = (cref*omega)/(2*U_Infinity)
k = [0 0.0469 0.0938 0.1875 0.375 0.75 1.5 3];

% Reduced frequencies = (Omega/U_Infinity)
kr = (2/MAC)*k;

%% Extract panel information

% Grab indices to panel nodes.
% Note: Viewing from top, if LE is toward north, left is towards west. 
% The panel data is [panel number, NW, NE, SW, SE nodes].
NWidx = PanelData(1,:);
NEidx = PanelData(2,:);
SWidx = PanelData(3,:);
SEidx = PanelData(4,:);

% Define downwash location (aka collocation points)
% The collocation points are defined at the 3/4 chord and half span of 
% the aero panels.

% Collocation points: x coordinate
% Computed as the average of the 3/4 chord point each side of the panel
x1 = 0.25*NodeData(1,NWidx) + 0.75*NodeData(1,SWidx);
x2 = 0.25*NodeData(1,NEidx) + 0.75*NodeData(1,SEidx);
P0x = (x1 + x2)/2;
P0y = (NodeData(2,NWidx) + NodeData(2,NEidx))/2;    % y coordinate
P0z = (NodeData(3,NWidx) + NodeData(3,NEidx))/2;    % z coordinate

% Define doublet / vortex locations 
% These are defined at the 1/4 chord and at the edges (0 and full span)
% of each aerodynamic panel.
% P1 contains points at the 1/4 chord and zero span of each panel
% The i-th row contains the (x,y,z) coordinate of the point.
P1x = NodeData(1,NWidx) + (NodeData(1,SWidx)-NodeData(1,NWidx))/4;
P1y = NodeData(2,NWidx);
P1z = NodeData(3,NWidx);

% P3 contains points at the 1/4 chord and full span
% The i-th row contains the (x,y,z) coordinate of the point.
P3x = NodeData(1,NEidx)+ (NodeData(1,SEidx)-NodeData(1,NEidx))/4; %xp3
P3y = NodeData(2,NEidx); %yp3
P3z = NodeData(3,NEidx); %zp3

% P2 contains the point at 1/4 chord and half-span location of the panel
P2x = (P1x+P3x)/2;
P2y = (P1y+P3y)/2;
P2z = (P1z+P3z)/2;

% root doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x01= (P0x-P1x')';
y01= (P0y-P1y')';
z01= (P0z-P1z')';

% semi-span doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x02= (P0x-P2x')';
y02= (P0y-P2y')';
z02= (P0z-P2z')';

% tip doublet
% define the x,y,z distance from this doublet to the recieving point at
% 1/2span 3/4chord location
x03= (P0x-P3x')';
y03= (P0y-P3y')';
z03= (P0z-P3z')';

% cos and sin angles of all doublet line dihedrals
% CosGamma = y/sqrt(z^2 + y^2)
cosGamma = (P3y - P1y)./(sqrt((P3z-P1z).^2 + (P3y-P1y).^2));
cosGamma = cosGamma';

% SinGamma = z/sqrt(z^2 + y^2)
sinGamma = (P3z - P1z)./(sqrt((P3z-P1z).^2 + (P3y-P1y).^2));
sinGamma = sinGamma';

DLMInfo.x01 = x01;
DLMInfo.y01 = y01;
DLMInfo.z01 = z01;
DLMInfo.x02 = x02;
DLMInfo.y02 = y02;
DLMInfo.z02 = z02;
DLMInfo.x03 = x03;
DLMInfo.y03 = y03;
DLMInfo.z03 = z03;
DLMInfo.cosGamma = cosGamma;
DLMInfo.sinGamma = sinGamma;

Mach = 0;
Kp_0.Ki = getKappa(x01,y01,z01,cosGamma,sinGamma,0,Mach);
Kp_0.Km = getKappa(x02,y02,z02,cosGamma,sinGamma,0,Mach);
Kp_0.K0 = getKappa(x03,y03,z03,cosGamma,sinGamma,0,Mach);

%% Obtain AIC
AIC = zeros(NumPanels,NumPanels,length(k));
for i = 1:length(k)
    AIC(:,:,i) = getDLM(NumPanels,DLMInfo,PanelSpan/2,PanelLength,kr(i),Mach,Kp_0,D0);
    fprintf('DLM calculations completed with kr = %4.3f\n',kr(i))
end

%% Rational function approximation
b_poles = [0.1100 0.2200];    % lag poles
[A_0,A_1,A_p,Alag,Blag,Clag] = ...
                    rational_approx(AIC,b_poles,k);

%% Check 
% for i = 1:length(k)
%     err(i) = norm(AIC(:,:,i) - (A_0 + A_1*1j*k(i) ...
%                 + A_p(:,:,1) * (1j*k(i))/(1j*k(i) + b_poles(1)) ...
%                 + A_p(:,:,2) * (1j*k(i))/(1j*k(i) + b_poles(2)) )) ...
%             /norm(AIC(:,:,i));
% end
% err

%% Output rational approximation of the AIC
DLMData.A0 = A_0;
DLMData.A1 = A_1;
DLMData.Ap = A_p;
DLMData.Alag = Alag;
DLMData.Blag = Blag;
DLMData.Clag = Clag;

end