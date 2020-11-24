function VLMData = getSteadyAero(NumPanels,PanelData,NodeData,N_hat0, ...
                                    N_hat_eta,rcg)
%% getSteadyAero  Obtain steady aerodynamic model
%
%   Obtain steady aerodynamic model based of vortex lattice method
%
%   Input:
%       - NumPanels: Total number of panels in the aerodynamic grid
%       - PanelData: A 4-by-NumPanels matrix where i-th column contains
%           the NW, NE, SW, SE node indexes of the i-th panel
%       - NodeData: A 3-by-NumNodes matrix where i-th column contains
%           [x;y;z] coordinates of the i-th node
%       - N_hat0: A 3-by-NumPanels vector where i-th column contains the
%           direction cosinse of the panel normal vector
%       - Nhat_eta: A 3-by-NumPanels-by-NumFlexModes matrix where
%           Nhat_eta(:,:,i) contains the the depenency of the panel normals
%           on the flexible deflections eta:
%           Nhat = Nhat_0 + Nhat_eta(:,:,1)*eta1 + Nhat_eta(:,:,2)*eta2 ...
%       - Rcg: Position vector of the structural node corresponing to the 
%           CG of the aircraft
%
%   Output:
%       - VLMData: A structure containing the following fields:
%           - D0: A (3*NumPanels)-by-NumPanels matrix which when multiplied
%               to the vortex strength vector (Gamma) results in 3DOF wash
%               at the collocation point
%           - Ainv0: An NumPanels-by-NumPanels matrix containing the AIC
%               obtained for zero flexible deflection
%           - Ainv_eta: A structure containing the following fields
%               - Ainv_eta_L: A NumFlexModes element cell where each cell
%                   contains NumPanel-by-x matrices
%               - Ainv_eta_L: A NumFlexModes element cell where each cell
%                   contains x-by-NumPanel matrices
%               such that AIC = Ainv0 + Ainv_eta_L(1)*Ainv_eta_R(1)*eta1 +
%               ...
%           - D_Vortex: A (3*NumPanels)-by-NumPanels matrix which when
%               multiplied to the vortex strength vector (Gamma) results in
%               3DOF wash at the vortex location of the panels
%           - r_cgcolloc: A 3-by-NumPanels matrix contining  the position
%               vector of the collocation points with origin as the CG of
%               the aircraft
%           - r_cgvortex: A 3-by-NumPanels matrix contining  the position
%               vector of the mid-point of the quarter chord line with 
%               origin as the CG of the aircraft
%           - r_colloc: A 3-by-NumPanels matrix contining  the position
%               vector of the collocation points with origin as the nose of
%               the aircraft
%           - r_vortex: A 3-by-NumPanels matrix contining  the position
%               vector of the mid-point of the quarter chord line with 
%               origin as the nose of the aircraft
%           - Gamma_len: A 3-by-NumPanels matrix where ith column contains
%               the vector joining the two ends of the finite part of the
%               horseshoe vortex

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

%% Calculate arm lengths
r_cgcolloc = P0 - rcg;
r_cgvortex = P2 - rcg;
r_colloc = P0;
r_vortex = P2;

%% Calculate vortex direction
Gamma_len = P3-P1;

%% Calculate linear approximation to AIC matrix
% AIC matrix is approximated as a Taylor series expansion as:
%
% AIC = A0inv + A1inv(:,:,1)*eta1 + A1inv(:,:,2)*eta2 .......
% This approximation requires that the modal deformation should be small

% Induced velocity due to vortices at collocation point
D_colloc = get6dofVLM(P0,P1,P3);

% For zero eta, the AIC is calculated as
D0 = reshapefordot(N_hat0,NumPanels) * D_colloc;
Ainv0 = inv(D0);

% NumModes = size(N_hat_eta,3);
% for i = 1:NumModes                    
% %     % Changes in normalwash 
%     A1 = reshapefordot(N_hat_eta(:,:,i),NumPanels)*D_colloc;
%     % Linear approximation of change in normalwash
%     A1inv_temp = -Ainv0*A1*Ainv0;
% 
%     % For faster processing in SIMULINK, split the matrices using SVD
%     [U,S,V] = svd(A1inv_temp);
%     
%     % Find index(j) where SVD reduces by factor of 100 and split there
%     [~,j] = min(abs(diag(S)-S(1,1)/100));
%     Srt = sqrt(S(1:j,1:j));
%     Ainv_eta.Ainv_eta_R{i} = Srt*V(:,1:j)';
%     Ainv_eta.Ainv_eta_L{i} = U(:,1:j)*Srt;
% end

%% Find induces velocity at vortex location due to Gamma
D_Vortex = get6dofVLM(P2,P1,P3);

%% Output VLMData and Grid Information
VLMData.D0 = D0;
VLMData.Ainv0 = Ainv0;
% VLMData.Ainv_eta = Ainv_eta;
VLMData.D_Vortex = D_Vortex;
VLMData.grid.r_cgcolloc = r_cgcolloc;
VLMData.grid.r_cgvortex = r_cgvortex;
VLMData.grid.r_colloc = r_colloc;
VLMData.grid.r_vortex = r_vortex;
VLMData.grid.Gamma_len = Gamma_len;

end

%% Function for dot product calculation
function [Nout] = reshapefordot(Nin,NumPanels)
    Nout = zeros(NumPanels,3*NumPanels);
    for i = 1:NumPanels
        Nout(i,(3*i-2:3*i)) = Nin(:,i)';
    end
end