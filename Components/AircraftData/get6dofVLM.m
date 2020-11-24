function D = get6dofVLM(PC,PA,PB)
% get6dofVLM  Obtains vortex lattice method based aerodynamic model
%
%   Inputs:
%       - PC: 3-by-NumPanels matrix of the collocation points. The i^th 
%           column of PC gives the (x,y,z) coordinations of the i^th 
%           collocation point
%
%       - PA,PB: 3-by-NumPanels matrices containing the points describing 
%           the horseshoe vortices for each panel. The i^th column of each 
%           matrix defines the (x,y,z) coordinates of the i^th point. The 
%           ith horseshoe vortex consists of:
%               *A semi-infinite vortex line from infinity to PA(:,i) along
%               the streamwise direction.
%               *A finite vortex line from PA(:,i) to PB(:,i) 
%               *A semi-infinite vortex line from PB(:,i) to infinity along
%               the streamwise direction.
%           See Figure 6-9 in the reference below. The unit is meters%
%
% Output:
%   - D: A (3*NumPanels)-by-NumPanels matrix which when multiplied by the
%           vector of strength of the horseshoe vortices, outputs the total
%           was at the quarter chord point
%
% Reference: 
% Chapter 6 entitled "Aerodynamics of 3D Lifting Surfaces through Vortex 
% Lattice Methods" in "Applied Computational Aerodynamics" by Cummings, 
% Mason, Morton and McDaniel.

%% Number of Panels
NumPanels = size(PC,2);

%% Threshold
% This threshold is used to determine if a collocation point is too
% close to the vortex or is collinear with the vortex direction.
epsilon = 10e-6;

%% Pre-compute vectors 
% See Figure 6.10 in the reference for diagram:
%  r0 := Vector from PA to PB
%  r1 := Vector from PA to PC 
%  r2 := Vector from PB to PC
% where PA and PB are the starting and ending points of the finite part of
% the horseshoe. PC are the collocation points.

% Compute r1:= Vector from PA to PC
% This gives the x-distance between every pair of vortex points PA and
% collocation points P1. The same is true for the y and z coords below.
r1x = PC(1,:)' - PA(1,:);
r1y = PC(2,:)' - PA(2,:);
r1z = PC(3,:)' - PA(3,:);

% Compute r2 := Vector from PB to PC
r2x = PC(1,:)' - PB(1,:);
r2y = PC(2,:)' - PB(2,:);
r2z = PC(3,:)' - PB(3,:);

% Compute r0 := Vector from PA to PB
% The code below computes r0x as a NumPanels-by-NumPanels matrix
% with each row identically the same, i.e.
%  r0x(i,j) = PB(j,1)-PA(j,1)
% This aligns with the matrices (r1x,r2x,etc) computed where
% the row "i" dimension corresponds to the collocation point. 
% This allows for more efficient vectorized code below.
r0x = repmat(PB(1,:) - PA(1,:), [NumPanels 1]);
r0y = repmat(PB(2,:) - PA(2,:), [NumPanels 1]);
r0z = repmat(PB(3,:) - PA(3,:), [NumPanels 1]);

%% Induced velocity due to finite vortex line from A to B
% This computation is based on Equation 6.48 in the reference. The result 
% is normalized by the vortex strength, i.e. it is the induced velocity
% divided by the vortex strength.

% Compute (x,y,z) components of r1 x r2
r1Xr2_x = (r1y.*r2z) - (r2y.*r1z);
r1Xr2_y = -(r1x.*r2z) + (r2x.*r1z);
r1Xr2_z = (r1x.*r2y) - (r2x.*r1y);

% Compute magnitudes of r1, r2, and r1 x r2
mag_r1Xr2 = sqrt((r1Xr2_x.^2)+(r1Xr2_y.^2)+(r1Xr2_z.^2));
mag_r1 = sqrt((r1x.^2) + (r1y.^2) + (r1z.^2));
mag_r2 = sqrt((r2x.^2) + (r2y.^2) + (r2z.^2));

% Compute dot (inner) products of <r0,r1> and <r0,r2>
r0r1 = (r0x.*r1x) + (r0y.*r1y) + (r0z.*r1z);
r0r2 = (r0x.*r2x) + (r0y.*r2y) + (r0z.*r2z);

% (D1_u,D1_v,D1_w) are NumPanels-by-NumPanels matrices of the induced 
% velocity components in the (x,y,z) directions. The (i,j) entry of each
% matrix provides the velocity at the i^th collocation point due to the 
% finite horseshoe segment on the j^th panel.  These are induced 
% velocities divided by the vortex strength.
D1_base = ( (r0r1./mag_r1) - (r0r2./mag_r2) ) ./ (mag_r1Xr2.^2)/(4*pi);
D1_u = r1Xr2_x.*D1_base;
D1_v = r1Xr2_y.*D1_base;
D1_w = r1Xr2_z.*D1_base;

% Correction: Adjust the D1 matrix if the collocation point is too 
% close to the vortex or is collinear with the vortex direction.
ind1 =  find( mag_r1<epsilon );
ind2 =  find( mag_r2<epsilon );
ind3 =  find( mag_r1Xr2<epsilon );
ind = unique( [ind1; ind2; ind3] );

D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

% Get final D1 matrix 
% D1 matrix contains the perpendicular component of induced velocities 
% at all panels. For wing panels, it's the z component of induced 
% velocities (D1_w) while for winglets, it's the y component of induced 
% velocities (D1_v)
% D1 = (D1_w.*n_hat_w) + (D1_v.*n_hat_wl);

%% Induced velocity due to semi-infinite vortex line from wake to A 
% Perpendicular distance of the collocation point from the vortex line
d2 = sqrt((r1y.^2) + (r1z.^2));

% Angle between collocation point and the vortex vertices
cosBB1 = 1;             % From the wake at x = +infinity  
cosBB2 = -r1x./mag_r1;  % From the bound vertex

% Equation 6.43 in the attached reference
D2_full = (1/(4*pi))*(cosBB1 - cosBB2)./d2;

% Find x,y,z component of the velocity
gamma = atan2(r1z,r1y);
D2_u = zeros(size(D2_full));
D2_v = sin(gamma).*D2_full;
D2_w = -cos(gamma).*D2_full;

% Adjust the D2 matrix if the collocation point is too close to the vortex 
% Note that as the vortex extends to infinity, therefore no need to check
% the other edge or cross product
ind = find( mag_r1<epsilon );
D2_u(ind) = 0;
D2_v(ind) = 0;
D2_w(ind) = 0;

% Get final D2 matrix (same as D1)
% D2 = (D2_w.*n_hat_w) + (D2_v.*n_hat_wl);

%% Induced velocity due to semi-infinite vortex line from B to wake  

% Perpendicular distance of the collocation point from the vortex line
d3 = sqrt((r2y.^2) + (r2z.^2)); 

% Angle between collocation point and the vortex vertices
cosBT1 = r2x./mag_r2;   % From the bound vertex
cosBT2 = -1;            % From the wake at x = +infinity

% Equation 6.43 in the attached reference
D3_base = (1/(4*pi))*(cosBT1 - cosBT2)./d3;

% Find x,y,z component of the velocity
gamma = atan2(r2z,r2y);
D3_u = zeros(size(D3_base));
D3_v = -sin(gamma).*D3_base;
D3_w = cos(gamma).*D3_base;

% Adjust the D3 matrix if the collocation point is too close to the vortex 
% Note that as the vortex extends to infinity, therefore no need to check
% the other edge or cross product
ind = find( mag_r2<epsilon );
D3_u(ind) = 0;
D3_v(ind) = 0;
D3_w(ind) = 0;

%% Combine effects of horseshoe vortices
D_u = D1_u + D2_u + D3_u;
D_v = D1_v + D2_v + D3_v;
D_w = D1_w + D2_w + D3_w;

% Obtain the total influcence matrix
D = zeros(3*NumPanels,NumPanels);
D(1:3:end,:) = D_u;
D(2:3:end,:) = D_v;
D(3:3:end,:) = D_w;

end