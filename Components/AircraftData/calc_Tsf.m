function Tfs = calc_Tsf(coords_fem,coords_spline,connections)
%calc_Tfs   Calculates the transformation matrix from FEM to spline grid 
%   Mode shapes known for the FEM grid have degrees of freedom
%   [Heave, Roll, Pitch]. This function creates a transformation matrix
%   which when multiplied by the FEM mode shape, outputs the corresponding 
%   heave deflection of the spline nodes. Note that it is assumed that both
%   FEM nodes and spline nodes lie in the x-y plane
% 
%   Inputs:
%       - coords_fem: 3-by-NumFEMNodes matrix where i-th column contains
%           [x;y;z] coordinates of the i-th FEM node. The coordinate frame is
%           such that:
%           A) Origin at the vehicle nose, 
%           B) +x pointing to rear or aircraft, and 
%           C) +y pointing starboard (to the right when looking forward)
%           D) +z pointing upward 
%       - coords_spline: 3-by-numSplineNodes matrix where i-th column 
%           contains [x;y;z] coordinates of the i-th spline node
%       - connections: A 1-by-NumSplineNodes vector where connections(i)=j
%           implies that i-th spline node is derived from the j-th FEM
%           node. For FEM nodes, connection(i) = i
% 
%   Output:
%       - Tfs: An NumSplineNodes-by-(3*NumFEMNodes) matrix which when
%           multiplied by the FEM mode shape (containing [Have,Roll,Pitch] 
%           of the FEM nodes) gives the resulting heave deflections of the 
%           spline nodes

%% Setup and checks
NumSplineNodes = size(coords_spline,2);    % Number of splinegrid nodes
xspline = coords_spline(1,:);              % x-coordinates of the spline nodes
yspline = coords_spline(2,:);              % y-coordinates of the spline nodes
zspline = coords_spline(3,:);              % z-coordinates of the spline nodes

NumFEMNodes = size(coords_fem,2);        % Number of FEM nodes
xstruc = coords_fem(1,:);                % x-coordinates of the FEM nodes
ystruc = coords_fem(2,:);                % y-coordinates of the FEM nodes
zstruc = coords_fem(3,:);                % z-coordinates of the FEM nodes

% Check if FEM noes and spline nodes lie in x-y plane
if any(zstruc)
    error('FEM nodes should lie in the x-y plane')
end

if any(zspline)
    error('Spline nodes should lie in the x-y plane')
end

%% Calculate transformation matrix
% X,Y displacements between spline nodes and corresponding FEM nodes
rx = xspline-xstruc(connections);
ry = yspline-ystruc(connections);

% Tfs multiplied by the FEM mode shape  (containing [Have,Roll,Pitch] of 
% the FEM nodes) gives the resulting heave deflections of the spline nodes
Tfs = zeros(NumSplineNodes,6*NumFEMNodes); 
for ii = 1:NumSplineNodes
    % Index to FEM connection point
    jj = connections(ii);
    
    % Effect of [surge,sway,heave,bend,twist,yaw] DOFs
    Tfs(ii,(6*jj-5):(6*jj)) = [0 0 1 ry(ii) -rx(ii) 0];
end

end