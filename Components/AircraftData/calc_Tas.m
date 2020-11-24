function Tsa = calc_Tas(coords_spline, aerogrid)
%calc_Tsa: Calculates the transformation from spline to aerodynamic grid
%   Calculates the transformation matrix which calculates the 
%   [heave,bend,twist] of the aerodynamic nodes from the heave of spline 
%   nodes by interpolation. This routine assumes that the spline nodes lie 
%   in the x-y plane. ([heave,bend,twist]) of the x-y projection of the 
%   aerodynamic nodes are calculated by interpolation using radial basis 
%   function on the assumed flat plate (x-y plane). ([heave,bend,twist]) of
%   the aerodynamic nodes are then calculated by assuming that the nodes 
%   are connected to their x-y plane projection by stiff rods
%
%   Reference for the radial basis interpolation is:
%   'MSC/NASTRAN Version 68, Aeroelastic Analysis User's Guide (Page 43)'
% 
%   Inputs:
%     - coords_spline: 3-by-numSplineNodes matrix where i-th column 
%           contains [x;y;z] coordinates of the i-th spline node
%     - aerogrid: NumAeroNodes-by-3 matrix where i-th row contains [x,y,z]
%           coordinates of the i-th aerodynamic node
%       The coordinate frame for both splinegrid and aerogrid is such that
%           A) Origin at the vehicle nose, 
%           B) +x pointing to rear or aircraft, and 
%           C) +y pointing starboard (to the right when looking forward)
%           D) +z pointing upward 
%
%   Output:
%     - Tsa: A (2*NumAeroNodes)-by-(NumSplineNodes) matrix which when
%           multiplied by the heaves of spline nodes outputs [Heave,Twist] 
%           of the aerodynamic nodes

%% Note about the stiffness term
% Note that there is a 1/pi/16/D term which is missing from this code in
% the calculations involving the radial spline evaluations K and Ka.
% That term includes D which is the proportionality matrix for Hook's law
% for a 2-dim plate with plane stress condition. The equations implemented
% below can be thought of as choosing D such that 1/pi/16/D = 1;
% Tests are needed to makes sure that changing D doesn't effect the
% dynamics too much.

%% Setup
% Spline Data
Nspline = size(coords_spline,2);    % # of splinegrid nodes
xspline = coords_spline(1,:);       % x-coordinates of the nodes as ROW vector
yspline = coords_spline(2,:);       % y-coordinates of the nodes as ROW vector
zspline = coords_spline(3,:);       % z-coordinates of the spline nodes

if any(zspline)
    error('Spline nodes should lie in the x-y plane')
end

% Aero grid data
Naero = size(aerogrid,2);       % # of aerogrid nodes
xaero = aerogrid(1,:);          % x-coordinates of the nodes as COL vector
yaero = aerogrid(2,:);          % y-coordinates of the nodes as COL vector
zaero = aerogrid(3,:);          % z-coordinates of the nodes as COL vector

%% Radial distances
% Square of distances between each pair of spline nodes
% R is a NumSplineNodes-by-NumSplineNodes symmetric matrix where R(i,j)=p
% implies that the square of distance between i-th and j-th spline node is
% p meters^2
R = (xspline' - xspline).^2 + (yspline' - yspline).^2;

% Square of distances between each pair aero and spline nodes
% Ra is a NumAeroNodes-by-NumSplineNodes matrix where Ra(i,j)=p implies 
% that the square of distance between i-th aero node and j-th spline node 
% is p meters^2
Ra = (xaero'-xspline).^2 + (yaero'-yspline).^2;

%% Radial basis functions
% Radial basis function (Evaluate at pair of spline points)
% Note: As noted before, 1/(16*pi*D) = 1
K = R.*log(R);
K(isnan(K)) = 0;

% Radial basis function (Evaluate at pair of aero and spline points)
% Note: As noted before, 1/(16*pi*D) = 1
Ka = Ra.*log(Ra);
Ka(isnan(Ka)) = 0;

%% Interpolation matrix for heave of aerodynamic nodes
% Equation (2.40) in Ref:
% [0;0;0; Heave at Spline 1; ... Heave at Spline N] = C*P
% where P =[a0; a1; a2; P1; ...; PN] are the unknown coefficients in the
% radial spline fit. Pi is interpreted as a point load at spline point i.
S = [ones(1,Nspline); xspline; yspline];
C = [zeros(3,3), S; S', K];

% Equation (2.41) in Ref: 
% Matrix Mh2h maps from spline heave to aero heave, i.e.
% [Heave at Aero 1; ... ; Heave at Aero N] 
%      = Mh2h * [Heave at Spline 1; ... Heave at Spline N] 
invC = inv(C);
invC = invC(:,4:end);
Ca = [ones(Naero,1), xaero' yaero' Ka];
Mh2h = Ca*invC;

%% Interpolation matrix for twist of aerodynamic nodes
% Note: Twist is the negative of the partial of heave with respect to x
% Equations (2.42-2.43) in Ref: 
% Matrix Mh2t maps from spline heave to the negative aero twist, i.e.
% [-Twist at Aero 1; ... ; -Twist at Aero N] 
%      = Mh2h * [Heave at Spline 1; ... Heave at Spline N] 
deltaxa = xaero'-xspline;
dKat = 2*deltaxa.*(1+log(Ra));
dKat(isnan(dKat)) = 0;

dCat = [zeros(Naero,1),ones(Naero,1),zeros(Naero,1),dKat];
Mh2t = -dCat*invC;

%% Interpolation matrix for bend of aerodynamic nodes
% Note: Bend is the partial of heave with respect to y
% Equations (2.42-2.43) in Ref modified for bend: 
% Matrix Mh2b maps from spline heave to the , i.e.
% [Bend at Aero 1; ... ; Bend at Aero N] 
%      = Mh2h * [Heave at Spline 1; ... Heave at Spline N] 
deltaya = yaero'-yspline;
dKab = 2*deltaya.*(1+log(Ra));
dKab(isnan(dKab)) = 0;

dCab = [zeros(Naero,1),ones(Naero,1),zeros(Naero,1),dKab];
Mh2b = dCab*invC;

%% Combine heave and twist calculation
% Compute matrix with heave of all nodes followed by twist of all nodes
% Tsa_tmp = [zeros(Naero,Nspline); zeros(Naero,Nspline); Mh2h; ...
%             Mh2b; Mh2t; zeros(Naero,Nspline)];

% Compute transfermotion from projection to the actual nodes
Tsa_tmp = [zaero'.*Mh2t; -zaero'.*Mh2b; Mh2h; ...
            Mh2b; Mh2t; zeros(Naero,Nspline)];
        
% Rearrange with [Heave,Twist] of each node one after another
Tsa = zeros(size(Tsa_tmp));
Tsa(1:6:end,:) = Tsa_tmp(1:Naero,:);
Tsa(2:6:end,:) = Tsa_tmp(Naero+1:2*Naero,:);
Tsa(3:6:end,:) = Tsa_tmp(2*Naero+1:3*Naero,:);
Tsa(4:6:end,:) = Tsa_tmp(3*Naero+1:4*Naero,:);
Tsa(5:6:end,:) = Tsa_tmp(4*Naero+1:5*Naero,:);
Tsa(6:6:end,:) = Tsa_tmp(5*Naero+1:6*Naero,:);

end