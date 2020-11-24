function [coords_spline,connections] = Geri_SplineGrid(coords_fem,BodySpan) 
%spline_grid   Creates a spline grid from the structural grid
%   Creates a spline grid by adding nodes forward and aft of some existing
%   FEM nodes. Note that the this function assumes that the FEM nodes lie
%   in the x-y plane. The resulting spline grid lies in the x-y plane as
%   well.
%
%   Inputs:
%       - coords_fem: A 3-by-numnodes matrix where the i-th column contains
%               [x;y;z] coordinate of the i-th node. The coordinate frame 
%               is such that:
%               A) Origin at the vehicle nose, 
%               B) +x pointing to rear or aircraft, and 
%               C) +y pointing starboard (to the right when looking
%                   forward)
%               D) +z pointing upward 
%       - BodySpan: Span of the centrebody defined by the outer edge of the 
%               body flaps
%
%   Outputs:
%       - coords_spline: A 3-by-NumSplineNodes matrix where i-th column
%               contains [x;y;z] coordinates of the spline grid. The spline
%               grids are defined using same coordiante frame as the FEM
%               nodes
%       - connections: A 1-by-NumSplineNodes vector where connections(i)=j
%               implies that the i-th spline node is derived From the
%               j-th FEM node. For FEM nodes, connection(i) = i;

%% Setup
displayfigs = false;     % Option to toggle figure display

%% Coordinates of FEM grid points
% coords is a 2-by-numnodes matrix where the i-th column contains [x;y] 
% coordinate of the i-th node. The coordinate frame is such that:
% A) Origin at the vehicle nose, 
% B) +x pointing to rear or aircraft, and 
% C) +y pointing starboard (to the right when looking forward)
% D) +z pointing upward 
xp = coords_fem(1,:);
yp = coords_fem(2,:);
zp = coords_fem(3,:);
Np = numel(xp);

% Plot FEM nodes
if displayfigs
    figure
    plot3(xp,yp,zp,'*')
    view(110,60)
    labels = cellstr(num2str((1:Np)'));
    grid on
    text(xp,yp,labels,'VerticalAlignment','bottom', ...
                      'HorizontalAlignment','right')
end

% Check if FEM nodes lie in x-y plane
if any(zp)
    error('FEM nodes should lie in the x-y plane')
end

%% Parameters
% Outlier Points: Do not make spline grid off of these gridpoints
% For Geri, these points include two nodes of the massless element
% connecting the servo to the spar for 6 wing control surfaces. They also
% include a nodes which are very close to other nodes.
nosplineidx = [52,51,7,50,43,44,45];

% Spline distance, m
% The spline grid adds points +/- sdist in front/behind (along the x
% direction) FEM grid points along the wing. Points listed in nosplineidx
% are not included for making these spline grid points.
sdist = 1e-3; 

%% Discard outlier and center body points
% Find centerbody points
cbidx = find(abs(yp)<=BodySpan/2);

% Remove centerbody and nospline points
idx = setdiff(1:Np, [nosplineidx(:); cbidx(:)]);

%% Spline Grid
% Aft Points
xaft = xp(idx) + sdist;
yaft = yp(idx);
zaft = zp(idx);
connectaft = idx;

% Forward Points
xfor = xp(idx) - sdist;
yfor = yp(idx);
zfor = zp(idx);
connectfor = idx;

% Spline Grid: FEM grid + Aft/Forward grid points
coords_spline = [xp xaft xfor;
                 yp yaft yfor;
                 zp zaft zfor]; 

% connections(i) = j implies that the i-th spline node is created using the
% j-th FEM node. For FEM nodes, connection(i) = i;             
connections = [1:Np, connectaft, connectfor]; 

end