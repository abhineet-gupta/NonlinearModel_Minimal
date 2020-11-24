function [offset_s,connections] = spline_grid(offset_g) 
% Creates a spline grid from the structural grid. 
%
% Input:
%       offset_g: Structural grid 
% Output:
%       offset_s: Spline grid
%       connections: Rigid connections between structural nodes and spline
%       grid.
%
%% Eliminate outlier points
% Do not make spline grid off of these gridpoints
Pt_discard = [52,30,51,28,8,50,26,43,15,44,45,3,17,19];

% Points on the FEM grid
xp = offset_g(1,:);
yp = offset_g(2,:);

% Number of FEM gridpoints
numpt=length(xp);

% Vector whos ith = j, where i is the index of spline grid point
% (which includes FEM grid points) and j is the correspoinding FEM
% gridpoint which is attached to the spline grid point.
% If i is in the FEM grid to begin with the connections(i) = i;
% Note that lenght of connections will change
connections = 1:numpt; 

j = 0;
% Wings grid
for i = 1:numpt
    if abs(yp(i))>0.36
        if ~ismember(i,Pt_discard)
            j=j+1;
            dx1(j) = xp(i) + 0.00001; 
            dy1(j) = yp(i);
            connections(numpt+j) = i;
        end
    end
end

k=0;
for i = 1:numpt
    if abs(yp(i))>0.36
        if ~ismember(i,Pt_discard)
            k=k+1;
            dx2(k) = xp(i) - 0.00001; % Spline_6
            dy2(k) = yp(i);
            connections(numpt+j+k) = i;
        end
    end
end

% winglet spline grid
offset_s = [xp dx1 dx2; yp dy1 dy2]; % Spline grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Body grid
% dx1(abs(yp) < 0.36) = xp(abs(yp) < 0.36) + 0.26;
% dx2(abs(yp) < 0.36) = xp(abs(yp) < 0.36) - 0.2895;


% dx1(abs(yp) < 0.36) = xp(abs(yp) < 0.36);
% dx2(abs(yp) < 0.36) = xp(abs(yp) < 0.36);
% 
% dy1 = yp;
% dy2 = yp;

% xp = -offset_g(1,:);
% yp = -offset_g(2,:);
% 
% % Wings grid
% % dx1(abs(yp) > 0.36) = xp(abs(yp) > 0.36) + 0.1616;
% % dx2(abs(yp) > 0.36) = xp(abs(yp) > 0.36) - 0.1344;
% 
% dx1(abs(yp) > 0.36) = xp(abs(yp) > 0.36) + 0.000001;
% dx2(abs(yp) > 0.36) = xp(abs(yp) > 0.36) - 0.000001;
% 
% % winglet spline grid
% 
% % Body grid
% % dx1(abs(yp) < 0.36) = xp(abs(yp) < 0.36) + 0.26;
% % dx2(abs(yp) < 0.36) = xp(abs(yp) < 0.36) - 0.2895;
% 
% 
% dx1(abs(yp) < 0.36) = xp(abs(yp) < 0.36) + 0.000001;
% dx2(abs(yp) < 0.36) = xp(abs(yp) < 0.36) - 0.000001;
% 
% dy1 = yp;
% dy2 = yp;
% 
% offset_s = [xp dx1 dx2; yp dy1 dy2];
% connections = repmat(1:length(xp),[1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%