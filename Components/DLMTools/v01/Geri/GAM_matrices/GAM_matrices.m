function [geom,modes] = GAM_matrices(griddata,cref,geom,modes)
%GAM matrices   Generates differentiation and transformation matrices
%   This script constructs the matrices required to carry out 
%   transformation of the FEM modes to the aerodynamic grid. It also 
%   calculates the differentiation matrices used to calculate downwash at 
%   the collocation points for a given deflection and velocity of the
%   aerodynamic panel
%
%   Inputs:
%       - griddata: A structure containing the following fields (might
%           contain other fields as well)
%           - Aeorgrid: An 3-by-NumPanels matrix where i-th column contains
%               the [x;y;z] coordinates of the centroid of the i-th panel.
%               The coordinate frame is such that
%               A) Origin at the vehicle nose, 
%               B) +x pointing to rear or aircraft, and 
%               C) +y pointing starboard (to the right when looking
%               forward)
%               D) +z pointing upward 
%           - PanelLength: A 1-by-NumPanels vector where i-th element is 
%               the average length of the panel in x-direction
%           - S: A (6*NumPanels)-by-NumPanels matrix where where blocks of 
%               six rows give the [Ax; Ay; Az; Mx; My; Mz] data.
%               Ax is the panel area projected perpendicular to the x-axis.
%               Ay and Az are defined similarly.
%               Mx is th moment of the area about the x-axis. My and Mz are
%               defined similarly. 
%           - BodySpan:  Span of the centrebody defined by the outer edge 
%               of the body flaps
%           - PhiCS: A (6*NumPanels)-by-NumCS matrix where 
%               PhiCS(6*i-5:6*i,k) contains [x,y,z,theta_x,theta_y,theta_z]
%               degrees of freedom for the i-the panel for 1 radian 
%               deflection mode of k-th control surface.
%       - cref: Mean aerodynamic chord
%       - geom: A structure with the following fields
%           - strcgrid: A field containing the following variables about 
%               the structural grid
%               - numnodes: Number of FEM nodes
%               - coords: A 2-by-numnodes matrix where the i-th column
%                   contains [x;y] coordinate of the i-th node
%                   The coordinate frame is such that:
%                   A) Origin at the vehicle nose, 
%                   B) +x pointing to rear or aircraft, and 
%                   C) +y pointing starboard (to the right when looking
%                        forward)
%                   D) +z pointing upward 
%
% Outputs:
%       - geom: A structure containing
%           - strcgrid: A field containing the following variables about 
%               the structural grid
%               - numnodes: Number of FEM nodes
%               - coords: A 2-by-numnodes matrix where the i-th column
%                   contains [x;y] coordinate of the i-th node
%                   The coordinate frame is such that:
%                   A) Origin at the vehicle nose, 
%                   B) +x pointing to rear or aircraft, and 
%                   C) +y pointing starboard (to the right when looking
%                       forward)
%                   D) +z pointing upward 
%           - c_ref: Mean aerodynamic chord
%           - splinegrid: A structure containing
%               - coords: A 2-by-NumSplineNodes matrix where i-th column 
%                   contains [x;y] coordinates of the spline nodes. The 
%                   spline nodes are defined using same coordinate frame as
%                   the FEM nodes
%               - connec: A 1-by-NumSplineNodes vector where 
%                   connections(i)=j implies that the i-th spline node is 
%                   derived From the j-th FEM node. For FEM nodes, 
%                   connection(i) = i
%           - Aeorgrid: An 3-by-NumPanels matrix where i-th column contains
%                   the [x;y;z] coordinates of the centroid of the i-th 
%                   panel.
%           - PhiCS: A (6*NumPanels)-by-NumCS matrix where 
%                   PhiCS(6*i-5:6*i,k) contains 
%                   [x,y,z,theta_x,theta_y,theta_z] degrees of freedom for 
%                   the i-the panel for 1 radian deflection mode of k-th
%                   control surface.
%           - S: A (6*NumPanels)-by-NumPanels matrix where where blocks of 
%                   six rows give the [Ax; Ay; Az; Mx; My; Mz] data.
%                   Ax is the panel area projected perpendicular to the 
%                   x-axis. Ay and Az are defined similarly.
%                   Mx is th moment of the area about the x-axis. My and Mz
%                   are defined similarly.
%           - D1,D2: NumAeroPanels-by-(6*NumAeroPanels) matrices which 
%                   multiplies by 6-dof displacement and velocities to give
%                   downwash on each panel
%           - Tfs: A NumSplineNodes-by-(3*NumFEMNodes) matrix which when
%                   multiplied by the FEM mode shape (containing
%                   [Have,Roll,Pitch] of the FEM nodes) gives the resulting 
%                   heave deflections of the spline nodes
%           - Tsa: A (2*NumAeroNodes)-by-(NumSplineNodes) matrix which when
%                   multiplied by the heaves of spline nodes outputs 
%                   [Heave,Twist] of the aerodynamic nodes
%           - Tfa: A (2*NumAeroNodes)-by-(3*NumFEMNodes) matrix which when 
%                   multiplied by the FEM mode shape (containing 
%                   [Have,Roll,Pitch] of the FEM nodes) outputs 
%                   [Heave,Twist] of the aerodynamic nodes

%% Setup
Aerogrid = griddata.Aerogrid;           % Aerodynamic grid
coords_fem = geom.strcgrid.coords;      % FEM coordinates
BodySpan = griddata.BodySpan;           % Span of the center body
displayfigs = true;                     % Option to toggle figures

%% Create spline grid
[coords_spline,connections] = spline_grid(coords_fem,BodySpan); 
% coords_spline is 3-by-NumSplineNodes matrix where i-th column contains
% [x;y;z] coordinates of i-th spline node.
% connections(i) = j implies that the i-th spline node is derived from the
% j-th FEM node. For FEM nodes, connection(i) = i;

%% Calculate transformation matrices
% Calculate transformation from structural grid to spline grid
% Tfs is a NumSplineNodes-by-(3*NumFEMNodes) matrix which when multiplied 
% by the FEM mode shape (containing [Have,Roll,Pitch] of the FEM nodes) 
% gives the resulting heave deflections of the spline nodes
Tfs = calc_Tfs(coords_fem,coords_spline,connections);

% Calculate transformation from spline to aero grid
% Tsa is a (2*NumAeroNodes)-by-(NumSplineNodes) matrix which when 
% multiplied by the heaves of spline nodes outputs [Heave,Twist] of the 
% aerodynamic nodes
Tsa = calc_Tsa(coords_spline,Aerogrid);

% Calculate transformation from structural to aero grid
% Tfa is a (2*NumAeroNodes)-by-(3*NumFEMNodes) matrix which when 
% multiplied by the FEM mode shape (containing [Have,Roll,Pitch] of the 
% FEM nodes) outputs [Heave,Twist] of the aerodynamic nodes
Tfa = Tsa*Tfs;

%% Plot Flex modes
% Obtain flexible modes for aerodynamic grid using transformation matrix
Phif_aero = Tfa*modes.Phif;

% Display flexible modes
% NOTE: The flexible modes are scaled for better visualization. This
% scaling does not affect the calculated modes
flexfigscale = 2;
if displayfigs
    for i = 1:size(Phif_aero,2)
        figure
        plot3(Aerogrid(1,:),Aerogrid(2,:),Aerogrid(3,:),'b*')
        hold on
        plot3(Aerogrid(1,:)+flexfigscale*Phif_aero(1:6:end,i)',...
                Aerogrid(2,:)+flexfigscale*Phif_aero(2:6:end,i)', ...
                Aerogrid(3,:)+flexfigscale*Phif_aero(3:6:end,i)','r*')
        xlabel('x'); ylabel('y'); zlabel('z'); grid on;
        legend('Undeflected nodes','Deflected nodes')
        axis equal
        view(110,60)
        title(['Aerodynic grid flexible mode: ', num2str(i)])
    end
end
%% Calculate differentiation matrix
% Differentiation matrices are Npanels-by-(6*Npanels) matrice which 
% multiplies by 6-dof modeshape and velocities to give downwash on each 
% aerodynamic panel
[D1,D2] = getDiff(griddata,cref);

%% Store geometric properties in geom (to be used in SIMULINK)
geom.c_ref = cref;
geom.splinegrid.coords = coords_spline;
geom.splinegrid.connec = connections;
geom.Aerogrid = Aerogrid;
geom.PhiCS_aero = griddata.PhiCS_aero;
geom.Phif_aero = Phif_aero;
geom.Phib_aero = griddata.Phib_aero;
geom.S = griddata.S;

geom.D1 = D1;
geom.D2 = D2;

geom.Tfa = Tfa;

save('Data_files/GAM_data.mat','geom');