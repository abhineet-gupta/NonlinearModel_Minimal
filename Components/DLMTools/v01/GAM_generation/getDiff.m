function [D1,D2] = getDiff(griddata,c_ref)
%getDiff   Calculates differentiation matrices
%   Function to get differentiation matrices. The matrices are geometric 
%   transformations which when multiplied by the 6-dof deflection and 
%   velocity of the aerodynamic panel, output the downwash at the 
%   collocation point
% 
%   Inputs:
%       - PanelLengths: Npanels-by-1 vector containing the average 
%           x-dimension length of each panel
%       - c_ref: Reference length (mean aerodynamic chord)
% 
%   Outputs:
%       - D1: An Npanels-by-(6*Npanels) matrix which multiplies by 6-dof 
%           modeshape to give downwash on each panel
%       - D2: An Npanels-by-(6*Npanels) matrix which multiplies by 6-dof 
%           modal velocity to give downwash on each panel
%
%   Note that the differentiation matrices can be truncated later if less 
%   than 6-dofs are used

%% Number of panels
PanelLength = griddata.PanelLength;
Npanels = griddata.NumPanels;
nhatPF = griddata.nhatPF;
nhatWL = griddata.nhatWL;

NumPFPanels = sum(abs(nhatPF));
NumWLPanels = sum(abs(nhatWL));
idx_PF = 1:NumPFPanels;
idx_WL = NumPFPanels+1:Npanels;

%% Calculate the differentiation matrices
% See section 4-4.1 of the reference.
% Multiplies by 6 dof modeshape to give downwash
D1 = zeros(Npanels,6*Npanels); 

% Multplies by 6 dof modeshape velocity to give downwash
D2 = zeros(Npanels,6*Npanels); 

idx_PFpitch = idx_PF*6 - 1;  % Index corresponding to pitch angle
idx_PFheave = idx_PF*6 - 3;  % Index corresponding to heave displacement

idx_WLyaw = idx_WL*6;  % Index corresponding to pitch angle
idx_WLsway = idx_WL*6-4;  % Index corresponding to pitch angle

% D1 and D2 matrices calculated for all six degrees of freedom 
% Note: These matrices can be truncated later if only some DOFs are used.
D1(idx_PF,idx_PFpitch) = eye(NumPFPanels); 
D2(idx_PF,idx_PFheave) = -eye(NumPFPanels)*(2/c_ref);
D2(idx_PF,idx_PFpitch) = (2/c_ref)*diag(PanelLength(idx_PF))/4;

% NOTE: The D1,D2 matrices for winglets is different from the planform as 
D1(idx_WL,idx_WLyaw) = eye(NumWLPanels); 
D2(idx_WL,idx_WLsway) = eye(NumWLPanels)*(2/c_ref);
D2(idx_WL,idx_WLyaw) = (2/c_ref)*diag(PanelLength(idx_WL))/4;

% Negative sign since VLM/DLM accepts normalwash
D1 = -D1;  
D2 = -D2;