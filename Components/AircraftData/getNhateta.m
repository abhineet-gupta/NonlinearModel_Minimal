function Nhat_eta = getNhateta(Nhat0,Phif_aero,NumPanels,NumModes)
%% getNhateta  Obtains dependence of panel normal on flexible deflection
%
%   Input:
%       - Nhat0: A 3-by-NumPanels vector where i-th column contains the
%           direction cosinse of the panel normal vector
%       - Phif_aero: A (6*NumPanels-by-NumFlexModes) matrix where i-th
%           column contains the [Rx,Ry,Rz,Th_x,Th_y,Th-z] degrees of
%           freedome for each panel for ith-flexible mode
%   Output:
%       - Nhat_eta: A 3-by-NumPanels-by-NumFlexModes matrix where
%           Nhat_eta(:,:,i) contains the the depenency of the panel normals
%           on the flexible deflections eta:
%           Nhat = Nhat_0 + Nhat_eta(:,:,1)*eta1 + Nhat_eta(:,:,2)*eta2 ...

%% Panel normal due to flexible deflection

% Obtain rotation due to a unit flexible deflection (eta_f=1)
Phi_rotx = Phif_aero(4:6:end,:)';
Phi_roty = Phif_aero(5:6:end,:)';
Phi_rotz = Phif_aero(6:6:end,:)';

% Obtain relationship between panel normal and flexible deflection as a
% vector rotation of the of normal vection (theta x n0);
Nhat_eta = zeros(3,NumPanels,NumModes);
for i = 1:NumModes
    % Rotation of panels
    Phi_rot = [Phi_rotx(i,:);Phi_roty(i,:);Phi_rotz(i,:)]; 

    % Change in Normal vector
    Nhat_eta(:,:,i) = cross(Phi_rot,Nhat0,1);
end

end