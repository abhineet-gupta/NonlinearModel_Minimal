function Dfinal = getVLM(P0,P1,P3,PAreas,n_hat_w,n_hat_wl)

% code developed using Katz & Plodkin as reference, to compute VLM solution 
% for a lifting surface. 
% returns the D matrix (downwash coeff) matrix, inverse of which gives the
% influence coeff matrix
% The AIC matrix as defined by Katz n plodkin provides vortex strength
% distribution across panels, further computation is carried out to get
% pressure differences as required
% multiply all y and z coordinates by by beta = sqrt(1-M^2) to account for
% compressibility



%% multiply y coordinates with beta
% beta = sqrt(1-(M^2));
% P0(:,2) = beta*P0(:,2);
% P1(:,2) = beta*P1(:,2);
% P3(:,2) = beta*P3(:,2);



%% get r1,r2,r0
epsilon = 10e-6;

% distance of vortex line beginning from control points
r1x= bsxfun(@minus,P0(:,1),P1(:,1)');  % r1x(i,j) = P0(i,1) - P1(j,1);
r1y= bsxfun(@minus,P0(:,2),P1(:,2)');
r1z= bsxfun(@minus,P0(:,3),P1(:,3)');

% distance of vortex line end from control points
r2x= bsxfun(@minus,P0(:,1),P3(:,1)');
r2y= bsxfun(@minus,P0(:,2),P3(:,2)');
r2z= bsxfun(@minus,P0(:,3),P3(:,3)');

% vortex line lengths
r0x = P3(:,1) - P1(:,1);
r0x = repmat(r0x,1,size(P0,1))';

r0y = P3(:,2) - P1(:,2);
r0y = repmat(r0y,1,size(P0,1))';

r0z = P3(:,3) - P1(:,3);
r0z = repmat(r0z,1,size(P0,1))';

% get normal vectors
n_hat_w = repmat(n_hat_w,1,size(P0,1));  % indicates cosine of panel dihedrals 
                                         % (1 for wing, 0 for winglet panels)

n_hat_wl = repmat(n_hat_wl,1,size(P0,1)); % sine of panel dihedrals 

%% induced velocity due to finite vortex line
r1Xr2_x = (r1y.*r2z) - (r2y.*r1z);
r1Xr2_y = -(r1x.*r2z) + (r2x.*r1z);
r1Xr2_z = (r1x.*r2y) - (r2x.*r1y);
mod_r1Xr2 = sqrt((r1Xr2_x.^2)+(r1Xr2_y.^2)+(r1Xr2_z.^2));

mod_r1 = sqrt((r1x.^2) + (r1y.^2) + (r1z.^2));
mod_r2 = sqrt((r2x.^2) + (r2y.^2) + (r2z.^2));

r0r1 = (r0x.*r1x) + (r0y.*r1y) + (r0z.*r1z);
r0r2 = (r0x.*r2x) + (r0y.*r2y) + (r0z.*r2z);

one = ones(size(r1Xr2_x));
D1_base = (1/(4*pi))*(one./(mod_r1Xr2.^2)).*((r0r1./mod_r1) - (r0r2./mod_r2));
D1_u = r1Xr2_x.*D1_base;
D1_v = r1Xr2_y.*D1_base;
D1_w = r1Xr2_z.*D1_base;

% adjust for non zero D matrix

ind =  find(mod_r1<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

ind =  find(mod_r2<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

ind =  find(mod_r1Xr2<epsilon);
D1_u(ind) = 0;
D1_v(ind) = 0;
D1_w(ind) = 0;

% get final D1 matrix 
% D1 matrix contains the perpendicular component of induced velocities at all panels.
% For wing panels, it's the z component of induced velocities (D1_w) while for
% winglets, it's the y component of induced velocities (D1_v)
D1 = (D1_w.*n_hat_w) + (D1_v.*n_hat_wl);

%% induced velocity due to inner semi-infinite vortex line

d2 = sqrt((r1y.^2) + (r1z.^2)); 
cosBB1 = 1;
cosBB2 = -r1x./mod_r1;
cosGamma = r1y./d2;
sinGamma = -r1z./d2;


D2_base = -(1/(4*pi))*(cosBB1 - cosBB2)./d2;
D2_u = zeros(size(D2_base));
D2_v = sinGamma.*D2_base;
D2_w = cosGamma.*D2_base;

ind = find(mod_r1<epsilon);
D2_u(ind) = 0;
D2_v(ind) = 0;
D2_w(ind) = 0;


% get final D2 matrix (same as D1)
D2 = (D2_w.*n_hat_w) + (D2_v.*n_hat_wl);

%% induced velocity due to outer semi-infinite vortex line
d3 = sqrt((r2y.^2) + (r2z.^2)); 
cosBT1 = r2x./mod_r2;
cosBT2 = -1;
cosGamma = -r2y./d3;
sinGamma = r2z./d3;

D3_base = -(1/(4*pi))*(cosBT1 - cosBT2)./d3;
D3_u = zeros(size(D3_base));
D3_v = sinGamma.*D3_base;
D3_w = cosGamma.*D3_base;

ind = find(mod_r2<epsilon);
D3_u(ind) = 0;
D3_v(ind) = 0;
D3_w(ind) = 0;

% get final D3 matrix (same as D1)
D3 = (D3_w.*n_hat_w) + (D3_v.*n_hat_wl);

%% total D
D = D1 + D2 + D3;

% Multiply additional factors to ensure that the D matrix maps pressure diff. 
% (rather than vortex strength) to downwash (ref Katz & Plodkin)  
deltaY = (P3(:,2) - P1(:,2)) + (P3(:,3) - P1(:,3));  % panel spans
F = 0.5*PAreas./deltaY;
F = repmat(F,1,size(P0,1))';
Dfinal = D.*F;   

return
