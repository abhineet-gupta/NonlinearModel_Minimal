function [xcoord_wingbody,ycoord_wingbody,xwl_u,zwl_u,xwl_l,zwl_l, ...
                                    Ycoord_sections] =  Geri_drawing()
%% BFFdrawing.m: Coordinates of the edges of the Geri aircraft
%
% Output:
% xcoord_wingbody: X-coords of LE and TE of body+wing
% ycoord_wingbody: Y-coords of LE or TE of body+wing
% xwl_u: x-coords of LE and TE of upper winglet
% zwl_u: z-coords of LE or TE of upper winglet
% xwl_l: x-coords of LE and TE of lower winglet
% zwl_l: z-coords of LE or TE of lower winglet
%
% The spacing is % 10^(-6). 
% All data is in SI

%% Geometric definitions
c_root = 33.55*0.0254;      % Root chord
b = 3.055;                  % Total span

% Ratio of distance of 1st kink in lead edge wrt half span
lead_break1 = 0.127;

% Ratio of distance of 2nd kink in lead edge wrt half span
lead_break2 = 0.958; 

% Ratio of distance of kink in the trailing edge wrt half span
trail_break = 0.243;  

theta_le_body = atan(13.03/7.63);       % Leading edge sweep angle
theta_te_body = atan(6/14.63);          % Trailing edge sweep angle
theta_wing = 22*pi/180;                 % Wing sweep angle
theta_wing_tip = atan(4/2.5);           % Wingtip sweep angle

xwl_LE_zend = 1.060;                    % xcoord of LE of winglet tip
wl_semiheight = 0.1397;                 % Winglet half height

Ycoord_sections = [0 0.0497 0.3405 0.3934 0.7437 1.0941 1.4444 b/2]; %

%% Centerbody coordinates (between root and kink1 in LE)
ycoord_wingbody = 0:0.00001:b/2;     % Divide halfspan  (note y is spanwise)
ind_ym = length(ycoord_wingbody);    % Number of spanwise coords
ind_LEkink1 = ceil(ind_ym*lead_break1);     % Index of the 1st kink in LE
ind_LEkink2 = ceil(ind_ym*lead_break2);     % Index of the 2nd kink in LE
ind_TEkink = ceil(ind_ym*trail_break);      % Index of the king in TE

% xcoord of the LE and TE of body
xbody_LE = -ycoord_wingbody(1:ind_LEkink1)*tan(theta_le_body);      
xbody_TE = ycoord_wingbody(1:ind_LEkink1)*tan(theta_te_body) - c_root; 

% xcoord of LE+TE of body
xbody_LE_TE = [xbody_LE;xbody_TE];  

%% wings coordinates (kink1 in LE to wingtip)
% xcord of LE between kink1 and 2
xwing_LE_kink1to2 = -(ycoord_wingbody(ind_LEkink1+1:ind_LEkink2) ... 
                    -ycoord_wingbody(ind_LEkink1))*tan(theta_wing) ...   % TODO you removed +1 from ind_LEkin1+1, check this
                    + xbody_LE(ind_LEkink1);

% xcord of LE between kink2 and wingtip
xwing_LE_kink2towingtip = -(ycoord_wingbody(ind_LEkink2+1:end) ...
                          -ycoord_wingbody(ind_LEkink2)) ...     % TODO  you removed +1 from ind_LEkink2+1, check this
                          *tan(theta_wing_tip) + ...
                          + xwing_LE_kink1to2(end);      

% xcord of TE between kink1 and kink in TE
xwing_TE_kink1tokinkTE = (ycoord_wingbody(ind_LEkink1+1:ind_TEkink) ... 
                         -ycoord_wingbody(ind_LEkink1)) ... 
                         *tan(theta_te_body) ...
                         + xbody_TE(ind_LEkink1);

% xcord of TE between kink in TE and wingtip
xwing_TE_kinkTEtowinglet = -(ycoord_wingbody(ind_TEkink+1:end) ...
                           -ycoord_wingbody(ind_TEkink)) ...
                           *tan(theta_wing) ...
                           +xwing_TE_kink1tokinkTE(end);

% xcoord of LE and TE of the wing
xwing_LE = [xwing_LE_kink1to2 xwing_LE_kink2towingtip];
xwing_TE = [xwing_TE_kink1tokinkTE xwing_TE_kinkTEtowinglet];

% xcoord of LE + TE of the wing
xwing_LE_TE = [xwing_LE;xwing_TE];

% xcoord of the LE and TE of full A/C (note the sign change)
xcoord_wingbody = -[xbody_LE_TE xwing_LE_TE];

% Note that at this point x = +ve from LE to TE and y is positive
% towards the right wing. Origin is nose

%% winglets definition (Note that upper and lower half has same x coord)

% Height of upper half of the winglet
zwl_u = 0:0.00001:wl_semiheight;
zwl_l = -zwl_u;

% xcoord of LE of winglet root
xwl_LE_z0 = -xwing_LE_kink2towingtip(end); 

% xcoord of TE of winglet root
xwl_TE_z0 = -xwing_TE_kinkTEtowinglet(end);

% xcoords of LE of winglet
xwl_LE = xwl_LE_z0 + zwl_u*(xwl_LE_zend - xwl_LE_z0)/max(zwl_u); 

% xcoords of TE of winglet
xwl_TE = ones(size(xwl_LE))*xwl_TE_z0; % winglet TE xcoord

xwl_u = [xwl_LE;xwl_TE];    % xcoord of upper half LE+TE
xwl_l = xwl_u;
