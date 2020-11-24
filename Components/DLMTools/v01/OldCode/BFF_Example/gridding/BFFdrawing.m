function [x_m,y_m,x_wingletu,z_wingletu,x_wingletl,z_wingletl] = BFFdrawing()
%% Draw BFF
% Returns coordinates of the BFF/miniMUTT aircraft geometry. The y coordinates are 
% evenly spaced between -b/2 to b/2 for the wing where b is the wingspan, the spacing is
% 10^(-6). x coordinates are a function of the corresponding to their y
% coordinates. z coordinates of the winglets are also spaced evenly spaced
% at 10^(-6).

%% geometric definitions
c_root = 33.55*0.0254;
b = 3.055;

% Ratio of distance of 1st kink in lead edge wrt half span
lead_break1 = 0.127;

% Ratio of distance of 2nd kink in lead edge wrt half span
lead_break2 = 0.958; 

% Ratio of distance of kink in the trailing edge wrt half span
trail_break = 0.243;  

theta_le_body = atan(13.03/7.63);
theta_te_body = atan(6/14.63);
theta_wing = 22*pi/180;
theta_wing_tip = atan(4/2.5);

%% airframe coordinates generation (center body)
y_m = 0:0.00001:b/2;   % Divide halfspan  (note y is spanwise)
l = length(y_m);            % Index of all the spanwise strips
l1 = ceil(l*lead_break1);   % Index of the 1st kink in LE
l2 = ceil(l*lead_break2);   % Index of the 2nd kink in LE
l3 = ceil(l*trail_break);   % Index of the king in TE

xbody_top = -y_m(1:l1)*tan(theta_le_body);      % xcoord of the LE of body
xbody_bot = y_m(1:l1)*tan(theta_te_body) - c_root; % xcoord of TE of body
xbody = [xbody_top;xbody_bot];  % xcoord of LE+TE of body


%% wings coordinates
% xcord of LE between kink1 and 2
xwing_top_1 = -(y_m(l1+1:l2)-y_m(l1+1))*tan(theta_wing) + xbody_top(l1);

% xcord of LE between kink2 and wingtip
xwing_top_2 = -(y_m(l2+1:end)-y_m(l2+1))*tan(theta_wing_tip) + ...
                xwing_top_1(end);

% xcord of TE between kink1 and kink in TE
xwing_bot_1 = (y_m(l1+1:l3)-y_m(l1+1))*tan(theta_te_body) + xbody_bot(l1);

% xcord of TE between kink in TE and wingtip
xwing_bot_2 = -(y_m(l3+1:end)-y_m(l3+1))*tan(theta_wing) + ...
                xwing_bot_1(end);

% xcoord of LE and TE of the wing
xwing_top = [xwing_top_1 xwing_top_2];
xwing_bot = [xwing_bot_1 xwing_bot_2];
xwing = [xwing_top;xwing_bot];

% xcoord of the LE and TE of full A/C (note the sign change)
x_m = -[xbody xwing];

% Note that at this point x = +ve from LE to TE and y is positive
% towards the right wing. origin is nose


%% winglets definition
% Upper half
z_m = 0:0.00001:0.14;   % Height of upper half of the winglet
xw_lead_z0 = -xwing_top_2(end); % xcoord of LE of winglet root
xw_trail_z0 = -xwing_bot_2(end); % xcoord of TE of winglet root
xw_lead_zend = 1.056;   % xcoord of LE of winglet tip
xw_lead = z_m*(xw_lead_zend - xw_lead_z0)/max(z_m); 
xw_lead = xw_lead + xw_lead_z0; % winglet LE xcoord

xw_trail = ones(size(xw_lead))*xw_trail_z0; % winglet TE xcoord
x_wingletu = [xw_lead;xw_trail];    % xcoord of upper half LE+TE
z_wingletu = z_m;                   % zcoord of upper half LE+TE

%lower half (Same as upper half)
z_ml = 0:0.00001:0.14; 
xw_lead_z0l = -xwing_top_2(end);
xw_trail_z0l = -xwing_bot_2(end);
% xw_lead_zend = 1.049;
xw_leadl = z_ml*(xw_lead_zend - xw_lead_z0l)/max(z_ml);
xw_leadl = xw_leadl + xw_lead_z0l;
xw_traill = ones(size(xw_leadl))*xw_trail_z0l;
x_wingletl = [xw_leadl;xw_traill];
z_wingletl = z_ml;