function [acgeom,LE,TE,LwlLE,LwlTE,RwlLE,RwlTE,CSVertices] =  Geri_Outline
%Drawing_Geri   Outline of the aircraft planform
%   This script defines the outline of the aircraft planform and winglets
%
%   All the outputs are in SI units
%   - acgeom: A structure containing the following field
%       - cRoot: Root chord
%       - WingSpan: Wing span
%       - BodySpan: Span of the centrebody defined by the outer edge of the
%           body flaps
%       - WingletHalfHeight: Half the height of the winglet
%       - BodyFlapChordSize: Body flap size along the chord
%       - WingFlapChordSize = Wing flap size along the chord
%       - AllFlapYEdges: NumFlaps-by-2 matrix where the first column
%           contains y-coodinate of the inner (close to root) edge of 
%           the control surfaces. The rows are correspond to the
%           control surfaces in the order [L1,L2,L3,L4,R1,R2,R3,R4]
%   - LE: A 3-by-N matrix where each column provides [x;y;z] coordinates 
%       of N points on the leading edge of the aircraft planform. The
%       columns of points are ordered from left to right.
%   - TE: Same as LE for M points on trailing edge
%   - LwlLE: Same for leading edge of the left winglet
%   - LwlTE: Same for trailing edge of the left winglet
%   - RwlLE: Same for leading edge of the right winglet
%   - RwlTE: Same for trailing edge of the right winglet
%   - CSVertices: 3-by-4-by-NFlaps array where CSVertices(:,:,i) describes
%       the [x,y,z] coordinates of the 4 vertices of the ith control 
%       surface in the order [L1,L2,L3,L4,R1,R2,R3,R4]. The vertices are
%       ordered as [inner TE, inner hinge-line, outer hinge-line, outer TE]
% 
%   The points are defined in a coordinate frame with:
%       A) Origin at the vehicle nose, 
%       B) +x pointing to rear or aircraft, and 
%       C) +y pointing starboard (to the right when looking forward)
%       D) +z pointing upward 

%% Setup
% Option to toggle figure displays
displayfigs = true;

%% Geometric definitions
cRoot = 0.8522;     % Root chord, m
b = 3.055;          % Total span, m
h = 0.1397;         % Winglet half height, m

% Flap Sizes
% Note: The actual flaps are aligned perpendicular to the trailing
% edge.  These flap dimensions are approximate and aligned with the
% flow (+x) direction.  See more detailed comments below.
cbf = 0.1112;            % Body flap size along chord, m 
cwf = 0.0757;            % Wing flap size along chord, m 

% Left/right edges for flaps on right side of the aircraft
% FlapY(i,:) contains the (approximate) y-coordinates for the left
% and right side of flap i. 
FlapY = [0.0497 0.3405; ...
         0.3934 0.7437; ...
         0.7437 1.0941; ...
         1.0941 1.4444];

% Define span of body (fuselage) by the outer edge of body flap
BodySpan = 2*FlapY(1,2);     
     
% Store basic geometric data in a structure
acgeom.cRoot = cRoot;
acgeom.WingSpan = b;
acgeom.BodySpan = BodySpan;
acgeom.WingletHalfHeight = h;
acgeom.BodyFlapChordSize = cbf;
acgeom.WingFlapChordSize = cwf;
acgeom.AllFlapYEdges = [-fliplr(FlapY); FlapY];

%% Leading edge of aircraft planform
% LE is a 3-by-N matrix where each column provides [x;y;z] coordinates 
% (in m) of a point on the leading edge of the aircraft planform.
% The columns of points are ordered from left to right wingtip, i.e.
% in order of increasing values of y.
% Note: The wing sweep is approximately 22degs.

% Leading edge nose
p0 = [0; 0; 0];                    

% Leading edge right 
p1R = [0.331; 0.194; 0];           % First kink in right LE
p2R = [0.842; 1.464; 0];           % Second kink in right LE
p3R = [0.943; b/2; 0];            % Wingtip on right LE

% Leading edge left (symmetric about centerline)
p1L = p1R; p1L(2)=-p1L(2);
p2L = p2R; p2L(2)=-p2L(2);
p3L = p3R; p3L(2)=-p3L(2);

% Leading Edge: Points ordered from left wingtip to right wingtip
LE = [p3L p2L p1L p0 p1R p2R p3R];

%% Trailing edge of aircraft planform
% TE is a 3-by-M matrix where each column provides [x;y;z] coordinates 
% (in m) of a point on the trailing edge of the aircraft planform.
% The columns of points are ordered from left to right wingtip, i.e.
% in order of increasing values of y.

% Trailing edge nose
p0 = [cRoot; 0; 0];                    

% Trailing edge right 
p1R = [0.700; 0.372; 0];        % Kink in right TE
p2R = [1.164; b/2; 0];          % Wingtip on right TE

% Leading edge left (symmetric about centerline)
p1L = p1R; p1L(2)=-p1L(2);
p2L = p2R; p2L(2)=-p2L(2);

% Leading Edge: Points ordered from left wingtip to right wingtip
TE = [p2L p1L p0 p1R p2R];

%% Full aircraft planform (LE and TE)
Planform = [LE fliplr(TE)];

%% Right winglet
% The winglet is symmetric above and below the planform.
p1R  = [1.060; b/2; -h];        % LE of right winglet at bottom
p2R  = [LE(1,end); b/2; 0];     % LE of right winglet at planform  
p3R  = [1.060; b/2; h];         % LE of right winglet at top

p4R  = [TE(1,end); b/2; -h];    % TE of right winglet at bottom
p5R  = [TE(1,end); b/2; h];     % TE of right winglet at top

RwlLE = [p1R p2R p3R];
RwlTE =[p4R p5R];

%% Left winglet
LwlLE = RwlLE;
LwlLE(2,:) = -LwlLE(2,:);

LwlTE = RwlTE;
LwlTE(2,:) = -LwlTE(2,:);

%% Control surfaces
% RCS is a 3-by-4-by-4 array where RCS(:,:,i) describes the outline
% of control surface i on the right side of the aircraft. Each
% column of RCS(:,:,i) provides [x;y;z] coordinates (in m) of a point on 
% the control surface.  The control surfaces are numbered from right
% body flap (i=1) to outer most wing surface (i=4).
%
% LCS is defined similarly for the control surfaces on the left side.
% These control surfaces are also numbered from left body flap (i=1) to 
% outer most wing surface (i=4).
%
% Note: The left/right edges of the control surfaces are perpendicular to
% the trailing edge of the aircraft.  These edges do not align with the
% flow direction. A slightly modified version of the control surfaces are 
% defined below with edges aligned with the flow direction. (This is 
% required for the DLM/VLM code).  The y-coordinate of the left edge is
% defined as the average y-coordinate of the actual left edge.  The right
% edge is defined similarly.  The top/bottom edges are defined parallel
% to the trailing edge with. The dimensions are defined to maintain the
% same area and centroid as the original flap.

% Interpolate to find locations of flap trailing edge x-coordinates
FlapX = interp1( TE(2,:), TE(1,:), FlapY);

% Right and Left Control Surfaces
RCS = zeros(3,4,4);
LCS = RCS;
for i=1:4    
    if i==1
        % Body Flap
        d = cbf;
    else
        % Wing Flaps
        d = cwf;
    end
       
    p1 = [FlapX(i,1); FlapY(i,1); 0];
    p2 = [FlapX(i,1)-d; FlapY(i,1); 0];
    p3 = [FlapX(i,2)-d; FlapY(i,2); 0];
    p4 = [FlapX(i,2); FlapY(i,2); 0];
    
    RCS(:,:,i) = [p1 p2 p3 p4];
    LCS(:,:,i) = fliplr(RCS(:,:,i));
    LCS(2,:,i) = -LCS(2,:,i);
end
CSVertices = cat(3,LCS,RCS);

%% Display aircraft outline
if displayfigs
    figure
    plot3(Planform(1,:),Planform(2,:),Planform(3,:),'k');
    axis([-1.75 1.75 -1.75 1.75 -1.75 1.75]);
    hold on;
    plot3(LwlLE(1,:),LwlLE(2,:),LwlLE(3,:),'r');
    plot3(LwlTE(1,:),LwlTE(2,:),LwlTE(3,:),'r');
    plot3(RwlLE(1,:),RwlLE(2,:),RwlLE(3,:),'r');
    plot3(RwlTE(1,:),RwlTE(2,:),RwlTE(3,:),'r');
    axis equal
    
    for i=1:4
        plot3(RCS(1,:,i), RCS(2,:,i), RCS(3,:,i),'c');
        plot3(LCS(1,:,i), LCS(2,:,i), LCS(3,:,i),'c');
    end    
    hold off
end
