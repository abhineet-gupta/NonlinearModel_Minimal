%% mAEWing1_points_Analysis.m: Calculates dimensions of the aircraft
% This script uses the points recorded in the file 'mAEWing1_points.pdf'
% (obtained from the cad model of the aircraft) to calculate the dimensions
% of the aircraft. This is used for gridding the aircaft for DLM/VLM based
% aerodynamic model. Note that the coordinate system used in the PDF file
% is different from the aerodynamic model coordinate. But this script is
% only used to get the dimensions of the aircraft so the coordinate system
% is irrelevent.

%% Clear workspace and close figures
clear
close all
clc

%% Set-up
in2m = 0.0254;      % Inches to meters

%% Planform points
PFpoints = in2m* [-13.03, -7.63
            -15.78, -14.63
            -33.15, -57.63
            -37.14, -60.13      
            -45.84, -60.13
            -44.93, -57.63
            -44.42, -56.35
            -38.85, -42.56
            -33.27, -28.76
            -27.70, -14.97
            -27.56, -14.63
            -27.86, -14.20
            -32.73, -2.75
            -33.55, 0
            -29.03, -1.16
            -24.16, -12.61
            -25.14, -16.01
            -30.71, -29.80
            -36.28, -43.59
            -41.86, -57.38];

%% Winglet points
WLpoints = in2m* [-37.14, -60.13, -0.78
                -41.75, -60.13, -6.28
                -45.64, -60.13, -6.28
                -45.64, -60.13, 4.72
                -41.75, -60.13, 4.72];

%% Calculate planform dimensions
cRoot = norm(PFpoints(14,1))                    % Root chord
b = 2*norm(PFpoints(14,2)-PFpoints(4,2))        % Total wingspan
h = norm(WLpoints(3,:) - WLpoints(4,:))/2       % Winglet half height

%% Control surface lengths
% Body flap angle from -axis
BFangle = atan2d(PFpoints(12,2) - PFpoints(13,2), ...   
                PFpoints(12,1) - PFpoints(13,1));
% Body flap length along chord
cbf = norm(PFpoints(13,:)-PFpoints(15,:))/-sind(BFangle)


% Body flap angle from -axis
WFangle = atan2d(PFpoints(9,2) - PFpoints(10,2), ...   
                PFpoints(9,1) - PFpoints(10,1));
% Body flap length along chord
cwf = norm(PFpoints(17,:)-PFpoints(10,:))/-sind(WFangle)

%% Control sufrace location
Flapy = [];

% Body flap
Flapy = -[
      (PFpoints(13,2)+PFpoints(15,2))/2 (PFpoints(12,2)+PFpoints(16,2))/2;
      (PFpoints(10,2)+PFpoints(17,2))/2 (PFpoints(9,2)+PFpoints(18,2))/2;
      (PFpoints(9,2)+PFpoints(18,2))/2 (PFpoints(8,2)+PFpoints(19,2))/2;  
      (PFpoints(8,2)+PFpoints(19,2))/2 (PFpoints(7,2)+PFpoints(20,2))/2;
       ]

%% Leading edge points
pLE = [ 0 0 0;
        -PFpoints(1,:),0;
        -PFpoints(3,:),0;
        -PFpoints(4,:),0;]
    
%% Trailing edge points
pTE = [-PFpoints(14,:),0;
       -PFpoints(11,:),0;
       -PFpoints(5,:),0;]
   
%% Winglet leading edge
% Note that winglet is displace around 2 cm in the z direction (upwards)
% this shift is ignored and the winglets are assumed to be symmetric about
% the x-y plane. Therefore the origin is moved
shift = [0,0,WLpoints(1,3)];
wlLE = [-(WLpoints(5,:)-shift);
        -(WLpoints(1,:)-shift);
        -(WLpoints(2,:)-shift);]
wlTE = [-(WLpoints(4,:)-shift);
        -(WLpoints(3,:)-shift);]