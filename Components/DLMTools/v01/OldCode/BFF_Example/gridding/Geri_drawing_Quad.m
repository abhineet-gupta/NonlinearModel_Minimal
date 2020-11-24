function [Quad_coord,Quad_points,Quad_edge_div,Quad_cs,Quad_cs_axis] = Geri_drawing()
% Geri_drawing.m : Describes the Geri aircraft surface with a set of
% quadrilaterals.
%
% The quadrilaterals should cover every surface of the aicraft.
% These quadrilaterals would further be divided into finer grid.
%
% Outputs:
% Trap_coord

%% Constants
% Unit conversion constants
in2m = 0.0254;


%% Number of divisions

div_root_chord = 12;
div_cs_chord = 3;
div_body_span1 = 1;
div_body_span2 = 3;
div_body_span3 = 2;
div_body_span4 = 1;
div_wing_span1 = 1;
div_wing_span2 = 4;
div_wing_span3 = 4;
div_wing_span4 = 4;
div_wing_span5 = 2;
div_wl_height = 3;

%% geometric definitions
% Ponits on the aircraft used for creating quadrilateral panels
Quad_coord =  [0, 0, 0;
               13.03, 7.63, 0;
               15.78, 14.63, 0;
               33.15, 57.63, 0;
               37.14, 60.13, 0;
               45.84, 60.13, 0;
               44.42, 56.35, 0;  
               38.85, 42.56, 0;
               33.27, 28.76, 0;
               27.70, 14.97, 0;
               27.56, 14.63, 0;
               27.86, 14.20, 0;
               32.73, 2.75, 0;
               33.55, 0, 0;
               29.03, 1.16, 0;
               24.16, 12.61, 0;
               25.14, 16.01, 0;
               30.71, 29.80, 0;
               36.28, 43.59, 0;
               41.86, 57.38, 0;
               41.75, 60.13, -6.28;
               45.64, 60.13, -6.28;
               45.64, 60.13, 4.72;
               41.75, 60.13, 4.72;
              ];
Quad_coord = [Quad_coord;
               13.03, -7.63, 0;
               15.78, -14.63, 0;
               33.15, -57.63, 0;
               37.14, -60.13, 0;
               45.84, -60.13, 0;
               44.42, -56.35, 0;  
               38.85, -42.56, 0;
               33.27, -28.76, 0;
               27.70, -14.97, 0;
               27.56, -14.63, 0;
               27.86, -14.20, 0;
               32.73, -2.75, 0;
               29.03, -1.16, 0;
               24.16, -12.61, 0;
               25.14, -16.01, 0;
               30.71, -29.80, 0;
               36.28, -43.59, 0;
               41.86, -57.38, 0;
               41.75, -60.13, -6.28;
               45.64, -60.13, -6.28;
               45.64, -60.13, 4.72;
               41.75, -60.13, 4.72;
              ];
Quad_coord = Quad_coord* 0.0254; % Convert to meters

% Correct offset in z coordinate of winglets
winglet_z_offset = 0.78*0.0254; % Amount of offset in the solidworks model
Quad_coord([21:24 43:46],3) = Quad_coord([21:24 43:46],3) + ...
                             repmat(winglet_z_offset,8,1);

% Plot the points
figure
plot3(Quad_coord(:,1),Quad_coord(:,2),Quad_coord(:,3),'*');
xlabel('x'); ylabel('y'); zlabel('z');
labels = num2str((1:size(Quad_coord,1))','%d');    %'
text(Quad_coord(:,1), Quad_coord(:,2),Quad_coord(:,3),...
    labels, 'horizontal','left', 'vertical','bottom');
axis equal


%% Obtain more points required to divide the surface into quadrilaterals
% Find the intersection of line joining p1,p2 and line joining p3,p4

% Add point 47
p1 = Quad_coord(15,1:2);
p2 = Quad_coord(16,1:2);
p3 = Quad_coord(1,1:2);
p4 = Quad_coord(14,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 48
p1 = Quad_coord(1,1:2);
p2 = Quad_coord(2,1:2);
p3 = Quad_coord(15,1:2);
p4 = Quad_coord(15,1:2) +  [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];


% Add point 49
p1 = Quad_coord(2,1:2);
p2 = Quad_coord(2,1:2) + [1,0];
p3 = Quad_coord(15,1:2);
p4 = Quad_coord(16,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 50
p1 = Quad_coord(49,1:2);
p2 = Quad_coord(49,1:2) + [-cosd(23.24), -sind(23.24)];
p3 = Quad_coord(12,1:2);
p4 = Quad_coord(13,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 51
p1 = Quad_coord(2,1:2);
p2 = Quad_coord(3,1:2);
p3 = Quad_coord(16,1:2);
p4 = Quad_coord(16,1:2)+[1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 52
p1 = Quad_coord(11,1:2);
p2 = Quad_coord(3,1:2);
p3 = Quad_coord(17,1:2);
p4 = Quad_coord(18,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 53
p1 = Quad_coord(3,1:2);
p2 = Quad_coord(4,1:2);
p3 = Quad_coord(17,1:2);
p4 = Quad_coord(17,1:2)+[1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 54
p1 = Quad_coord(3,1:2);
p2 = Quad_coord(4,1:2);
p3 = Quad_coord(18,1:2);
p4 = Quad_coord(18,1:2) + [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 55
p1 = Quad_coord(3,1:2);
p2 = Quad_coord(4,1:2);
p3 = Quad_coord(19,1:2);
p4 = Quad_coord(19,1:2) + [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 56
p1 = Quad_coord(19,1:2);
p2 = Quad_coord(20,1:2);
p3 = Quad_coord(5,1:2);
p4 = Quad_coord(6,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 57
p1 = Quad_coord(23,:);
p2 = Quad_coord(24,:);
Quad_coord = [Quad_coord; p1+0.3*(p2-p1)];

% Add point 58
p1 = Quad_coord(22,:);
p2 = Quad_coord(21,:);
Quad_coord = [Quad_coord; p1+0.3*(p2-p1)];


% Add point 59
p1 = Quad_coord(1,1:2);
p2 = Quad_coord(25,1:2);
p3 = Quad_coord(37,1:2);
p4 = Quad_coord(37,1:2) +  [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];


% Add point 60
p1 = Quad_coord(25,1:2);
p2 = Quad_coord(25,1:2) + [1,0];
p3 = Quad_coord(37,1:2);
p4 = Quad_coord(38,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 61
p1 = Quad_coord(60,1:2);
p2 = Quad_coord(60,1:2) + [-cosd(23.24), sind(23.24)];
p3 = Quad_coord(35,1:2);
p4 = Quad_coord(36,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 62
p1 = Quad_coord(25,1:2);
p2 = Quad_coord(26,1:2);
p3 = Quad_coord(38,1:2);
p4 = Quad_coord(38,1:2)+[1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 63
p1 = Quad_coord(34,1:2);
p2 = Quad_coord(26,1:2);
p3 = Quad_coord(39,1:2);
p4 = Quad_coord(40,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 64
p1 = Quad_coord(26,1:2);
p2 = Quad_coord(27,1:2);
p3 = Quad_coord(39,1:2);
p4 = Quad_coord(39,1:2)+[1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 65
p1 = Quad_coord(26,1:2);
p2 = Quad_coord(27,1:2);
p3 = Quad_coord(40,1:2);
p4 = Quad_coord(40,1:2) + [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 66
p1 = Quad_coord(26,1:2);
p2 = Quad_coord(27,1:2);
p3 = Quad_coord(41,1:2);
p4 = Quad_coord(41,1:2) + [1,0];
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 67
p1 = Quad_coord(41,1:2);
p2 = Quad_coord(42,1:2);
p3 = Quad_coord(28,1:2);
p4 = Quad_coord(29,1:2);
[x_temp,y_temp] = findintersection(p1,p2,p3,p4);
Quad_coord = [Quad_coord; x_temp,y_temp,0];

% Add point 68
p1 = Quad_coord(44,:);
p2 = Quad_coord(43,:);
Quad_coord = [Quad_coord; p1+0.3*(p2-p1)];

% Add point 69
p1 = Quad_coord(45,:);
p2 = Quad_coord(46,:);
Quad_coord = [Quad_coord; p1+0.3*(p2-p1)];

num_pts_added = 23;

% Plot new points
hold on
plot3(Quad_coord(end-num_pts_added+1:end,1),...
    Quad_coord(end-num_pts_added+1:end,2),...
    Quad_coord(end-num_pts_added+1:end,3),'*');
labels = num2str((size(Quad_coord,1)- ...
    num_pts_added+1:(size(Quad_coord,1)))','%d');
text(Quad_coord(end-num_pts_added+1:end,1),...
    Quad_coord(end-num_pts_added+1:end,2),...
    Quad_coord(end-num_pts_added+1:end,3),...
    labels, 'horizontal','left', 'vertical','bottom');

%% Define quadrilaterals vertices
% Vertices of each quadrilateral
Quad_points = [28 27 42 67;
               27 66 41 42;
               66 65 40 41;
               65 64 39 40;
               64 26 63 39;
               26 62 38 63;
               62 25 60 38;
               25 59 37 60;
               59 1 47 37;
               1 48 15 47;
               48 2 49 15;
               2 51 16 49;
               51 3 52 16;
               3 53 17 52;
               53 54 18 17;
               54 55 19 18;
               55 4 20 19;
               4 5 56 20;
               67 42 30 29;
               42 41 31 30;
               41 40 32 31;
               40 39 33 32;
               39 63 34 33;
               63 38 35 34;
               38 60 61 35;
               60 37 36 61;
               37 47 14 36;
               47 15 13 14;
               15 49 50 13;
               49 16 12 50;
               16 52 11 12;
               52 17 10 11;
               17 18 9 10;
               18 19 8 9;
               19 20 7 8;
               20 56 6 7;
               58 56 6 22
               56 57 23 6;
               21 5 56 58;
               5 24 57 56;
               69 67 29 45;
               67 68 44 29;
               46 28 67 69;
               28 43 68 67;
               ];

%% Identify quadrilaterals assigned to control surfaces

% Quad_cs(i) == j: ith quadrilateral is part of jth control
% surface
Quad_cs = zeros(1,size(Quad_points,1)); 
Quad_cs(20:22) = 8:-1:6;
Quad_cs(25:26) = [5 5];
Quad_cs(29:30) = [1 1];
Quad_cs(33:35) = 2:4;

% Quad_cs_axis(i,:) = [j,k]: ith control panel has point j,k as it's axis
Quad_cs_axis = [15 16;
                17 18;
                18 19;
                19 20;
                38 37;
                40 39;
                41 40;
                42 41];

% Quad_edge_div(i,:) = [j,k] : ith qudrilateral should be divided into j
% panels along it's point(1-2) edge and in k panels along point(2-3) edge
Quad_edge_div = [div_wing_span5, div_root_chord;
                 div_wing_span4, div_root_chord;
                 div_wing_span3, div_root_chord;
                 div_wing_span2, div_root_chord;
                 div_wing_span1, div_root_chord;
                 div_body_span4, div_root_chord;
                 div_body_span3, div_root_chord;
                 div_body_span2, div_root_chord;
                 div_body_span1, div_root_chord;
                 div_body_span1, div_root_chord;
                 div_body_span2, div_root_chord;
                 div_body_span3, div_root_chord;
                 div_body_span4, div_root_chord;
                 div_wing_span1, div_root_chord;
                 div_wing_span2, div_root_chord;
                 div_wing_span3, div_root_chord;
                 div_wing_span4, div_root_chord;
                 div_wing_span5, div_root_chord;
                 div_wing_span5, div_cs_chord;
                 div_wing_span4, div_cs_chord;
                 div_wing_span3, div_cs_chord;
                 div_wing_span2, div_cs_chord;
                 div_wing_span1, div_cs_chord;
                 div_body_span4, div_cs_chord;
                 div_body_span3, div_cs_chord;
                 div_body_span2, div_cs_chord;
                 div_body_span1, div_cs_chord;
                 div_body_span1, div_cs_chord;
                 div_body_span2, div_cs_chord;
                 div_body_span3, div_cs_chord;
                 div_body_span4, div_cs_chord;
                 div_wing_span1, div_cs_chord;
                 div_wing_span2, div_cs_chord;
                 div_wing_span3, div_cs_chord;
                 div_wing_span4, div_cs_chord;
                 div_wing_span5, div_cs_chord;
                 div_wl_height, div_cs_chord;
                 div_wl_height, div_cs_chord;
                 div_wl_height, div_root_chord;
                 div_wl_height, div_root_chord;
                 div_wl_height, div_cs_chord;
                 div_wl_height, div_cs_chord;
                 div_wl_height, div_root_chord;
                 div_wl_height, div_root_chord;
    ];
end

function [x,y] = findintersection(p1,p2,p3,p4)
% [x,y] are the coordinates of intersection of lines (p1-p2) and (p3-p4)
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);
    x4 = p4(1); y4 = p4(2);
    A = [-(y2-y1) (x2-x1);
         -(y4-y3) (x4-x3)];
    B = [-x1*(y2-y1) + y1*(x2-x1);
         -x3*(y4-y3) + y3*(x4-x3)];
    pt = A\B;
    x = pt(1); y = pt(2);
end
       
