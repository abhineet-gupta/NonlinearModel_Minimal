function griddata = GridData_Geri_Quad(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate grid for the Geri and return information required by
% the getAIC function. The aircraft geometry is given by function
% BFFdrawing() .Note that BFFdrawing gives only half the wing.

% Input Variables: 
% NStr - no. of spanwise divisions on the wing (strips)
% N: no. of chordwise divisions of airframe exculding control surfaces
% C: no. of chordwise divisions of control surfaces
% Nz: no. of divisions along z axis for winglets

% Developed by:
% Abhineet Gupta
% Graduate assistant
% University of Minnesota
% Date: 11/22/2017

switch nargin
    case 1
        NStr = 70;
        N = 9; 
        C = 3;
        Nz = 15;
        Components = varargin{1};
    case 5
        NStr = varargin{1};
        N = varargin{2}; 
        C = varargin{3};
        Nz = varargin{4};
        Components = varargin{5};        
    otherwise
        error('Check number of arguments to GridData.m')
end

saveloc = Components.path_data;
filename = Components.Grid_file;

%% Define varaibles;
NodeData_temp = [];
PanelData_temp = [];
PanelData_temp_cs = [];
npwings = 0;
npwinglets = 0;


%% Obtain data of quadrilaterals on the Geri surface
[Quad_coord,Quad_points,Quad_edge_div,Quad_cs,Quad_cs_axis] = ...
                                                        Geri_drawing_Quad;

Num_Quads = size(Quad_points,1);    % Number of quadrilaterals

nodei = 1; % Initialize current node number
for i = 1:Num_Quads
    P1 = Quad_coord(Quad_points(i,1),:)';       % Identify 4 vertices
    P2 = Quad_coord(Quad_points(i,2),:)';
    P3 = Quad_coord(Quad_points(i,3),:)';
    P4 = Quad_coord(Quad_points(i,4),:)';
    
    len_e1 = norm(P2-P1);       % Identify 4 edge lengths
    len_e2 = norm(P3-P2);
    len_e3 = norm(P4-P3);
    len_e4 = norm(P1-P4);
    
    dc_e1 = (P2-P1)/len_e1';    % Identify direction cosines of the edges
    dc_e2 = (P3-P2)/len_e2';
    dc_e3 = (P4-P3)/len_e3';
    dc_e4 = (P1-P4)/len_e4';
    
    div1 = Quad_edge_div(i,1);  % Identify number of divisions required
    div2 = Quad_edge_div(i,2);
    
    
    % Identify coordinates of all the grid points on edge 1
    Pts_e1 = repmat(P1,1,div1+1) + dc_e1*(0:div1)*len_e1/div1;
    
    dc_panel_e2e4 = repmat(-dc_e4,1,div1+1) + (dc_e2 + dc_e4)/div1*(0:div1);
    dc_panel_e2e4 = repmat(dc_panel_e2e4,[1,div2+1]);
    
    len_panel_e2e4 = (len_e4 + (len_e2 - len_e4)/div1*(0:div1))'*(0:div2)/div2;
    len_panel_e2e4 = reshape(len_panel_e2e4,[1,(div1+1)*(div2+1)]);
    
    Pts_panel = repmat(Pts_e1,[1,div2+1]) + ...
                bsxfun(@times,dc_panel_e2e4,len_panel_e2e4);
    
    % Asign to NodeData
    NodeData_temp = [NodeData_temp;
                (nodei:(nodei-1+size(Pts_panel,2)))' Pts_panel'];
    
    % Asign PanelData
    for j = 1:div2
        for k = 1:div1
            PanelData_temp = [PanelData_temp; ...
                nodei + (j-1)*(div1+1) + k-1 , ...
                nodei + (j-1)*(div1+1) + k, ...
                nodei + j*(div1+1) + k-1,...
                nodei + j*(div1+1) + k];
            
            % PanelData_temp_cs(i) = j: ith panel is asigned to CS j
            PanelData_temp_cs = [PanelData_temp_cs;
                                 Quad_cs(i)];
        end
    end
    
    nodei = nodei+size(Pts_panel,2);
    
    if i<41
        npwings = npwings + (div1*div2);        % Num panels in wings
    else
        npwinglets = npwinglets + (div1*div2);  % Num panels in winglets
    end
end


%% Control surface modes.
Npanel = size(PanelData_temp,1);
phiCS = zeros(Npanel*6,max(Quad_cs));
CSind = zeros(Npanel,max(Quad_cs));

% Identify unique nodes and asign to NodeData
[NodeData,~,ni_new] = uniquetol(NodeData_temp(:,2:end),1e-6,'ByRows',true);
NodeData = [(1:size(NodeData,1))', NodeData];

for i = 1:size(PanelData_temp,1)
    % Asign new node index to PanelData_temp
    for j = 1:4
        k = NodeData_temp(:,1)==PanelData_temp(i,j);
        PanelData_temp(i,j) = ni_new(k);
    end
    
    % Find centroid of each quadrilateral
    aerogrid(:,i) = (NodeData(PanelData_temp(i,1),2:end) + ...
                     NodeData(PanelData_temp(i,2),2:end) + ...
                     NodeData(PanelData_temp(i,3),2:end) + ...
                     NodeData(PanelData_temp(i,4),2:end))'/4;
    
    % Find x length of each panel
    Plength(1,i) = (norm(NodeData(PanelData_temp(i,1),2) - ...
                   NodeData(PanelData_temp(i,4),2)) + ...
                   norm(NodeData(PanelData_temp(i,2),2) - ...
                   NodeData(PanelData_temp(i,3),2)))/2;
               
	% Vectors for each edge
    Pl1 = NodeData(PanelData_temp(i,2),2:end) - ...
                   NodeData(PanelData_temp(i,1),2:end);
    Pl2 = NodeData(PanelData_temp(i,4),2:end) - ...
                   NodeData(PanelData_temp(i,2),2:end);
    Pl3 = NodeData(PanelData_temp(i,3),2:end) - ...
                   NodeData(PanelData_temp(i,4),2:end);
    Pl4 = NodeData(PanelData_temp(i,1),2:end) - ...
                   NodeData(PanelData_temp(i,3),2:end);
               
    % Area of each panel
    PAreas(1,i) = norm(cross(Pl1,Pl2))/2 + norm(cross(Pl3,Pl4))/2;
    
    % Asigen control surface modes.
    cs_no =  PanelData_temp_cs(i);
    if cs_no >0
        % Coordinates of end points of the axis
        a = Quad_coord(Quad_cs_axis(cs_no,1),:);
        b = Quad_coord(Quad_cs_axis(cs_no,2),:);
        p = aerogrid(:,i)';
        dc_ab = (b-a)/norm(b-a);
        dist_p = norm(cross(p-a,p-b)/norm(b-a));
        phiCS((i*6 - 5):(i*6),cs_no) = ...
            [0;0;-dist_p; -dc_ab(1); dc_ab(2); 0];
        CSind(i,cs_no) = 1;
    end
    
    PSpan(i) = norm(Pl1(2:3));
end


    n_hat_w = [ones(npwings,1); zeros(npwinglets,1)];
    n_hat_wl = [zeros(npwings,1); -ones(npwinglets,1)];

PanelData = [(1:size(PanelData_temp,1))', PanelData_temp];

% Npanel = size(PanelData,1);
S = zeros(6*Npanel,Npanel);
rows1 = 6*[1:Npanel] - 3;
rows2 = 6*[1:Npanel] - 1;
S(rows1,:) = diag(PAreas(1,:));
S(rows2,:) = diag(PAreas.*Plength*0.25);

% xxx matrix containing aeras and momnt of areas about(quarter chord point)
% in rows, (3,9,15, ... ) a and (5,11,17,...) etc
S = sparse(S);



%% Plot Nodes
figure
hold on
for i=1:size(PanelData,1)
    P1 = NodeData(PanelData(i,2),2:end);
    P2 = NodeData(PanelData(i,3),2:end);
    P3 = NodeData(PanelData(i,4),2:end);
    P4 = NodeData(PanelData(i,5),2:end);
    
    L1 = [P1;P2]; L2 = [P2;P4]; L3 = [P4;P3]; L4 = [P3;P1]; 
    
    if PanelData_temp_cs(i) <1
        plot3(L1(:,1),L1(:,2),L1(:,3),'b');
        plot3(L2(:,1),L2(:,2),L2(:,3),'b');
        plot3(L3(:,1),L3(:,2),L3(:,3),'b');
        plot3(L4(:,1),L4(:,2),L4(:,3),'b');
    else
        plot3(L1(:,1),L1(:,2),L1(:,3),'r');
        plot3(L2(:,1),L2(:,2),L2(:,3),'r');
        plot3(L3(:,1),L3(:,2),L3(:,3),'r');
        plot3(L4(:,1),L4(:,2),L4(:,3),'r');
    end
end
axis equal


%% Save data
griddata.NodeData = NodeData;
griddata.PanelData = PanelData;
griddata.num_of_nodes = size(NodeData,1);
griddata.num_of_panels = size(PanelData,1);
griddata.aerogrid = aerogrid;
griddata.PanelLength = Plength;
griddata.PanelArea = PAreas;
griddata.S = S;
griddata.phiCS = phiCS;
griddata.PSpan = PSpan;
griddata.n_hat_w = n_hat_w;
griddata.n_hat_wl = n_hat_wl;
griddata.CSPanel = CSind;

save([saveloc,filename],'griddata'); 