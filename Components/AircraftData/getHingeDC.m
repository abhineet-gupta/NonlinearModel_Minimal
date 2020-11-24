function CSData = getHingeDC(CSVertices,PanelIndexCS)
%% getHingeDC  Direction cosine of the hinge line of control surfaces
%
%   Obtains the direction cosines of the hinge line of the control surfaces
%
%   Input:
%       - CSVertices: 3-by-4-by-NFlaps array where CSVertices(:,:,i) 
%           describes the [x,y,z] coordinates of the 4 vertices of the ith
%           control surface in the order [L1,L2,L3,L4,R1,R2,R3,R4]. The 
%           vertices are ordered as [inner TE, inner hinge-line, outer 
%           hinge-line, outer TE]
%       - PanelIndexCS: A (xCSDiv+xPFDiv)-by-(NyPF-1) matrix. This is 
%           the dimension of the panels of the planform. If 
%           PanelIndexCS(i,j) == k then the i,j-th panel belongs to the 
%           control surface indexed k. The order of the control surfaces is
%           [L1,L2,L3,L4,R1,R2,R3,R4]
%
%   Output:
%       - CSData: A structure containg the following fields:
%           - CSHingeDC: A 3-by-NumFlaps matrix where ith column contains
%               the direction-cosines of the hinge line of th
%           - CSPanelIndex_all:
%           - DT_CS_Tau:

%%
NumCS = size(CSVertices,3);
CSHingeDC = zeros(3,NumCS);
CSPanelIndex_all = [];

for i = 1:NumCS
    HingePoint1 = CSVertices(:,2,i);    % Find ends of the hingeline
    HingePoint2 = CSVertices(:,3,i);
    CSHingeDC(:,i) = (HingePoint2-HingePoint1)/ ...
                        (norm(HingePoint2-HingePoint1));
    CSPanelIndex_all = [CSPanelIndex_all find(PanelIndexCS==i)'];
end

CSData.CSHingeDC = CSHingeDC;
CSData.CSPanelIndex_all = CSPanelIndex_all;
CSData.DT_CS_Tau = 300;

end