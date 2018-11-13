% specify gmsh file
% read it in
% get nodes and elements
% (for now, just manually making up nodes and elements)
% This is problem 1 part d from homework 3
digits(8)
node1 = Node(1, 0, 1000, 0);
node2 = Node(2, 0, 0, 0);
node3 = Node(3, 1000, 1000, 0);
node4 = Node(4, 1000, 0, 0);
node5 = Node(5, 2000, 0, 0);

nodes = [node1, node2, node3, node4, node5];

elements(1:7) = Linear2DElement(node1, node2);
elements(1) = Linear2DElement(node1, node2);
elements(2) = Linear2DElement(node1, node3);
elements(3) = Linear2DElement(node1, node4);
elements(4) = Linear2DElement(node2, node4);
elements(5) = Linear2DElement(node3, node4);
elements(6) = Linear2DElement(node3, node5);
elements(7) = Linear2DElement(node4, node5);


kGlo = GlobalStiffnessMatrix(5, 2);
rGlo = GlobalLoadVector(5, 2);

for elm = elements
    localMat = elm.LocalStiffnessMatrix(68*500/elm.LengthSymbol);    
    kGlo.Add2NodeElementStiffnessMatrix(elm, localMat);    
    rGlo.Add2NodeElementLoad(elm, elm.CreateGravityLoadVector(-0.00981*(2.698e-6)*500*elm.LengthSymbol));
end
rGlo.Add2DConcentratedLoad(node5, 0, -50)
kGlo.K
rGlo.R
fixedNodes = [node1, node2];
kGlo.ApplyZeroBoundaryConditionToNode(fixedNodes)
rGlo.ApplyZeroBoundaryConditionToNode(fixedNodes)
kinv=inv(kGlo.K);
kinv*rGlo.R

  % nodes are a map of index to node data
  % elements is a list of elements
% this script knows what kind of problem it is, will call the approprate
% ...function on the particular elements to get the load vector (maybe a
% ...helper for particular types of problems)

% run through routine to get geometry portion of K
% run through routine to get geometry portion of load vector
% apply boundary conditions (to K AND Load Vector)
% invert remaining K an solve for displacements

