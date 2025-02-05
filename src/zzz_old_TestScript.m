% specify gmsh file
% read it in
% get nodes and elements
% (for now, just manually making up nodes and elements)

% This is problem 1 part d from homework 3
% first, makr our nodes
% DISTRICTRIZE
nodes(1:5) = Node(0,0,0,0);
node1 = Node(1, 0, 1000, 0);
node2 = Node(2, 0, 0, 0);
node3 = Node(3, 1000, 1000, 0);
node4 = Node(4, 1000, 0, 0);
node5 = Node(5, 2000, 0, 0);

nodes(1) = node1;
nodes(2) = node2;
nodes(3) = node3;
nodes(4) = node4;
nodes(5) = node5;

% then our elements
% ELEMENT PROPERTIES
elements(1:7) = Linear2DElement(node1, node2);
elements(1) = Linear2DElement(node1, node2);
elements(2) = Linear2DElement(node1, node3);
elements(3) = Linear2DElement(node1, node4);
elements(4) = Linear2DElement(node2, node4);
elements(5) = Linear2DElement(node3, node4);
elements(6) = Linear2DElement(node3, node5);
elements(7) = Linear2DElement(node4, node5);

% initialize our global items
kGlo = GlobalStiffnessMatrix(5, 2);
kGloAlt = GlobalStiffnessMatrix(5, 2);
rGlo = GlobalLoadVector(5, 2);
elements(1).LocalStiffnessMatrixNumerical(68*500)
elements(1).LocalStiffnessMatrix(68*500)

% fill in our global items
% ASSEMBLE
for elm = elements
    localMat = elm.LocalStiffnessMatrix(68*500);    
    kGlo.AddElementStiffnessMatrix(elm, localMat);        
    rGlo.Add2NodeElementLoad(elm, elm.CreateGravityLoadVector(-0.00981*(2.698e-6)*500));
end
%
% add in the concentrated load
rGlo.Add2DConcentratedLoad(node5, 0, -50)

% print out our values
%kGlo.K
%rGlo.R;

% apply boundary conditions
% BOUNDARY CONDITIONS
fixedNodes = [];
for node =nodes
    if(node.X == 0) 
        % I know this is a bad way of doing an expanding array, but for the
        % 2 nodes it is fine
        fixedNodes = [fixedNodes, node];
    end
end
kGlo.ApplyZeroBoundaryConditionToNode(fixedNodes)
rGlo.ApplyZeroBoundaryConditionToNode(fixedNodes)

% solve the equation
% SOLVE
kinv=inv(kGlo.K);
kinv*rGlo.R

% Still need to POST PROCESS

  % nodes are a map of index to node data
  % elements is a list of elements
% this script knows what kind of problem it is, will call the approprate
% ...function on the particular elements to get the load vector (maybe a
% ...helper for particular types of problems)

% run through routine to get geometry portion of K
% run through routine to get geometry portion of load vector
% apply boundary conditions (to K AND Load Vector)
% invert remaining K an solve for displacements

