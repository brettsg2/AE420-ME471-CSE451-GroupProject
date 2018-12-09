function femSystem = ApplyStressToMesh(mesh, E, poison, pressure, pressureNodes, xBcs, yBcs) 
%TODO, pass in which nodes are gasket

%APPLYSTRESSTOMESH Summary of this function goes here

    N_n      = mesh.NumNodes;
    nodes       = mesh.Nodes;
    %----------------------------------------------
    % create the Node instances
    myElements(1:length(mesh.TwoDElements)) = Triangular3Node2DElement(nodes(1), nodes(2), nodes(3), poison, E, 1);

    % create the element instances
    for i = 1:length(mesh.TwoDElements)
        myElements(i) = Triangular3Node2DElement(mesh.Nodes(mesh.TwoDElements(i,1)), mesh.Nodes(mesh.TwoDElements(i,2)), mesh.Nodes(mesh.TwoDElements(i,3)), poison, E, 1);
    end

    % initalize the global items
    kGlo = GlobalStiffnessMatrix(N_n, 2);
    rGlo = GlobalLoadVector(N_n, 2);
    for i = 1:length(mesh.TwoDElements)
        localMat = myElements(i).LocalStiffnessMatrix(1);    
        kGlo.AddElementStiffnessMatrix(myElements(i), localMat);        
    end
    
    % apply pressure
%      for element = myElements
%         for pni1 = 1:length(pressureNodes)
%             pNode1 = pressureNodes(pni1);
%             if(element.ContainsNodes([pNode1]) == 1)
%                 for pni2 = (pni1+1):length(pressureNodes)
%                     pNode2 = pressureNodes(pni2);
% 
%                     if(element.ContainsNodes([pNode2]) == 1)
%                         localPressure = element.ComputeLocalPressureVector(pNode1, pNode2, pressure);
%                         rGlo.Add3NodeElementLoad(element, localPressure);
%                         break;
%                     end
%                 end
%             end
%         end
%     end

    for pni1 = 1:length(pressureNodes)
          pNode1 = pressureNodes(pni1);
          rGlo.Add2DConcentratedLoad(pNode1, 0,pressure);
    end

    % populate final answer
    finalAnswer = NaN([1,N_n*2]);
    % create copies of the entire stiffness matrix and load vector 
    % so post processing can be done later
    kBck = kGlo.K;
    rBck = rGlo.R;
   
    % deal with boundary conditions
    xIndices = zeros([1,length(xBcs) + length(yBcs)]);
    i=0;
    for fNode = xBcs
        i=i+1;
        finalAnswer((fNode.Index-1)*2+1) = 0;
        xIndices(i) = (fNode.Index-1)*2+1;
    end
    
    for fNode = yBcs
        i=i+1;
        finalAnswer((fNode.Index-1)*2+2) = 0;
        xIndices(i) = (fNode.Index-1)*2+2;
    end
    
    kGlo.ApplyZeroBoundaryConditionToIndices(xIndices)
    rGlo.ApplyZeroBoundaryConditionToIndices(xIndices)
    
    % solve!
    kinv=inv(kGlo.K);
    answer = kinv*rGlo.R;

    otherCounter = 1;
    for i = 1:N_n*2
        if(finalAnswer(i) ~= 0)
            finalAnswer(i) = answer(otherCounter);
            otherCounter=otherCounter+1;
        end
    end
    femSystem = FemSystem(kBck, rBck, finalAnswer, myElements);

end

