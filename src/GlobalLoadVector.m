classdef GlobalLoadVector < handle
    %GlobalLoadVector is a helper type to assist in creating a global load
    %vector.  
    
    properties (GetAccess = public, SetAccess = private)
        R
        % The number of variables per node (2 for 2D X and Y, 3 for 3D X, Y
        % and Z).
        NumberOfFreeElementsPerNode
    end
    
    methods
        function obj = GlobalLoadVector(numNodes, numFreeElements)
            obj.R = zeros(numNodes*numFreeElements, 1);
            obj.NumberOfFreeElementsPerNode=numFreeElements;
        end
        
        function Add2NodeElementLoad(obj, element, localLoad)
            node1 = element.Nodes(1);
            node2 = element.Nodes(2);
            obj.R((node1.Index-1)*2+1) = obj.R((node1.Index-1)*2+1) + localLoad(1);
            obj.R((node1.Index-1)*2+2) = obj.R((node1.Index-1)*2+2) + localLoad(2);
            
            obj.R((node2.Index-1)*2+1) = obj.R((node2.Index-1)*2+1) + localLoad(3);
            obj.R((node2.Index-1)*2+2) = obj.R((node2.Index-1)*2+2) + localLoad(4);
        end        
        
        function Add2DConcentratedLoad(obj, node, x, y)
            obj.R((node.Index-1)*2+1) = obj.R((node.Index-1)*2+1) + x;
            obj.R((node.Index-1)*2+2) = obj.R((node.Index-1)*2+2) + y;
        end
        
        function ApplyZeroBoundaryConditionToNode(obj, fixedNodes)
            [~, ind] = sort([fixedNodes.Index], 'descend');
            for node = fixedNodes(ind)
                obj.R((node.Index-1)*2+1) =[];
                obj.R((node.Index-1)*2+2) =[];
            end
        end
    end    
end


