classdef GlobalLoadVector < handle
    %GlobalLoadVector is a helper type to assist in creating a global load
    %vector.  
    
    properties %(GetAccess = public, SetAccess = private)
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
            n = obj.NumberOfFreeElementsPerNode;
            obj.R((node1.Index-1)*n+1) = obj.R((node1.Index-1)*n+1) + localLoad(1);
            obj.R((node1.Index-1)*n+2) = obj.R((node1.Index-1)*n+2) + localLoad(2);
            
            obj.R((node2.Index-1)*n+1) = obj.R((node2.Index-1)*n+1) + localLoad(3);
            obj.R((node2.Index-1)*n+2) = obj.R((node2.Index-1)*n+2) + localLoad(4);
        end  
        
        function Add3NodeElementLoad(obj, element, localLoad)
            node1 = element.Nodes(1);
            node2 = element.Nodes(2);
            node3 = element.Nodes(3);
            n = obj.NumberOfFreeElementsPerNode;
            obj.R((node1.Index-1)*n+1) = obj.R((node1.Index-1)*n+1) + localLoad(1);
            obj.R((node1.Index-1)*n+2) = obj.R((node1.Index-1)*n+2) + localLoad(2);
            
            obj.R((node2.Index-1)*n+1) = obj.R((node2.Index-1)*n+1) + localLoad(3);
            obj.R((node2.Index-1)*n+2) = obj.R((node2.Index-1)*n+2) + localLoad(4);
            
            obj.R((node3.Index-1)*n+1) = obj.R((node3.Index-1)*n+1) + localLoad(5);
            obj.R((node3.Index-1)*n+2) = obj.R((node3.Index-1)*n+2) + localLoad(6);
        end  
        
        function Add3NodeElementScalarLoad(obj, element, localLoad)
            node1 = element.Nodes(1);
            node2 = element.Nodes(2);
            node3 = element.Nodes(3);
            obj.R(node1.Index) = obj.R(node1.Index) + localLoad(1);
            obj.R(node2.Index) = obj.R(node2.Index) + localLoad(2);
            obj.R(node3.Index) = obj.R(node3.Index) + localLoad(3);
        end   
        
        function Add2DConcentratedLoad(obj, node, x, y)
            n = obj.NumberOfFreeElementsPerNode;
            obj.R((node.Index-1)*n+1) = obj.R((node.Index-1)*n+1) + x;
            obj.R((node.Index-1)*n+2) = obj.R((node.Index-1)*n+2) + y;
        end
        
        function ApplyZeroBoundaryConditionToNode(obj, fixedNodes)
            [~, ind] = sort([fixedNodes.Index], 'descend');
            n = obj.NumberOfFreeElementsPerNode;
            for node = fixedNodes(ind)
                obj.R((node.Index-1)*n+2) =[];                
                obj.R((node.Index-1)*n+1) =[];                
            end
        end
        
        function ApplyZeroBoundaryConditionTo1DNode(obj, fixedNodes)
            [~, ind] = sort([fixedNodes.Index], 'descend');
            n = obj.NumberOfFreeElementsPerNode;
            for node = fixedNodes(ind)
                obj.R(node.Index) =[];
            end
        end
        
        function ApplyZeroBoundaryConditionToIndices(obj, indices)
            [~, ind] = sort(indices, 'descend');
            for index = indices(ind)
                obj.R(index) =[];
            end
        end
    end    
end


