classdef (Abstract) Element
    % Element is the base class of all element types.  There will
    % eventually be many sub-classes. 
    
    properties (GetAccess = public, SetAccess = private)
        % A list of all the nodes in this element.  Note that the order 
        % does matter
        Nodes
    end
    
    methods
        function obj = Element(nodes)            
            obj.Nodes = nodes;
        end
        
        function trueOrFalse = ContainsNodes(obj, nodes)
            trueOrFalse = 1;
            for node = nodes
                localFound = 0;
                for myNode = obj.Nodes
                    if(myNode.Index == node.Index)
                        localFound = 1;
                        break;
                    end
                end
                if(localFound == 0)
                    trueOrFalse = 0;
                    break;
                end
            end
        end
    end

    methods (Static)
        function leng = ComputeLengthBetweenNodes(node1, node2)
            % Computes the mangitude between two nodes.
            diffX = node1.X-node2.X;
            diffY = node1.Y-node2.Y;
            diffZ = node1.Z-node2.Z;
            leng = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
        end
        
        function leng = Compute2DLengthBetweenNodes(node1, node2)
            % Computes the mangitude between just the X and Y components of 
            % two nodes.
            diffX = node1.X-node2.X;
            diffY = node1.Y-node2.Y;
            leng = sqrt(diffX*diffX+diffY*diffY);
        end
        
        
    end
    
    methods (Abstract)        
        LocalStiffnessMatrix(obj, coefficient)
    end
end

