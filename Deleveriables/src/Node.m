classdef Node < handle
    % A Node is simply 3 coordinates and an index. 
    % It is simply data, it does not do anything.
    
    properties (GetAccess = public, SetAccess = private)
        X
        Y
        Z
        Index
    end
    
    methods
        function obj = Node(i, x, y, z)
           obj.X = x;
           obj.Y = y;
           obj.Z = z;
           obj.Index = i;           
        end       
    end    
    
    methods (Static)
        function nodesInRange = FindNodesInXRange(nodes, minX, maxX)
            nodesInRange = [];
            for node = nodes
                if(node.X >= minX && node.X <= maxX)
                    nodesInRange = [nodesInRange, node];
                end
            end
        end
        
        
        function nodesInRange = FindNodesInYRange(nodes, minY, maxY)
            nodesInRange = [];
            for node = nodes
                if(node.Y >= minX && node.Y <= maxX)
                    nodesInRange = [nodesInRange, node];
                end
            end
        end
        
        function nodesInRange = FindNodesInRange(nodes, minX, maxX, minY, maxY)
            nodesInRange = [];
            for node = nodes
                if(node.Y >= minY && node.Y <= maxY && node.X >= minX && node.X <= maxX)
                    nodesInRange = [nodesInRange, node];
                end
            end
        end
    end
end

