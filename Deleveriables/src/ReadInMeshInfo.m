classdef ReadInMeshInfo < handle
    properties
        NodeIds
        NodeValues
        Nodes
        Elements
        TwoDNodes
        TwoDElements
        ElementTypes
        NumNodes
        NumElements
    end
    
    methods
        function obj = ReadInMeshInfo(filePathString, nodeScale, displacements)
            file    =  (filePathString);
            obj.NumNodes      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
            obj.NumElements      = dlmread(file,'',[7+obj.NumNodes 0 7+obj.NumNodes 0]);
            obj.NodeIds     = dlmread(file,'',[5 0 4+obj.NumNodes 0]);
            obj.NodeValues       = dlmread(file,'',[5 1 4+obj.NumNodes 3])*nodeScale;
            obj.Elements    = dlmread(file,'',[8+obj.NumNodes 0 7+obj.NumNodes+obj.NumElements 7]);
            obj.TwoDNodes = obj.NodeValues(:,1:2);
            obj.ElementTypes   = obj.Elements(:,2);        
                        
            %--- find the starting indices of 2D elements
            two_ind = 1;
            for i = 1:obj.NumElements
                if(obj.ElementTypes(i) ~= 2)
                    two_ind = two_ind+1;
                end
            end

            obj.TwoDElements(1:obj.NumElements-two_ind,1:3) = 0;
            k = 1;
            for i = two_ind:obj.NumElements
               obj.TwoDElements(k,1:3) = obj.Elements(i,6:8);
               k = k+1;
            end
            
            if(length(displacements) == obj.NumNodes*2)
                for n = 1:obj.NumNodes
                    obj.NodeValues(n,1) = obj.NodeValues(n,1)+displacements((n-1)*2+1);
                    obj.NodeValues(n,2) = obj.NodeValues(n,2)+displacements((n-1)*2+2);
                    
                    obj.TwoDNodes(n,1) = obj.TwoDNodes(n,1)+displacements((n-1)*2+1);
                    obj.TwoDNodes(n,2) = obj.TwoDNodes(n,2)+displacements((n-1)*2+2);
                end
            end
                        
            % create the Node instances
            myNodes(1:obj.NumNodes) = Node(0,0,0,0);
            nodeCount = 0;
            for tempNode = 1:obj.NumNodes
                nodeCount = nodeCount+1;
                myNodes(nodeCount) = Node(nodeCount, obj.NodeValues(nodeCount, 1), obj.NodeValues(nodeCount, 2), 0);
            end
            obj.Nodes = myNodes;
        end
    end    
end

