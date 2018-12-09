classdef GlobalStiffnessMatrix < handle
    %GlobalStiffnessMatrix is a type that can be used to assemble the
    %global stiffness matrix from a set of Elements.  
    
    properties (GetAccess = public, SetAccess = private)
        % The stiffness matrix.  Note that this does mutate as elements are
        % added to the matrix.
        K
        % The number of variables per node (2 for 2D X and Y, 3 for 3D X, Y
        % and Z).
        NumberOfFreeElementsPerNode
    end
    
    methods
        function obj = GlobalStiffnessMatrix(numNodes, numFreeElements)
            % switch to the sparse matrix when we are happy with the
            % assembly part of FEM
            obj.K=sparse ( numNodes*numFreeElements , numNodes*numFreeElements );
            %obj.K = zeros(numNodes*numFreeElements, numNodes*numFreeElements);
            obj.NumberOfFreeElementsPerNode=numFreeElements;
        end
                
        function AddElementStiffnessMatrix(obj, element, kLocal)
            % Combine a local stiffness matrix into the global stiffness
            % matrix.
            n=obj.NumberOfFreeElementsPerNode;       
            for node1 = element.Nodes
                for node2 = element.Nodes
                    startRow = n*(find(element.Nodes==node1)-1);
                    startColumn = n*(find(element.Nodes==node2)-1);                    
                    for r = 1:n
                        for c = 1:n                    
                            kgRow = (node1.Index-1)*n+r;
                            kgCol = (node2.Index-1)*n+c;
                            loRow = startRow+r;
                            loCol = startColumn+c;
                            obj.K(kgRow, kgCol) = obj.K(kgRow, kgCol)+kLocal(loRow, loCol);
                        end
                    end
                end
            end
        end
        
        function ApplyZeroBoundaryConditionToNode(obj, fixedNodes)
            [~, ind] = sort([fixedNodes.Index], 'descend');
            for node = fixedNodes(ind)
                for el = 1:obj.NumberOfFreeElementsPerNode
                    obj.K((node.Index-1)*obj.NumberOfFreeElementsPerNode+1, :) =[];
                    obj.K(:, (node.Index-1)*obj.NumberOfFreeElementsPerNode+1) =[];
                end
            end
        end
        
        function ApplyZeroBoundaryConditionToIndices(obj, indices)
            [~, ind] = sort(indices, 'descend');
            for index = indices(ind)
                    obj.K(index,:) =[];
                    obj.K(:,index) =[];
            end
        end
    end    
end


