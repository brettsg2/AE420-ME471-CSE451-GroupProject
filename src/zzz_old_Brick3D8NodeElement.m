classdef Brick3D8NodeElement < Element
    %BRICK3D8NODEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Node1
        Node2
        Node3
        Node4
        Node5
        Node6
        Node7
        Node8
        LocalEtas
        LocalNus
        LocalZetas
        N
        E
        Nu
    end
    
    methods
        function obj = Brick3D8NodeElement(node1, node2, node3, node4, node5, node6, node7, node8, e, nu)
            obj@Element([node1, node2, node3, node4, node5, node6, node7, node8]);
            obj.Node1 = node1;
            obj.Node2 = node2;
            obj.Node3 = node3;
            obj.Node4 = node4;
            obj.Node5 = node5;
            obj.Node6 = node6;
            obj.Node7 = node7;
            obj.Node8 = node8;
            
            
%             obj.LocalNode1 = [-1, -1, -1];
%             obj.LocalNode2 = [1, -1, -1];
%             obj.LocalNode3 = [1, 1, -1];
%             obj.LocalNode4 = [-1, 1, -1];
%             obj.LocalNode5 = [-1, -1, 1];
%             obj.LocalNode6 = [1, -1, 1];
%             obj.LocalNode7 = [1, 1, 1];
%             obj.LocalNode8 = [-1, 1, 1];
            
            obj.LocalEtas = [-1, 1, 1, -1, -1, 1, 1, -1];
            obj.LocalNus = [-1, -1, 1, 1, -1, -1, 1, 1];
            obj.LocalZetas = [-1, -1, -1, -1, 1, 1, 1, 1];
            
            obj.N = 8;
            obj.E = e;
            obj.Nu=nu;
        end
        
        function n = ShapeFunctionForNode(obj, nodeIndex, eta, nu, zeta)
            n = (1/8)*(1-eta*obj.LocalEtas(nodeIndex))*(1-nu*obj.LocalNus(nodeIndex))*(1-zeta*obj.LocalZetas(nodeIndex));
        end
        
        function n = ShapeFunction(obj, xi, eta, zeta)
            n=0;
            for a = 1:8
                n = n + obj.ShapeFunctionForNode(xi, eta, zeta);
            end
        end                
        
        function x = LocalToX(obj, eta, nu, zeta, nodeIndex)
            x=0;
            for a = 1:8
                x = x + obj.ShapeFunctionForNode(eta, nodeIndex, nu, zeta);
            end
            theNode = obj.Nodes(nodeIndex);
            x = x*theNode.X;
        end
        
        function y = LocalToY(obj, eta, nu, zeta, nodeIndex)
            y=0;
            for a = 1:8
                y = y + obj.ShapeFunctionForNode(eta, nodeIndex, nu, zeta);
            end
            theNode = obj.Nodes(nodeIndex);
            y = y*theNode.Y;
        end
        
        function z = LocalToZ(obj, eta, nu, zeta, nodeIndex)
            z=0;
            for a = 1:8
                z = z + obj.ShapeFunctionForNode(eta, nodeIndex, nu, zeta);
            end
            theNode = obj.Nodes(nodeIndex);
            z = z*theNode.Z;
        end
        
        function dNdLoc = ComputeShapeFunctionDerivatives(obj, dEta, dNu, dZeta, eta, nu, zeta)
            dNdLoc = zeros(8, 3);
            for n = 1:obj.N
                dNdLoc(n, 1) = (obj.ShapeFunctionForNode(eta+dEta, nu, zeta, n)- obj.ShapeFunctionForNode(eta-dEta, nu, zeta, n))/dEta;
                dNdLoc(n, 2) = (obj.ShapeFunctionForNode(eta, nu + dNu, zeta, n)- obj.ShapeFunctionForNode(eta, nu - dNu, zeta, n))/dNu;
                dNdLoc(n, 3) = (obj.ShapeFunctionForNode(eta, nu, zeta+dZeta, n)- obj.ShapeFunctionForNode(eta, nu, zeta - dZeta, n))/dZeta;
            end
        end
        
        function J = ComputeJacobianOfShapeFunction(obj, dEta, dNu, dZeta, eta, nu, zeta)
            J = zeros(8,3);
            for n = 1:obj.N
                J(n, 1) = (obj.LocalToX(eta+dEta, nu, zeta, n)- obj.LocalToX(eta-dEta, nu, zeta, n))/dEta;
                J(n, 2) = (obj.LocalToY(eta, nu + dNu, zeta, n)- obj.LocalToY(eta, nu - dNu, zeta, n))/dNu;
                J(n, 3) = (obj.LocalToZ(eta, nu, zeta+dZeta, n)- obj.LocalToZ(eta, nu, zeta - dZeta, n))/dZeta;
            end
        end
        
        function kLocal = LocalStiffnessMatrixNumerically(obj, coefficientMatrix)
            A = coefficientMatrix;
            results = ode45(@(innerZeta, y) obj.FuncToIntegrateWrtZeta(A, innerZeta), [-1, 1], 0);
            kLocal = results.y(end);
        end
        
        function tbi = FuncToIntegrateWrtZeta(obj, A, zeta)
            results = ode45(@(innerEta, y) obj.FuncToIntegrateWrtNu(A, innerEta, zeta), [-1, 1], 0);
            tbi = results.y(end);
        end
        
        function tbi = FuncToIntegrateWrtNu(obj, A, nu, zeta)
            results = ode45(@(innerEta, y) obj.FuncToIntegrateWrtEta(A, innerEta, nu, zeta), [-1, 1], 0);
            tbi = results.y(end);
        end
        
        function tbi = FuncToIntegrateWrtEta(obj, A, eta, nu, zeta)
            J = obj.ComputeJacobianOfShapeFunction(obj, 0.05, 0.05, 0.05, eta, nu, zeta);
            dNdLoc = ComputeShapeFunctionDerivatives(obj, 0.05, 0.05, 0.05, eta, nu, zeta);
            %jInv = inv(J);
            B = dNdLoc/J;
            tbi = transpose(B)*A*B*det(J);            
        end
        
        function kLocal = LocalStiffnessMatrix(obj, coefficient)
%             cMag = obj.E./((1+obj.Nu)*(1-2*obj.nu))
%             cVec = [1-nu, nu, nu, 0, 0, 0;...
%                     nu, nu, 1-nu, 0, 0, 0;...
%                     0, 0, 0, (1-2*nu)/2, 0, 0;...
%                     0, 0, 0, 0, (1-2*nu)/2, 0;...
%                     0, 0, 0, 0, 0, (1-2*nu)/2];
%             ghostPoints = [-1/sqrt(3), 1/sqrt(3)];
%             coordinatesCenteredAtOrigin = zeros(8, 3);
            
            kLocal = LocalStiffnessMatrixNumerically(obj, coefficient);
                
        end
                
        function localLoad = LocalLoadVector(obj, Q)
            localLoad = obj.FinalLocalLoadIntegration(Q);
        end

        function integratedThird = FinalLocalLoadIntegration(obj, Q)
            results = ode45(@(innerZeta, y) obj.IntegrateLocalLoadTwice(Q, innerZeta), [-1, 1], 0);
            integratedThird = results.y(end);
        end
        
        function integratedTwice = IntegrateLocalLoadTwice(obj, Q, zeta)
            results = ode45(@(innerEta, y) obj.IntegrateLocalLoadOnce(Q, innerEta, zeta), [-1, 1], 0);
            integratedTwice = results.y(end);
        end
        
        function integratedOnce = IntegrateLocalLoadOnce(obj, Q, eta, zeta)
            results = ode45(@(innerXi, y) obj.LocalLoadFuncToIntegrate(Q, innerXi, eta, zeta), [-1, 1], 0);
            integratedOnce = results.y(end);
        end
        
        function innerLoadVecFunc = LocalLoadFuncToIntegrate(obj, Q, xi, eta, zeta)
            J = obj.ComputeJacobianOfShapeFunction(obj, 0.05, 0.05, 0.05, eta, nu, zeta);
            detJ = determinant(J);
            innerLoadVecFunc = Q*obj.ShapeFunction(xi, eta, zeta)*detJ;
        end
        
        
        
        
        function localLoadOnEdge = LocalLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta)
            localLoadOnEdge = obj.ThirdIntegrationOfLoadOnSide(q, fixedXi, fixedEta, fixedZeta);
        end
        
        function integrateLoadOnSideThird = ThirdIntegrationOfLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta)
            if(fixedEta ~= 0)
                zeta = fixedEta;
                integrateLoadOnSideThird = obj.SecondIntegrationOfLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta, xi, eta, zeta);
            else            
                results = ode45(@(innerZeta, y) obj.SecondIntegrationOfLoadOnSide(q, fixedXi, fixedEta, fixedZeta, xi, innerZeta), [-1, 1], 0);
                integrateLoadOnSideThird = results.y(end);
            end
        end
        
        function integrateLoadOnSideSecond = SecondIntegrationOfLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta, zeta)
            if(fixedEta ~= 0)
                eta = fixedEta;
                integrateLoadOnSideSecond = obj.FirstIntegrationOfLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta, xi, eta, zeta);
            else            
                results = ode45(@(innerEta, y) obj.FirstIntegrationOfLoadOnSide(q, fixedXi, fixedEta, fixedZeta, xi, innerEta, zeta), [-1, 1], 0);
                integrateLoadOnSideSecond = results.y(end);
            end
        end
        
        function integrateLoadOnSideFirst = FirstIntegrationOfLoadOnSide(obj, q, fixedXi, fixedEta, fixedZeta, eta, zeta)
            if(fixedXi ~= 0)
                xi = fixedXi;
                integrateLoadOnSideFirst = obj.InnerLoadOnSideToIntegrate(obj, q, fixedXi, fixedEta, fixedZeta, xi, eta, zeta);
            else            
                results = ode45(@(innerXi, y) obj.InnerLoadOnSideToIntegrate(q, fixedXi, fixedEta, fixedZeta, innerXi, eta, zeta), [-1, 1], 0);
                integrateLoadOnSideFirst = results.y(end);
            end
        end
        
        function innerLoadOnSideToIntegrate = InnerLoadOnSideToIntegrate(obj, q, fixedXi, fixedEta, fixedZeta, xi, eta, zeta)
            if(fixedXi ~= 0)
                xi = fixedXi;
            end
            if(fixedEta ~= 0)
                eta = fixedEta;
            end
            if(fixedZeta ~= 0)
                zeta = fixedZeta;
            end
            
            J = obj.ComputeJacobianOfShapeFunction(obj, 0.05, 0.05, 0.05, xi, eta, zeta);
            dNdLoc = ComputeShapeFunctionDerivatives(obj, 0.05, 0.05, 0.05, xi, eta, zeta);
            %jInv = inv(J);
            B = dNdLoc/J;
            magTerm = sqrt(B(1,1)*B(1,1)+B(1,2)*B(1,2)+B(1,3)*B(1,3));
            nVal = obj.ShapeFunction(xi, eta, zeta);
            innerLoadOnSideToIntegrate = q*nVal*magTerm;
        end
    end
    
end

