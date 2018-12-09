classdef Linear2DElement < Element
    % Linear2DElement is the simplest type of element that I plan on making.
    % This can be used for truss style problems.  Note that the Z property
    % on the nodes will be ignored.
    
    properties (GetAccess = public, SetAccess = private)
        Node1
        Node2
        Length
    end
    
    methods
        function obj = Linear2DElement(node1, node2)            
            obj@Element([node1, node2]);
            obj.Node1 = node1;
            obj.Node2 = node2;
            obj.Length = obj.Compute2DLengthBetweenNodes(node1, node2);
        end
        
        function kLocal = LocalStiffnessMatrix(obj, coefficient)
            % Creates the local stiffness matrix.
            node1 = obj.Node1;
            node2 = obj.Node2;
            angle = atan2(node1.Y-node2.Y, node1.X - node2.X);
            c = cos(angle);
            s = sin(angle);
            % rounding to 1e-13, the less than impressive end of double precision in other languages
            c2 = round(c*c * 1e13)/1e13;
            s2=round(s*s* 1e13)/1e13;
            cs=round(c*s* 1e13)/1e13;
            realLength = obj.Length;            
            kLocal = (coefficient/realLength) * [c2, cs, -c2, -cs; cs, s2, -cs, -s2; -c2, -cs, c2, cs; -cs, -s2, cs, s2];            
        end
                
        function gravityLoadVector = CreateGravityLoadVector(obj, coef)
            node1 = obj.Node1;
            node2 = obj.Node2;
            realLength = obj.Compute2DLengthBetweenNodes(node1, node2);
            
            gravityLoadVector = coef *realLength* [0; 1 /2; 0; 1 /2];
        end
        
        function thermalLoadVector = CreateThermalLoadVector(obj, coef)
            node1 = obj.Node1;
            node2 = obj.Node2;
            angle = atan2(node1.Y-node2.Y, node1.X - none2.X);
            c = cos(angle);
            s = sin(angle);
            thermalLoadVector= coef*[-c; -s; c; s];
        end
        
        
        % in an effort to better understand shape functions I'm making this
        % set of functions to evaluate it all numerically.  This code is very 
        % inefficient, I know. I will improve it!
        function na = ShapeFunctionA(obj, s)
            l = obj.Length;
            na = 1-s/l;            
        end
        
        function nb = ShapeFunctionB(obj, s)
            l = obj.Length;
            nb = s/l;            
        end
        
        function val = func01ToIntegrate(obj, coef, s)
            l = obj.Length;
            dl = l/100.0;
            dnads = (obj.ShapeFunctionA(s+dl) - obj.ShapeFunctionA(s-dl))/(2*dl);
            dnbds = (obj.ShapeFunctionB(s+dl) - obj.ShapeFunctionB(s-dl))/(2*dl);
            val = coef*dnads*dnbds;
        end
        
        function val = func00ToIntegrate(obj, coef, s)
            l = obj.Length;
            dl = l/100.0;
            dnads = (obj.ShapeFunctionA(s+dl) -obj.ShapeFunctionA(s-dl))/(2*dl);            
            val = coef*dnads*dnads;
        end
        
        function val = func11ToIntegrate(obj, coef, s)
            l = obj.Length;
            dl = l/100;
            dnbds = (obj.ShapeFunctionB(s+dl) -obj.ShapeFunctionB(s-dl))/(2*dl);
            val = coef*dnbds*dnbds;
        end
        
        function altLocalK = LocalStiffnessMatrixNumerical(obj, coef)
            % NOTE that the coef should NOT have the 1/l in it as that
            % comes from the shape functions
            node1 = obj.Node1;
            node2 = obj.Node2;
            angle = atan2(node1.Y-node2.Y, node1.X - node2.X);
            c = round(cos(angle) * 1e13)/1e13;
            s = round(sin(angle) * 1e13)/1e13;
            tMat = [c, s, 0, 0; 0, 0, c, s];
            
            realLength = obj.Length;
            k00 = ode45(@(s, y) obj.func00ToIntegrate(coef, s), [0 realLength], 0);
            k01 = ode45(@(s, y) obj.func01ToIntegrate(coef, s), [0 realLength], 0);
            k11 = ode45(@(s, y) obj.func11ToIntegrate(coef, s), [0 realLength], 0);
            altLocalK = transpose(tMat)*[k00.y(end), k01.y(end); k01.y(end), k11.y(end)] * tMat;
        end
    end    
end

