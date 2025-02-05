classdef Triangular3Node2DElement < Element
    %TRIANGULAR3NODE2DELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    % again, Z is ignored
    
    properties
        Node1
        Node2
        Node3
        PoisonsRatio
        E
        Thickness
    end
    
    methods
        function obj = Triangular3Node2DElement(node1, node2, node3, v, e, thickness)
            obj@Element([node1, node2, node3]);
            obj.Node1 = node1;
            obj.Node2 = node2;
            obj.Node3 = node3;
            obj.PoisonsRatio = v;
            obj.E = e;
            obj.Thickness = thickness;
        end
        
        
        function kLocal = LocalStiffnessMatrix3x3(obj, kxx, kyy)
            x1 = obj.Node1.X;
            x2 = obj.Node2.X;            
            x3 = obj.Node3.X;
            
            y1 = obj.Node1.Y;
            y2 = obj.Node2.Y;
            y3 = obj.Node3.Y;
            
            a1 = x2*y3-x3*y2;
            a2 = x3*y1-x1*y3;
            a3 = x1*y2-x2*y1;
            
            b1 = y2-y3;
            b2 = y3-y1;
            b3 = y1-y2;
            
            c1 = x3-x2;
            c2 = x1-x3;
            c3 = x2-x1;
            
            
            
            area = obj.AreaOfElement();
            kl11 = kxx*b1*b1+kyy*c1*c1;
            kl12 = kxx*b1*b2+kyy*c1*c2;
            kl13 = kxx*b1*b3+kyy*c1*c3;
            kl22 = kxx*b2*b2+kyy*c2*c2;
            kl23 = kxx*b2*b3+kyy*c2*c3;
            kl33 = kxx*b3*b3+kyy*c3*c3;
            
            kLocal = (1/(4*area))*[kl11, kl12, kl13;
                                   kl12, kl22, kl23;
                                   kl13, kl23, kl33];
            %det(kLocal)
            
        end
        
        function kLocal = LocalStiffnessMatrix(obj, coefficient)
            % From this link:  http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf
            x1 = obj.Node1.X;
            x2 = obj.Node2.X;            
            x3 = obj.Node3.X;
            
            y1 = obj.Node1.Y;
            y2 = obj.Node2.Y;
            y3 = obj.Node3.Y;
            
            x13 = x1-x3;
            x23 = x2-x3;
            x12 = x1-x2;
            x31 = -1*x13;
            x32 = -1*x23;
            x21 = -1*x12;
            
            y13 = y1-y3;
            y23 = y2-y3;
            y12 = y1-y2;
            y31 = -1*y13;
            y32 = -1*y23;
            y21 = -1*y12;
            
            
            detJ = x13*y23-y13*x23;
            area = (1/2)*abs(detJ);
            B = (1/detJ)*[y23, 0, y31, 0, y12, 0; 0, x32, 0, x13, 0, x21; x32, y23, x13, y31, x21, y12];
            v = obj.PoisonsRatio;
            D = obj.E/(1-v*v)*[1, v, 0; v, 1, 0; 0, 0, (1-v)/2];
            
            kLocal = obj.Thickness*area*transpose(B)*D*B;
        end
        
        function area = AreaOfElement(obj)
            area = 0.5*((obj.Node1.X-obj.Node3.X)*(obj.Node2.Y-obj.Node1.Y) - (obj.Node1.X-obj.Node2.X)*(obj.Node3.Y-obj.Node1.Y));
        end   
        
        function localPressureLoad = ComputeLocalPressureVector(obj, node1, node2, pressure)
            dist = Element.Compute2DLengthBetweenNodes(node1, node2);
            rVec = [-1, 0, -1, 0, 0, 0];
            if(node1.Index == obj.Nodes(1).Index && node2.Index == obj.Nodes(3).Index)
                rVec = [-1, 0, 0, 0, -1, 0];
            elseif(node1.Index == obj.Nodes(2).Index && node2.Index == obj.Nodes(3).Index)
                rVec = [0, 0, -1, 0, -1, 0];
            end
            localPressureLoad = rVec * pressure*dist/2;
            %TODO: rVec is not correct
        end
        
        function strain = ComputeStrain(obj, deltas)
            x1 = obj.Node1.X;
            x2 = obj.Node2.X;
            x3 = obj.Node3.X;
            
            y1 = obj.Node1.Y;
            y2 = obj.Node2.Y;
            y3 = obj.Node3.Y;
            bMat = [y2-y3, 0, y3-y1, 0, y1-y2, 0;
                    0, x3-x2, 0, x1-x3, 0, x2-x1;
                    x3-x2, y2-y3, x1-x3, y3-y1, x2-x1, y1-y2];
            A = obj.AreaOfElement();
            bMat = (1/2*A)*bMat;
            strain = bMat*transpose(deltas);
        end
        
        function stress = ComputeStress(obj, strain)
            v=obj.PoisonsRatio;
            eSc = obj.E/(1-v*v);
            eMat = eSc*[1, v, 0; 
                        v, 1, 0;
                        0, 0, (1-v)/2];
            stress = eMat*strain;
        end
    end    
end

