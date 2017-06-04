classdef CurvedEdge
    %CURVEDEDGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodes
        xnodes
        ynodes
        order
        x_coord
        y_coord
    end
    
    properties (Dependent)
        sj
    end
    
    methods 
        function obj = CurvedEdge(nodes,xnodes,ynodes,order)
            obj.nodes = nodes;
            obj.xnodes = xnodes;
            obj.ynodes = ynodes;
            obj.order = order;
        
            obj.x_coord = zeros(1,obj.order);
            obj.y_coord = zeros(1,obj.order);
        end
        
        function obj = SetCircularEdge(obj,R,xc)
%
%           Compute each edge angles
%           ------------------------
            theta1 = atan2(obj.ynodes(1)-xc(2) , obj.xnodes(1)-xc(1));
            theta2 = atan2(obj.ynodes(2)-xc(2) , obj.xnodes(2)-xc(1));

            if ( theta2-theta1 > pi )
                theta2 = theta2 - 2*pi;
            elseif ( theta1 - theta2 > pi )
                theta1 = theta1 - 2*pi;
            end
                

%
%           Compute the curve angle
%           -----------------------
            theta = theta1 + 0.5*(obj.sj+1) * (theta2-theta1);
%
%           Compute the coordinates
%           -----------------------
            obj.x_coord = xc(1) + R * cos(theta);
            obj.y_coord = xc(2) + R * sin(theta);
            
        end
        
        function Plot(obj)
           plot(obj.x_coord,obj.y_coord,'-r'); 
        end
        
        function value = get.sj(obj)
            value = cos(pi-pi*double(0:obj.order)/double(obj.order));
        end
    end
    
end

