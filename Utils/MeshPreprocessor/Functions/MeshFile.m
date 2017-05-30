classdef MeshFile
    %MESHFILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        no_of_elements
        no_of_nodes
        no_of_curvedEdges
        no_of_bdryEdges
        points_of_elements
        xp
        yp
        points_of_bdryEdges
        bdrymarker_of_edges
        curvesOrder
        curvedEdges
        p_refinementZoneConditions
    end
    
    properties (Dependent)
        sj
        points_of_curvedEdges
        x_curvilinear_coords
        y_curvilinear_coords
        pRefinement_zones        
        no_of_pRefinementzones
    end
    
    methods
        function self = MeshFile
           
            self.no_of_elements = 0;
            self.no_of_nodes = 0;
            self.no_of_curvedEdges = 0;
            self.no_of_bdryEdges = 0;
            self.curvesOrder = 0;            
            self.points_of_elements = [];
            self.xp = [];
            self.yp = [];
            self.points_of_bdryEdges = [];
            self.bdrymarker_of_edges = [];
            self.p_refinementZoneConditions = {};
            self.curvedEdges = {};
            
        end
        
        function self = Load(self,points_of_elements , x , y , points_of_bdryEdges , bdrymarker_of_edges , curvesOrder )
%
%           Load points of elements
%           -----------------------
            if ( size(points_of_elements,1) == 4 ) 
                self.points_of_elements = points_of_elements';
            elseif ( size(points_of_elements,2) == 4 ) 
                self.points_of_elements = points_of_elements;
            end
%
%           Compute the number of elements
%           ------------------------------
            self.no_of_elements = size(self.points_of_elements,1);
%           
%           Compute the number of nodes
%           ---------------------------
            self.no_of_nodes = length(x);
%
%           Load the nodes coordinates
%           --------------------------
            self.xp = zeros(self.no_of_nodes,1);
            self.xp = x;            
            self.yp = zeros(self.no_of_nodes,1);            
            self.yp = y;
%
%           Load the points of boundary edges
%           ---------------------------------
            if ( size(points_of_bdryEdges,1) == 2 )
                self.points_of_bdryEdges = points_of_bdryEdges';
            elseif ( size(points_of_bdryEdges,2) == 2 ) 
                self.points_of_bdryEdges = points_of_bdryEdges;
            end
%
%           Compute the number of boundary edges
%           ------------------------------------
            self.no_of_bdryEdges = size(self.points_of_bdryEdges,1);
%
%           Load the boundary markers
%           -------------------------
            self.bdrymarker_of_edges = zeros(1,self.no_of_bdryEdges);
            self.bdrymarker_of_edges = bdrymarker_of_edges;
%
%           Load the boundary curves order
%           ------------------------------
            self.curvesOrder = curvesOrder;
        end
%
%/////////////////////////////////////////////////////////////////////////////////////
%
%               NON CONFORMING GEOMETRIES PROCEDURES
%               ------------------------------------
%/////////////////////////////////////////////////////////////////////////////////////
%
        function self = SubdivideElement(self,eID)
%
%           Gather old nodes coordinates
%           ----------------------------
            x1 = [self.xp(self.points_of_elements(eID,1)) , self.yp(self.points_of_elements(eID,1)) ];
            x2 = [self.xp(self.points_of_elements(eID,2)) , self.yp(self.points_of_elements(eID,2)) ];
            x3 = [self.xp(self.points_of_elements(eID,3)) , self.yp(self.points_of_elements(eID,3)) ];
            x4 = [self.xp(self.points_of_elements(eID,4)) , self.yp(self.points_of_elements(eID,4)) ];
%
%           Get the new nodes coordinates
%           -----------------------------
            xn0 = (x1+x2+x3+x4)/4;
            xn{1} = (x1+x2)/2;
            xn{2} = (x2+x3)/2;
            xn{3} = (x3+x4)/2;
            xn{4} = (x4+x1)/2;
%
%           Get whether the new nodes already exist or not
%           ----------------------------------------------
            ID1 = self.points_of_elements(eID,1);
            ID2 = self.points_of_elements(eID,2);
            ID3 = self.points_of_elements(eID,3);
            ID4 = self.points_of_elements(eID,4);
            
            self.no_of_nodes = self.no_of_nodes + 1 ;
            IDn0 = self.no_of_nodes ;
            self.xp(IDn0) = xn0(1);
            self.yp(IDn0) = xn0(2);
            
            for i = 1 : 4
                create = true;
                for j = 1 : self.no_of_nodes
                    x = [self.xp(j),self.yp(j)];
                    
                    if ( self.almostEqual(x,xn{i}) )
                        IDn{i} = j;
                        create = false;
                        break;
                    end
                end
                
                if (create)
                    self.no_of_nodes = self.no_of_nodes + 1 ;
                    IDn{i} = self.no_of_nodes;
                    self.xp(IDn{i}) = xn{i}(1);
                    self.yp(IDn{i}) = xn{i}(2);
                end
            end
%
%           Add the elements
%           ----------------
            self.points_of_elements(eID,:) = [ID1,IDn{1},IDn0,IDn{4}];
            Nel = self.no_of_elements;
            self.points_of_elements(Nel+1,:) = [IDn{1},ID2,IDn{2},IDn0];
            self.points_of_elements(Nel+2,:) = [IDn0,IDn{2},ID3,IDn{3}];
            self.points_of_elements(Nel+3,:) = [IDn{4},IDn0,IDn{3},ID4];
            self.no_of_elements = self.no_of_elements + 3 ;
%
%           Search for points of edges and replace
%           --------------------------------------
            for i = 1 : self.no_of_bdryEdges
                if ( self.sameEdges ( self.points_of_bdryEdges(i,:) , [ID1,ID2] ) ) 
                    self.points_of_bdryEdges(i,:) = [ID1,IDn{1}];
                    self.no_of_bdryEdges = self.no_of_bdryEdges + 1 ;
                    self.points_of_bdryEdges(self.no_of_bdryEdges,:) = [IDn{1},ID2];
                    self.bdrymarker_of_edges(self.no_of_bdryEdges) = self.bdrymarker_of_edges(i);
                end
                
                if ( self.sameEdges ( self.points_of_bdryEdges(i,:) , [ID2,ID3] ) ) 
                    self.points_of_bdryEdges(i,:) = [ID2,IDn{2}];
                    self.no_of_bdryEdges = self.no_of_bdryEdges + 1 ;
                    self.points_of_bdryEdges(self.no_of_bdryEdges,:) = [IDn{2},ID3];
                    self.bdrymarker_of_edges(self.no_of_bdryEdges) = self.bdrymarker_of_edges(i);
                end
                
                if ( self.sameEdges ( self.points_of_bdryEdges(i,:) , [ID3,ID4] ) )
                    self.points_of_bdryEdges(i,:) = [ID3,IDn{3}];
                    self.no_of_bdryEdges = self.no_of_bdryEdges + 1 ;
                    self.points_of_bdryEdges(self.no_of_bdryEdges,:) = [IDn{3},ID4];
                    self.bdrymarker_of_edges(self.no_of_bdryEdges) = self.bdrymarker_of_edges(i);
                end
                
                if ( self.sameEdges ( self.points_of_bdryEdges(i,:) , [ID4,ID1] ) )
                    self.points_of_bdryEdges(i,:) = [ID4,IDn{4}];
                    self.no_of_bdryEdges = self.no_of_bdryEdges + 1 ;
                    self.points_of_bdryEdges(self.no_of_bdryEdges,:) = [IDn{4},ID1];
                    self.bdrymarker_of_edges(self.no_of_bdryEdges) = self.bdrymarker_of_edges(i);
                end                                                
            end
        end
        
        function [self,newEdges] = InflationLayerForMarker(self,marker)

            counter = 0;
            for i = 1 : self.no_of_bdryEdges
                if ( self.bdrymarker_of_edges(i) ~= marker )
                    continue;
                end
%
%               Just edges with the desired marker pass
%               ---------------------------------------
                edge = self.points_of_bdryEdges(i,:);
%
%               Look for the edge in the element list
%               -------------------------------------
                for el = 1 : self.no_of_elements
                    if ( any( self.points_of_elements(el,:) == edge(1) ) && ...
                         any( self.points_of_elements(el,:) == edge(2) ) )
                            break;
                    end
                    
                    if ( el == self.no_of_elements)
                        error('Element for marker was not found!');
                    end
                end
%
%               Once the element is found, divide it according to the edge
%               position
%               -----------------------------------------------------------
                element = self.points_of_elements(el,:);
                type = self.edgeType(edge , element );
                
                ID1 = element(1);
                ID2 = element(2);
                ID3 = element(3);
                ID4 = element(4);
                
                x1 = [self.xp(ID1) , self.yp(ID1) ];
                x2 = [self.xp(ID2) , self.yp(ID2) ];
                x3 = [self.xp(ID3) , self.yp(ID3) ];
                x4 = [self.xp(ID4) , self.yp(ID4) ];
                
                if ( type == 1 ) 
%
%                   Check whether the nodes should be created
%                   -----------------------------------------
                    xn{1} = 0.5*(x1+x4);
                    xn{2} = 0.5*(x2+x3);
                    for k = 1 : 2
                        create = true;
                        for j = 1 : self.no_of_nodes
                            x = [self.xp(j),self.yp(j)];
                    
                            if ( self.almostEqual(x,xn{k}) )
                                IDn{k} = j;
                                create = false;
                                break;
                            end
                        end
                
                        if (create)
                            self.no_of_nodes = self.no_of_nodes + 1 ;
                            IDn{k} = self.no_of_nodes;
                            self.xp(IDn{k}) = xn{k}(1);
                            self.yp(IDn{k}) = xn{k}(2);
                        end
                    end     
%
%                   Add the elements
%                   ----------------
                    self.points_of_elements(el,:) = [ID1,ID2,IDn{2},IDn{1}];
                    Nel = self.no_of_elements;
                    self.points_of_elements(Nel+1,:) = [IDn{1},IDn{2},ID3,ID4];
                    self.no_of_elements = self.no_of_elements + 1 ;                    

                    
                elseif ( type == 2 )
%
%                   Check whether the nodes should be created
%                   -----------------------------------------
                    xn{1} = 0.5*(x3+x4);
                    xn{2} = 0.5*(x1+x2);
                    for k = 1 : 2
                        create = true;
                        for j = 1 : self.no_of_nodes
                            x = [self.xp(j),self.yp(j)];
                    
                            if ( self.almostEqual(x,xn{k}) )
                                IDn{k} = j;
                                create = false;
                                break;
                            end
                        end
                
                        if (create)
                            self.no_of_nodes = self.no_of_nodes + 1 ;
                            IDn{k} = self.no_of_nodes;
                            self.xp(IDn{k}) = xn{k}(1);
                            self.yp(IDn{k}) = xn{k}(2);
                        end
                    end     
%
%                   Add the elements
%                   ----------------
                    self.points_of_elements(el,:) = [ID1,IDn{1},IDn{2},ID4];
                    Nel = self.no_of_elements;
                    self.points_of_elements(Nel+1,:) = [IDn{1},ID2,ID3,IDn{2}];
                    self.no_of_elements = self.no_of_elements + 1 ;                          
                end

                 counter = counter + 1 ;
                 newEdges(counter,:) = [IDn{1},IDn{2}];
            end
        end        
%
%/////////////////////////////////////////////////////////////////////////////////////
%
%               SET CURVES PROCEDURES
%               ---------------------
%/////////////////////////////////////////////////////////////////////////////////////
%
        function self = SetCircularCurveForMarker(self,marker,xc,R)
%
%           Loop in all boundary edges
%           --------------------------
            for i = 1 : self.no_of_bdryEdges
                if ( self.bdrymarker_of_edges(i) == marker ) 
                    self = self.SetCircularCurve(self.points_of_bdryEdges(i,:),xc,R);
                end
            end
        end
        
        function self = SetCircularCurve(self,nodes,xc,R)
%
%           First, force the nodes to belong to the curve
%           ---------------------------------------------
            x1 = self.xp(nodes(1));
            y1 = self.yp(nodes(1));
            theta1 = atan2(y1-xc(2),x1-xc(1));
            
            self.xp(nodes(1)) = R * cos(theta1)+xc(1);
            self.yp(nodes(1)) = R * sin(theta1)+xc(2);
            
            x2 = self.xp(nodes(2));
            y2 = self.yp(nodes(2));
            theta2 = atan2(y2-xc(2),x2-xc(1));
            
            self.xp(nodes(2)) = R * cos(theta2)+xc(1);
            self.yp(nodes(2)) = R * sin(theta2)+xc(2);
%
%           Next, construct a new curved edge selfect
%           ----------------------------------------
            self.no_of_curvedEdges = self.no_of_curvedEdges + 1 ;
            self.curvedEdges{self.no_of_curvedEdges} = CurvedEdge(nodes , self.xp(nodes) , self.yp(nodes) , self.curvesOrder);
            current = self.no_of_curvedEdges;
%
%           Build the curve
%           ---------------
            self.curvedEdges{current} = self.curvedEdges{current}.SetCircularEdge( R , xc );
                                               
        end
        
        function value = get.sj(self)
            value = cos(pi-pi*double(0:self.curvesOrder)/double(self.curvesOrder));
        end        
        
        function value = get.points_of_curvedEdges(self)
            value = zeros(self.no_of_curvedEdges , 2);
            
            for i = 1 : self.no_of_curvedEdges
                value(i,:) = self.curvedEdges{i}.nodes;
            end
        end
        
        function value = get.x_curvilinear_coords(self)
            value = zeros(self.no_of_curvedEdges , self.curvesOrder + 1 );
            
            for i = 1 : self.no_of_curvedEdges
                value(i,:) = self.curvedEdges{i}.x_coord;
            end
        end
        
        function value = get.y_curvilinear_coords(self)
            value = zeros(self.no_of_curvedEdges , self.curvesOrder + 1 );
            
            for i = 1 : self.no_of_curvedEdges
                value(i,:) = self.curvedEdges{i}.y_coord;
            end
        end
        
        function Plot(self)
%
%           Create new figure
%           -----------------
            h = figure; hold on; axis equal;
            dcm = datacursormode(h);
            dcm.UpdateFcn = @PointDataCursor;
%
%           Plot elements
%           -------------
            colors = {'-m','-g','-b','-r','-c','-y'};
            pZones = self.pRefinement_zones;
            for i = 1 : self.no_of_elements
                patch(self.xp(self.points_of_elements(i,:)),self.yp(self.points_of_elements(i,:)),colors{pZones(i)+1},'FaceAlpha',0.2);
                xc(i) = mean(self.xp(self.points_of_elements(i,:)));
                yc(i) = mean(self.yp(self.points_of_elements(i,:)));
            end
%
%           Plot nodes
%           ----------
            scatter(self.xp,self.yp,40,'ok','LineWidth',2,'MarkerFaceColor',[1,1,1]);
%
%           Plot element centers
%           --------------------
          %  scatter(xc,yc,40,'ok','LineWidth',2,'MarkerFaceColor',[0,0,0]);
%
%           Plot curved edges
%           -----------------
            for i = 1 : self.no_of_curvedEdges
                self.curvedEdges{i}.Plot;
            end
            grid off
        end

        function Save(self,fileName)
            
            ncid        = netcdf.create(fileName,'CLOBBER');
            Nel_ID      = netcdf.defDim(ncid,'no_of_elements',self.no_of_elements);
            Np_ID       = netcdf.defDim(ncid,'no_of_nodes' , self.no_of_nodes);
            Nedges_ID   = netcdf.defDim(ncid,'no_of_bdryedges',self.no_of_bdryEdges);
            if ( self.no_of_curvedEdges ~= 0 )
                Ncurvilinear_ID   = netcdf.defDim(ncid,'no_of_curvilinearedges',self.no_of_curvedEdges);
                Np1_ID            = netcdf.defDim(ncid,'Np1',self.curvesOrder+1);
            end
            
            N1_ID       = netcdf.defDim(ncid,'one',1);
            N2_ID       = netcdf.defDim(ncid,'two',2);
            N3_ID       = netcdf.defDim(ncid,'three',3);
            N4_ID       = netcdf.defDim(ncid,'four',4);    
            NChar_ID    = netcdf.defDim(ncid,'character_length',20);
            Nmarkers_ID = netcdf.defDim(ncid,'no_of_markers',5);
            NPzones_ID  = netcdf.defDim(ncid,'no_of_pRefinementzones',self.no_of_pRefinementzones);
%           *********
%           Variables
%           *********
            Nelements   = netcdf.defVar(ncid,'points_of_quads','NC_INT',[N4_ID , Nel_ID]);
            Ncoord      = netcdf.defVar(ncid,'points','NC_DOUBLE',[N2_ID , Np_ID]);
            edges_ID    = netcdf.defVar(ncid,'points_of_bdryedges','NC_INT',[N2_ID , Nedges_ID]);
            bdrymarker_ID = netcdf.defVar(ncid,'bdrymarker_of_edges','NC_INT',Nedges_ID);
            if ( self.no_of_curvedEdges ~= 0 ) 
                curvilinear_ID = netcdf.defVar(ncid,'curvilinear_edges','NC_INT',[N2_ID , Ncurvilinear_ID]);
                x_curvilinear_ID = netcdf.defVar(ncid , 'x_curvilinear_edges','NC_DOUBLE',[Np1_ID , Ncurvilinear_ID]);
                y_curvilinear_ID = netcdf.defVar(ncid , 'y_curvilinear_edges','NC_DOUBLE',[Np1_ID , Ncurvilinear_ID]);
            end
            pRefinementZones_ID = netcdf.defVar(ncid,'pRefinement_zones','NC_INT',Nel_ID);
    
            marker1_ID = netcdf.defVar(ncid , 'marker1','NC_CHAR',NChar_ID);
            marker2_ID = netcdf.defVar(ncid , 'marker2','NC_CHAR',NChar_ID);
            marker3_ID = netcdf.defVar(ncid , 'marker3','NC_CHAR',NChar_ID);
            marker4_ID = netcdf.defVar(ncid , 'marker4','NC_CHAR',NChar_ID);
            marker5_ID = netcdf.defVar(ncid , 'marker5','NC_CHAR',NChar_ID);
    
%           ===================    
            netcdf.endDef(ncid);    
%           ===================
    
            netcdf.putVar(ncid , Nelements , self.points_of_elements');
            netcdf.putVar(ncid , Ncoord , [self.xp,self.yp]')
            netcdf.putVar(ncid , edges_ID , self.points_of_bdryEdges');
            netcdf.putVar(ncid , bdrymarker_ID , self.bdrymarker_of_edges);            
            if ( self.no_of_curvedEdges ~= 0 ) 
                netcdf.putVar(ncid , curvilinear_ID , self.points_of_curvedEdges');
                netcdf.putVar(ncid , x_curvilinear_ID , self.x_curvilinear_coords');
                netcdf.putVar(ncid , y_curvilinear_ID , self.y_curvilinear_coords');
            end
            netcdf.putVar(ncid , pRefinementZones_ID , self.pRefinement_zones);
            netcdf.putVar(ncid , marker1_ID , ['bottom',blanks(20-6)]);
            netcdf.putVar(ncid , marker2_ID , ['top',blanks(20-3)]);
            netcdf.putVar(ncid , marker3_ID , ['left',blanks(20-4)]);
            netcdf.putVar(ncid , marker4_ID , ['right',blanks(20-5)]);
            netcdf.putVar(ncid , marker5_ID , ['circle',blanks(20-6)]);
            netcdf.close(ncid);
            
            
        end
        
        function val = get.pRefinement_zones(self)
            val = zeros(1,self.no_of_elements);    % Set the default value
            
            for i = 1 : self.no_of_elements
                xc = self.ElementCentroid(i);
                
                for j = 1 : self.no_of_pRefinementzones
                    if ( self.p_refinementZoneConditions{j}(xc) ) 
                        val(i) = j;
                        break;
                    end
                end
            end
            
        end        
        
        function val = get.no_of_pRefinementzones(self)
            val = length(self.p_refinementZoneConditions);
        end
        
        function xc = ElementCentroid(self,eID)
            x1 = [self.xp(self.points_of_elements(eID,1)) , self.yp(self.points_of_elements(eID,1)) ];
            x2 = [self.xp(self.points_of_elements(eID,2)) , self.yp(self.points_of_elements(eID,2)) ];
            x3 = [self.xp(self.points_of_elements(eID,3)) , self.yp(self.points_of_elements(eID,3)) ];
            x4 = [self.xp(self.points_of_elements(eID,4)) , self.yp(self.points_of_elements(eID,4)) ];            
           
            xc = 0.25*(x1+x2+x3+x4);
        end
        
    end     %/*  Methods  */
    
    methods ( Static )
        function val = almostEqual(a,b)
            
           for i = 1 : length(a)
              if ( abs(a(i)-b(i)) < 10*eps(1.0) )
                val = true;
              else
                val = false;
                return;
              end
           end
        end
        
        function val = sameEdges(ed1 , ed2)
            
            if ( (ed1(1) == ed2(1)) && (ed1(2) == ed2(2)) )
                val = true;
                
            elseif ( (ed1(1) == ed2(2)) && (ed1(2) == ed2(1)) )
                val = true;
                
            else
                val = false;
                
            end
            
        end
        
        function val = edgeType(ed,element)
            
            if ( MeshFile.sameEdges(ed,element([1,2])) )
                val = 1;
                
            elseif ( MeshFile.sameEdges(ed,element([2,3])))
                val = 2;
                
            elseif ( MeshFile.sameEdges(ed,element([3,4]))) 
                val = 1;
                
            elseif ( MeshFile.sameEdges(ed,element([4,1])))
                val = 2;
                
            end
            
        end
        

    end
end %/*  Class  */

