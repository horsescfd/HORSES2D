function mesh = importGMSH(file)

EDGE_TYPE = 1;
QUAD_TYPE = 3;

fid = fopen(file);

%   1/ Read until $Nodes is find.
%      --------------------------
       flag = true;
       while (flag)
         str = fgetl(fid);
         
         if ( strcmpi(str , '$Nodes') ) 
             flag = false;
         end
       end 

%      Read nodes:
       no_of_nodes = str2num(fgetl(fid));
      
       aux = textscan(fid,'%d %f %f %f',no_of_nodes);
       pcoord(:,1) = aux{2};
       pcoord(:,2) = aux{3};
       
fclose(fid);       
fid = fopen(file);
%   2/ Read until $Elements is found.
%      -----------------------------
        flag = true;
        while (flag)
          str = fgetl(fid);
          
          if ( strcmpi(str , '$Elements') ) 
              flag = false;
          end
        end 
 
 %      Read elements:
        no_of_entities = str2num(fgetl(fid));
        no_of_elements = 0;
        no_of_bdryedges    = 0;
        
        aux = textscan(fid ,'%d %d %d %d %d %d %d %d %d',no_of_entities);
        
        type = aux{2};
        location = aux{4};
        IDs  = [aux{6:9}];
        
        for ent = 1 : no_of_entities
            
            if ( type(ent) == QUAD_TYPE )
                no_of_elements = no_of_elements + 1;
                elements(no_of_elements,1:4) = IDs(ent,1:4);
            
            elseif ( type(ent) == EDGE_TYPE )
                no_of_bdryedges = no_of_bdryedges + 1;
                bdryedges(no_of_bdryedges,1:2) = IDs(ent,1:2);
                marker_of_bdryedges(no_of_bdryedges) = location(ent);
                
            end
       
            
        end
        
fclose(fid);
fid = fopen(file);
 %   3/ Read until $PhysicalNames is found.
 %      ----------------------------------
        flag = true;
         while (flag)
           str = fgetl(fid);
           
           if ( strcmpi(str , '$PhysicalNames') ) 
               flag = false;
           end
         end 
         
         no_of_markers = str2num(fgetl(fid))-1;     % Always set the last marker to "interior"
         
         aux = textscan(fid,'%d %d %s',no_of_markers);
         markerNames = aux{3};
         for i = 1 : length(markerNames)
             markerNames{i} = markerNames{i}(2 : length(markerNames{i})-1);
         end
       
 %   Close file       
     fclose(fid);
     
%   Save the mesh in a MeshFile object.
%   ----------------------------------
    mesh = MeshFile;
    mesh = mesh.Load(elements,pcoord(:,1),pcoord(:,2),bdryedges,marker_of_bdryedges,10);
    mesh.markerNames = markerNames;
end
      