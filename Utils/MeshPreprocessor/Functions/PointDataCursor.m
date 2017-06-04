function output_txt = PointDataCursor(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

% Import x and y
x = get(get(event_obj,'Target'),'XData');
y = get(get(event_obj,'Target'),'YData');

% Find index
index_x = find(x == pos(1));
index_y = find(y == pos(2));
index = intersect(index_x,index_y);

% Set output text
output_txt = {['X: ',num2str(pos(1),4)], ...
              ['Y: ',num2str(pos(2),4)], ...
              ['Index: ', num2str(index)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

fid = fopen('./Data/Indexes.dat','a');

fprintf(fid,'%d\n',index);
fclose(fid);