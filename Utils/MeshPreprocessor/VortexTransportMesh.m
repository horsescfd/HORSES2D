%
clear all
clc
close all
%
addpath('Functions');

% Generate a simple linear mesh
x = [-0.2,0.2,0.2,-0.2,-0.5,0.5,0.5,-0.5,-1,1,1,-1]';
y = [-0.2,-0.2,0.2,0.2,-0.5,-0.5,0.5,0.5,-1,-1,1,1]';

el = [1,2,3,4 ; ...
      5,6,2,1 ; ...
      2,6,7,3 ; ...
      4,3,7,8 ; ...
      1,4,8,5 ; ...
      9,10,6,5;...
      6,10,11,7;...
      8,7,11,12;...
      9,5,8,12];

  
  
mesh = MeshFile;
mesh = mesh.Load(el,x,y,[],[],10);

%
%   Subdivide all elements
for i = 1 : size(el,1)
    mesh = mesh.SubdivideElement(i);
end

elementsToSub = [12,10,2,4,6,8,13,14,15,19,20,21,25,26,27,31,32,33,37,42,5,16,17,22,34,35,9,28,29,36,96,91,1,11,11,134];

for i = 1 : length(elementsToSub)
    mesh = mesh.SubdivideElement(elementsToSub(i));
end
%
%   Set curves
interior_curve_points = [8    67    26    87     7   23      6    72    19    53     5    29     8];
for i = 1 : length(interior_curve_points)-1
    mesh = mesh.SetCircularCurve(interior_curve_points(i:i+1),[0,0],0.5);
end

if(true)
%TOP
circNodes = [-0.5,0.5;0.5,0.5];
straightNodes = [-1,1;1,1];
mesh = mesh.SetCircularStraightMappedCurve([39,69],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([69,37],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([37,102],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([102,36],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([70,66],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([66,68],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([68,100],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([100,101],[0,0],0.5,circNodes,straightNodes,0.25);

straightNodes = [-0.2,0.2;0.2,0.2];
mesh = mesh.SetCircularStraightMappedCurve([90,89],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([89,88],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([88,85],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([85,86],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([27,59],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([59,25],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([25,84],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([84,24],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([60,57],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([57,58],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([58,81],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([81,83],[0,0],0.5,circNodes,straightNodes,0.75);

%RIGHT
circNodes = [0.5,0.5;0.5,-0.5];
straightNodes = [1,1;1,-1];
mesh = mesh.SetCircularStraightMappedCurve([36,34],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([34,32],[0,0],0.5,circNodes,straightNodes,0.5);

straightNodes = [0.2,0.2;0.2,-0.2];
mesh = mesh.SetCircularStraightMappedCurve([86,128],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([128,126],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([126,124],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([124,73],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([24,22],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([22,20],[0,0],0.5,circNodes,straightNodes,0.5);

%LEFT
circNodes = [-0.5,0.5;-0.5,-0.5];
straightNodes = [-0.2,0.2;-0.2,-0.2];
mesh = mesh.SetCircularStraightMappedCurve([27,28],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([28,21],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([60,131],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([131,122],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([122,120],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([120,80],[0,0],0.5,circNodes,straightNodes,0.75);


%BOTTOM
circNodes = [-0.5,-0.5;0.5,-0.5];
straightNodes = [-0.2,-0.2;0.2,-0.2];
mesh = mesh.SetCircularStraightMappedCurve([80,78],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([78,77],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([77,75],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([75,76],[0,0],0.5,circNodes,straightNodes,0.75);
mesh = mesh.SetCircularStraightMappedCurve([56,52],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([52,54],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([54,71],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([71,73],[0,0],0.5,circNodes,straightNodes,0.25);
mesh = mesh.SetCircularStraightMappedCurve([21,55],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([55,18],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([18,74],[0,0],0.5,circNodes,straightNodes,0.5);
mesh = mesh.SetCircularStraightMappedCurve([74,20],[0,0],0.5,circNodes,straightNodes,0.5);

end

mesh.mode = 'Elements';
mesh.p_refinementElements = zeros(1,mesh.no_of_elements);
mesh.p_refinementElements([3,18]) = 1;
mesh.p_refinementElements([7,30]) = 2;
mesh.p_refinementElements([23,24]) = 3;
%
%   Boundary markers
mesh.no_of_bdryEdges = 16;
mesh.points_of_bdryEdges = [ 9,62;62,31;31,92;92,10;10,144;144,35;35,148;148,11;11,105;105,38;38,108;108,12;12,151;151,41;41,142;142,9];
mesh.bdrymarker_of_edges = [1,1,1,1,4,4,4,4,2,2,2,2,3,3,3,3];

%mesh.Plot;
mesh.Save('VortexTransport.HiOMesh');