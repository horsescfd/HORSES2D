%
clear all
clc
close all
%
addpath('./Functions');
%
load('./Data/Cylinder.mat');

mesh = MeshFile;
mesh = mesh.Load(el,x,y,points_of_bdryface,bdryMarkers,P);

load('./Data/ElementsToSubdivide');
for i = 1 : length(elements_to_subdivide)
    mesh = mesh.SubdivideElement(elements_to_subdivide(i));
end

[mesh,inf1] = mesh.InflationLayerForMarker(5);
[mesh,inf2] = mesh.InflationLayerForMarker(5);
[mesh,inf3] = mesh.InflationLayerForMarker(5);


mesh = mesh.SetCircularCurveForMarker(5,[0,0],0.5);
for i = 1 : size(inf3,1)
    mesh = mesh.SetCircularCurve(inf3(i,:),[0,0],0.52);
end

for i = 1 : size(inf2,1)
    mesh = mesh.SetCircularCurve(inf2(i,:),[0,0],0.55);
end

for i = 1 : size(inf1,1)
    mesh = mesh.SetCircularCurve(inf1(i,:),[0,0],0.6);
end

load('./Data/InteriorCurvePoints');
for i = 1 : size(interior_curve_points,1)
    mesh = mesh.SetCircularCurve(interior_curve_points(i,:),[0,0],0.65);
end

load('./Data/ExteriorCurvePoints');
for i = 1 : size(exterior_curve_points,1)
    mesh = mesh.SetCircularCurve(exterior_curve_points(i,:),[0,0],0.75);
end
%
%   Set p-refinement zones
%   ----------------------
    mesh.p_refinementZoneConditions{1} = @(x)(sqrt(x(1).^2 + x(2).^2) < 0.55);
    mesh.p_refinementZoneConditions{2} = @(x)(sqrt(x(1).^2 + x(2).^2) < 0.75);
    mesh.p_refinementZoneConditions{3} = @(x)(gt(x(1),0) && lt(x(2),0.6) && gt(x(2),-0.6));
    mesh.p_refinementZoneConditions{4} = @(x)((gt(x(1),0) && lt(x(2),1.2) && gt(x(2),-1.2)) || (sqrt(x(1).^2 + x(2).^2) < 1.2));

%   Plot the mesh
%   -------------
    mesh.Plot;

%
%   Save the mesh
%   -------------
    mesh.Save('CylinderRefined.HiOMesh');
