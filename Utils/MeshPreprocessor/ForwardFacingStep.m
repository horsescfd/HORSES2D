%
clear all
clc
close all
%
addpath('./Functions');

file = '../../TestCases/SupersonicTests/ForwardFacingStep/MESH/forward-facing-step.msh';

mesh = importGMSH(file);
mesh.Save('../../TestCases/SupersonicTests/ForwardFacingStep/MESH/forward-facing-step.HiOMesh');
