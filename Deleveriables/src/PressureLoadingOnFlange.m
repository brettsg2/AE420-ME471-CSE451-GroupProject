%========================================================================
% This is the script that will run the cases we care about
%========================================================================

clc
close all
clear 
currentFolder = pwd;
file    =  strcat(currentFolder, '\2DSPIRAL_HALFWOCAP.msh');
% these constants are their natural SI value then multiplied/divided
% to use a consistent set of units (the numbers at this point are simply
% made up/plausable)
E_fl = 204000; % from abaqus
E_gsk = 18000; % from abaqus
poison = 0.3;
poison_gsk = 0.26;
% old_pressure = 310.3;
pressure = 199.544*8; % pressure times distance pressure is being applied on
% read in the mesh
originalMesh = ReadInMeshInfo(file, 100, []);
mesh = originalMesh;
% number of nodes gets used a lot, store it
N_n=originalMesh.NumNodes;

% get nodes that are fixed (just making up numbers here)
%fixedXNodes = Node.FindNodesInRange(originalMesh.Nodes, 2.98, 4.02, 0, 5); % holding inside wall constant
% and where we stop modeling the pipe (cause why not)
fixedXNodes = Node.FindNodesInRange(originalMesh.Nodes, 0.144*100, 0.1465*100, -0.25*100, 0.02*100);
fixedYNodes = Node.FindNodesInRange(originalMesh.Nodes, 0.16*100, 0.184*100, -0.003*100, -0.002*100);

% get nodes where there is a pressure
pressureNodes = Node.FindNodesInRange(originalMesh.Nodes, (0.305-0.017245)*100, (0.305-0.011255)*100, (-0.088)*100, (-0.0875)*100);
pressureNodes = [pressureNodes, Node.FindNodesInRange(originalMesh.Nodes, (0.305-0.06535)*100, (0.305-0.05736)*100, (-0.088)*100, (-0.0875)*100)];

% do it all!
numSteps = 10;
finalAnswer = zeros(1,N_n*2);
pressures = [pressure];
maxDisplacements = [1,length(pressures)];
femSystems(1:length(numSteps)) = FemSystem(1,2,3, 4);
pc = 0;
%for pressureSte = 1:numSteps
for pc = 1:numSteps
    %pc=pc+1;
    % always adding the delta of the pressure 
    pressureToLoad = pressure * 1/numSteps;
    femSystem = ApplyStressToMesh(mesh, E_fl, poison, pressureToLoad, pressureNodes, fixedXNodes, fixedYNodes);
    femSystems(pc)=femSystem;
    mesh = ReadInMeshInfo(file, 100, femSystem.D);
    finalAnswer = finalAnswer + femSystem.D;
    maxDisplacements(pc) = max(abs(finalAnswer));
end

maxDisplacements
% post processing
magnitudes = 1:N_n;
xs = 1:N_n;
ys = 1:N_n;
stresses = zeros(3,N_n);
strains = zeros(3,N_n);
stressesPerNode = zeros([1,N_n]);
strainsPerNode = zeros(1,N_n);
elc = 0;
for el = femSystem.Elements
    elc=elc+1;
    x1 = finalAnswer((el.Node1.Index-1)*2+1);
    y1 = finalAnswer((el.Node1.Index-1)*2+2);
    x2 = finalAnswer((el.Node2.Index-1)*2+1);
    y2 = finalAnswer((el.Node2.Index-1)*2+2);
    x3 = finalAnswer((el.Node3.Index-1)*2+1);
    y3 = finalAnswer((el.Node3.Index-1)*2+2);
    
    thisStrain = el.ComputeStrain([x1, y1, x2, y2, x3, y3]);
    strains(1,elc) = thisStrain(1);
    strains(2,elc) = thisStrain(2);
    strains(2,elc) = thisStrain(2);
    thisStress = el.ComputeStress(thisStrain);
    stresses(1, elc) = thisStress(1);
    stresses(2, elc) = thisStress(2);
    stresses(3, elc) = thisStress(3);

    stressesPerNode(el.Node1.Index) = sqrt(thisStress(1)^2 + stressesPerNode(el.Node1.Index)^2);
    stressesPerNode(el.Node2.Index) =sqrt(thisStress(2)^2 + stressesPerNode(el.Node2.Index)^2);
    stressesPerNode(el.Node3.Index) =sqrt(thisStress(3)^2 + stressesPerNode(el.Node3.Index)^2);
    
    strainsPerNode(el.Node1.Index) =sqrt(thisStrain(1)^2 + strainsPerNode(el.Node1.Index)^2);
    strainsPerNode(el.Node2.Index) =sqrt(thisStrain(2)^2 + strainsPerNode(el.Node2.Index)^2);
    strainsPerNode(el.Node3.Index) = sqrt(thisStrain(3)^2 + strainsPerNode(el.Node3.Index)^2);
%I Don't know how to properly compare the stresses and strains here with
%the von mises stress in Abaqus...
%      stressesPerNode(el.Node1.Index) = abs(thisStress(1)) + stressesPerNode(el.Node1.Index);
%      stressesPerNode(el.Node2.Index) =abs(thisStress(2)) + stressesPerNode(el.Node2.Index);
%      stressesPerNode(el.Node3.Index) =abs(thisStress(3)) + stressesPerNode(el.Node3.Index);
%      
%      strainsPerNode(el.Node1.Index) =abs(thisStrain(1)) + strainsPerNode(el.Node1.Index);
%      strainsPerNode(el.Node2.Index) =abs(thisStrain(2)) + strainsPerNode(el.Node2.Index);
%      strainsPerNode(el.Node3.Index) =abs(thisStrain(3)) + strainsPerNode(el.Node3.Index);

end

for nC= 1:N_n
    x=finalAnswer((nC-1)*2+1);
    xs(nC) = x;
    y=finalAnswer((nC-1)*2+2);
    ys(nC) = y;
    magnitudes(nC) = sqrt(x*x+y*y);
end

figure(1)
trimesh(originalMesh.TwoDElements, originalMesh.TwoDNodes(:,1)*100,originalMesh.TwoDNodes(:,2)*100)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Original Mesh','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white');

figure(2)
trimesh(mesh.TwoDElements, mesh.TwoDNodes(:,1),mesh.TwoDNodes(:,2), stressesPerNode)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Stresses','fontsize',14)
fh = figure(2);
set(fh, 'color', 'white'); 

figure(3)
trimesh(mesh.TwoDElements, mesh.TwoDNodes(:,1),mesh.TwoDNodes(:,2), strainsPerNode)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Strains','fontsize',14)
fh = figure(3);
set(fh, 'color', 'white'); 

figure(4)
trimesh(mesh.TwoDElements, mesh.TwoDNodes(:,1),mesh.TwoDNodes(:,2), magnitudes*100)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('New Mesh','fontsize',14)
fh = figure(4);
set(fh, 'color', 'white'); 

figure(5)
exMesh = ReadInMeshInfo(file, 100, femSystem.D*100);
trimesh(exMesh .TwoDElements, exMesh.TwoDNodes(:,1),exMesh .TwoDNodes(:,2), magnitudes)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Exagerated Mesh','fontsize',14)
fh = figure(5);
set(fh, 'color', 'white'); 

%-------------------------------------------------------------------------