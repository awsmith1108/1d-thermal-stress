% 1D Thermal Stress Finite Element Program

% clean slate
close all
clear
clc

% formating data
format longE

% loading data
fileIn  = 'sample.dat';
fileID = fopen(fileIn,'r');
   
% reading data from input file
refData = fscanf(fileID,'%g',[1 3]); 
numNodes = refData(1);
numElements = refData(2);
refTemp = refData(3);
nodeData = fscanf(fileID,'%g',[4 numNodes])';      
elementData = fscanf(fileID,'%g',[7 numElements])';

% extracting node information
node = nodeData(:,1);
xPosition = nodeData(:,2);
conditionType = nodeData(:,3);
conditionVal = nodeData(:,4);
element = elementData(:,1);
node1 = elementData(:,2);
node2 = elementData(:,3);
D = elementData(:,4);
E = elementData(:,5);
temp = elementData(:,6);
alpha = elementData(:,7);

% setting up stiffness matrix
K = zeros(numNodes);
G = zeros(numNodes, 1);
Kindv = [1 -1;-1 1];
Gindv = [-1;1];
for n = 1:numElements
    n1 = node1(n);
    n2 = node2(n);
    L(n) = abs(xPosition(n2)-xPosition(n1));
    A(n) = pi*(D(n)/2)^2;
    dTemp(n) = temp(n) - refTemp;
    if xPosition(n1)>xPosition(n2)
        G([n2,n1],:) = G([n2,n1],:) + E(n)*A(n)*alpha(n)*dTemp(n)*Gindv;
        K([n2,n1],[n2,n1]) = K([n2,n1],[n2,n1]) + E(n)*A(n)/L(n)*Kindv;
    else
        G([n1,n2],:) = G([n1,n2],:) + E(n)*A(n)*alpha(n)*dTemp(n)*Gindv;
        K([n1,n2],[n1,n2]) = K([n1,n2],[n1,n2]) + E(n)*A(n)/L(n)*Kindv;
    end
end

% setting up field and boundary variables
for n = 1:numNodes
    if conditionType(n) == 0
        F(n) = conditionVal(n);
        u(n) = NaN;
    else
        F(n) = NaN;
        u(n) = conditionVal(n);
    end
end

% call function to perform modify in place operations
[u,F] = modInPlace(K,G,u',F');

% calculating stresses and strains in each element
for n = 1:numElements
    n1 = node1(n);
    n2 = node2(n);
    if xPosition(n1)>xPosition(n2)
        du(n) = u(n1)-u(n2);
    else
        du(n) = u(n2)-u(n1);
    end
    strain(n) = du(n)/L(n);
    stress(n) = E(n)*strain(n)-E(n)*alpha(n)*dTemp(n);
end

% Printing results
nodeResults = [node,u,F];
elementResults = [element,strain',stress'];
fprintf('%20s %20s %20s \n', 'node','displacement (m)', 'external force (N)');
fprintf('%20d %20d %20d \n', nodeResults.');
fprintf('\n')
fprintf('%20s %20s %20s \n', 'element','strain (m/m)', 'stress (Pa)');
fprintf('%20d %20d %20d \n', elementResults.');


%% functions

function [u,F] = modInPlace(K,G,u,F)
    Ki = K;
    Gi = G;
    for n = 1:length(u)
        if isnan(u(n)) == false
            kn = K(n,n);
            F(n) = 0;
            G(n) = 0;
            F = F - u(n)*K(:,n);
            F(n) = -F(n);
            K(n,:) = 0;
            K(:,n) = 0;
            K(n,n) = kn;
        end
    end
    u = K\(F+G);
    F = Ki*u-Gi;
end

                                                   
   
