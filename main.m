% 1D Thermal Stress Finite Element Program

% clean slate
close all
clear
clc

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

% extracting node and elelment information
node = nodeData(:,1);           % node number
xPos = nodeData(:,2);           % node x position (m)
conType = nodeData(:,3);        % node condition type
conVal = nodeData(:,4);         % node condition value (N) or (m)
element = elementData(:,1);     % element number
node1 = elementData(:,2);       % node 1 of element
node2 = elementData(:,3);       % node 2 of element
D = elementData(:,4);           % element diameter (m)
E = elementData(:,5);           % element elastic modulus (N/m^2)
temp = elementData(:,6);        % element temperature (deg C)
alpha = elementData(:,7);       % element thermal expansion coefficient (1/deg C)

% setting up stiffness matrix and body force vector
K = zeros(numNodes);            % initializing global stiffness matrix
G = zeros(numNodes, 1);         % initializing global body force vector
Kindv = [1 -1;-1 1];            % individual element stiffness matrix form
Gindv = [-1;1];                 % individual element body force vector form
for n = 1:numElements
    n1 = node1(n);
    n2 = node2(n);
    L(n) = abs(xPos(n2)-xPos(n1));  % element length
    A(n) = pi*(D(n)/2)^2;           % element area
    dTemp(n) = temp(n) - refTemp;   % element temperature change
    if xPos(n1)>xPos(n2) % fixing node position issues
        G([n2,n1],:) = G([n2,n1],:) + E(n)*A(n)*alpha(n)*dTemp(n)*Gindv;
        K([n2,n1],[n2,n1]) = K([n2,n1],[n2,n1]) + E(n)*A(n)/L(n)*Kindv;
    else % fixing node position issues
        G([n1,n2],:) = G([n1,n2],:) + E(n)*A(n)*alpha(n)*dTemp(n)*Gindv;
        K([n1,n2],[n1,n2]) = K([n1,n2],[n1,n2]) + E(n)*A(n)/L(n)*Kindv;
    end
end

% setting up field and boundary variables
for n = 1:numNodes
    if conType(n) == 0 % 0 for force and 1 for displacement
        F(n) = conVal(n);
        u(n) = NaN;
    else
        F(n) = NaN;
        u(n) = conVal(n);
    end
end

% call function to perform modify in place operations
[u,F] = modInPlace(K,G,u',F');

% calculating stresses and strains in each element
for n = 1:numElements
    n1 = node1(n);
    n2 = node2(n);
    if xPos(n1)>xPos(n2) % fixing node position issues
        du(n) = u(n1)-u(n2); % calculating expansion
    else
        du(n) = u(n2)-u(n1); % calculating expansion
    end
    strain(n) = du(n)/L(n); % element strain
    stress(n) = E(n)*strain(n)-E(n)*alpha(n)*dTemp(n); % element stress
end

% formating data
format longE

% printing results
nodeResults = [node,u,F];
elementResults = [element,strain',stress'];
fprintf('%20s %20s %20s \n', 'node','displacement (m)', 'external force (N)');
fprintf('%20d %20d %20d \n', nodeResults.');
fprintf('\n')
fprintf('%20s %20s %20s \n', 'element','strain (m/m)', 'stress (Pa)');
fprintf('%20d %20d %20d \n', elementResults.');


%% functions

% modify in place
function [u,F] = modInPlace(K,G,u,F)
    Ki = K; % original stiffness matrix
    Gi = G; % original body force vector
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
    u = K\(F+G); % solving for displacement vector
    F = Ki*u-Gi; % solving for force vector
end