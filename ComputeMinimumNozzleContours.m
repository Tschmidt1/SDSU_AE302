%{ 
This script computes the (x,y) coordinates of a minimum length supersonic
nozzle contour using the method of characteristics. This script is intended
to be used in an educational setting and should NOT be used as an actual
engineering tool. 

This is a 2D computation, so variables such as "AStar" and "AExit" are not
true area's, they are unit areas (a.k.a. lengths). 

The coordinate system is set such that (0,0) is located at the nozzle
throat, along the nozzle centerline. 

The necessary input variables are:
    Me      -- Mach number of exit flow
    AStar   -- "area" at the throat. For convienence, this is 1.0 by default
    numChar -- number of characteristic linesto compute

More "realistic" design parameters may also be specified, and used to
compute Me and/or AStar instead of directly defining them. Such variables
may be but are not limited to:
    Pa   -- ambient pressure
    Po   -- chamber pressure
    To   -- chamber temperature
    mDot -- mass flow rate
Note: The use will need to implement this themselves, this script is as
barebones as possible to aid basic understanding of the method of
characteristics.

%}

%% Basic Setup %%

close all;
clear all;
clc;

% Add path to Matlab that are called in this script here


% Set default interpreters
set(0,'defaulttextinterpreter','latex')
set(0,'defaulttextfontname','times')
set(0,'defaultaxesfontname','times')

warning('OFF','ALL');

%% Basic Variables %%

% Geometric variables
AStar = 1.0;    % "area" of the throat

% Flow variables
Me = 2.4;       % Mach number at exit plane
G  = 1.4;       % Ratio of specific heats
R  = 287;       % Gas constant
% Po =            % Chamber pressure
% To =            % Chamber temperature
% Ro = Po/R/To;   % Chamber density

% Method of Characteristic specific variables
numCharLines = 16;    % Number of characteristic lines to compute

% Plotting variables
fontSize        = 18;
fontName        = 'times';
lineWidth       = 2;
markerSize      = 6;
charLineStyle   = 'k--';
nozzleLineStyle = 'k-';
charNodeStyle   = 'ro';

% Figure setup variables
figPos = [680 338 1030 640];

%% Set Up -- Create Data Structures %% 

% Set the initial number of nodes per characteristic line, including the
% node located on the nozzle wall
nodesPerCharLine = numCharLines + 1;

% Compute the total number of characteristic nodes there will be, including
% those along the nozzle wall
numNodesTot = cumsum(1:numCharLines) + numCharLines;
numNodesTot = numNodesTot(end);

% Create data structures for nodes generated at each characteristic line
% intersection
for i = 1:numNodesTot
    % Create structure for C- characteristic lines
    Node{i}.ID      = zeros(1,1);
    Node{i}.KMinus  = zeros(1,1);
    Node{i}.KPlus   = zeros(1,1);
    Node{i}.Theta   = zeros(1,1);
    Node{i}.Nu      = zeros(1,1);
    Node{i}.Mach    = zeros(1,1);
    Node{i}.Mu      = zeros(1,1);
    Node{i}.Coords  = zeros(1,2);
    Node{i}.R       = zeros(1,1);
    Node{i}.U       = zeros(1,1);
    Node{i}.V       = zeros(1,1);
    Node{i}.P       = zeros(1,1);
    Node{i}.T       = zeros(1,1);
%     Node{i}.RhoStag = Ro;
%     Node{i}.TStag   = To;
%     Node{i}.PStag   = Po;
end

% Create data structure for nozzle contour
NozzleNodes = cell(numCharLines,1);

% Initialize an array that will store the IDs of nodes along the wall
nozzleWallNodeIDs         = [];

%% Set up -- Initiliaze Variables %%

% COMPUTE FLOW VARIABLES FROM CHAMBER AND AMBIENT CONDITIONS HERE

% Compute the deflection of the flow across each expansion wave
% (characteristic line) at the throat. This is in degrees.
thetaMax = prandtlMeyer(Me,G)/2;
thetaMat = linspace(0.375,thetaMax,numCharLines);

%% Set Up -- Compute Initial Expansion Fan %%

% Compute initial expansion fan and assign CMinus lines
for i = 1:numCharLines
    % Compute Theta, Nu, Mach, Mu across each expansion wave
    Theta = thetaMat(i);
    Nu    = Theta;
    Mach  = inversePrandtlMeyer(Nu);
    Mu    = asind(1/Mach);

    % Compute the KMinus values
    KMinus = Theta + Nu;

    % Assign nodes to CMinus lines
    Node{i}.ID     = i;
    Node{i}.Theta  = Theta;
    Node{i}.Nu     = Nu;
    Node{i}.Mach   = Mach;
    Node{i}.Mu     = Mu;
    Node{i}.KMinus = KMinus;

    % Compute state variables and assign too nodes here
%     [R, T, P]     = calculateStateFromStag(Mach,G,Po,To,Ro);
%     UMag          = Mach*sqrt(G*R*T);
%     U             = UMag*cosd(Theta);
%     V             = UMag*sind(Theta);
%     Node{i}.R = R;
%     Node{i}.U = U;
%     Node{i}.V = V;
%     Node{i}.T = T;
%     Node{i}.P = P;
    
    
    % Compute coordinates and KPlus values for each node
    if i == 1
         Node{i}.KPlus  = Theta - Nu;
         Node{i}.Coords = computeIntersectionPoint([0,AStar/2], ...
                                                         Theta-Mu,    ...
                                                         [0, 0],      ...
                                                         0.0);
    else
         Node{i}.KPlus  = Node{i-1}.KPlus;
         Node{i}.Coords = computeIntersectionPoint([0,AStar/2], ...
                                                         Theta-Mu,    ...
                                                         Node{i-1}.Coords,      ...
                                                         Node{i-1}.Theta + Node{i-1}.Mu);

    end
end

% Compute the final "initial" node that is located on the nozzle wall
Node{numCharLines+1}.ID     = numCharLines + 1;
Node{numCharLines+1}.Coords = computeIntersectionPoint([0,AStar/2], ...
                                                         Theta,    ...
                                                         Node{numCharLines}.Coords,      ...
                                                         Node{numCharLines}.Theta + Node{i-1}.Mu);
Node{numCharLines+1}.Theta  = Theta;
Node{numCharLines+1}.Nu     = Nu;
Node{numCharLines+1}.Mach   = Mach;
Node{numCharLines+1}.Mu     = Mu;
Node{numCharLines+1}.KPlus  = Node{numCharLines}.KPlus;
Node{numCharLines+1}.KMinus = Theta + Nu;
% Node{numCharLines+1}.R      = R;
% Node{numCharLines+1}.U      = U;
% Node{numCharLines+1}.V      = V;
% Node{numCharLines+1}.T      = T;
% Node{numCharLines+1}.P      = P;

% Add this nodeID to the list of nodes that make up the nozzle wall
nozzleWallNodeIDs = [nozzleWallNodeIDs, numCharLines+1];

% Decrease the number of nodes per line by 1 every time a complete
% characteristic is computed
nodesPerCharLine = nodesPerCharLine - 1;

% Indicate the ID of the node at the centerline
centerlineNodeID = numCharLines + 2;

%% Main -- Compute Remaining Characteristic Nodes %%

% Create a switch to continue running this loop or not
keepComputing = true;

while keepComputing
    for i = centerlineNodeID:(centerlineNodeID + nodesPerCharLine-1)
       if i == centerlineNodeID
           % Compute nodal properties
           Theta  = 0.0;                              % Axially aligned along the centerline
           KMinus = Node{i-nodesPerCharLine}.KMinus;  % Inherit KMinus from internal node
           Nu     = KMinus - Theta;
           KPlus  = Theta - Nu;
           Mach   = inversePrandtlMeyer(Nu);
           Mu     = asind(1/Mach);
%            [R, T, P]     = calculateStateFromStag(Mach,G,Po,To,Ro);
%            UMag          = Mach*sqrt(G*R*T);
%            U             = UMag*cosd(Theta);
%            V             = UMag*sind(Theta);
           point1 = Node{i-nodesPerCharLine}.Coords;
           angle1 = Node{i-nodesPerCharLine}.Theta - Node{i-nodesPerCharLine}.Mu;
           point2 = [0, 0];
           angle2 = 0.0;

        % Assign Nodal Properties
        Node{i}.ID     = i;
        Node{i}.Theta  = Theta;
        Node{i}.Nu     = Nu;
        Node{i}.Mach   = Mach;
        Node{i}.Mu     = Mu;
        Node{i}.KPlus  = KPlus;
        Node{i}.KMinus = KMinus;
        Node{i}.Coords = computeIntersectionPoint(point1, angle1, point2, angle2);
%         Node{i}.R = R;
%         Node{i}.U = U;
%         Node{i}.V = V;
%         Node{i}.T = T;
%         Node{i}.P = P;

        
       elseif i == centerlineNodeID + nodesPerCharLine - 1
           % Mode nodal properties are inherited, compute intersection with
           % wall
           point1 = Node{i-nodesPerCharLine}.Coords;
           angle1 = Node{i-nodesPerCharLine}.Theta;
           point2 = Node{i-1}.Coords;
           angle2 = Node{i-1}.Theta + Node{i-1}.Mu;

        % Assign Nodal Properties
        Node{i}.ID     = i;
        Node{i}.Theta  = Theta;
        Node{i}.Nu     = Nu;
        Node{i}.Mach   = Mach;
        Node{i}.Mu     = Mu;
        Node{i}.KPlus  = KPlus;
        Node{i}.KMinus = Theta + Nu;
        Node{i}.Coords = computeIntersectionPoint(point1, angle1, point2, angle2);
%         Node{i}.R = R;
%         Node{i}.U = U;
%         Node{i}.V = V;
%         Node{i}.T = T;
%         Node{i}.P = P;

       % Add this nodeID to the list of nodes that make up the nozzle wall
       nozzleWallNodeIDs = [nozzleWallNodeIDs, i];
       else
           % Compute nodal properties
           KPlus  = Node{i-1}.KPlus;                  % Inherit KPlus from internal node
           KMinus = Node{i-nodesPerCharLine}.KMinus;  % Inherit KMinus from internal node
           Nu     = 0.5*(KMinus - KPlus);
           Theta  = 0.5*(KMinus + KPlus);
           Mach   = inversePrandtlMeyer(Nu);
           Mu     = asind(1/Mach);
%            [R, T, P]     = calculateStateFromStag(Mach,G,Po,To,Ro);
%            UMag          = Mach*sqrt(G*R*T);
%            U             = UMag*cosd(Theta);
%            V             = UMag*sind(Theta);
           point1 = Node{i-nodesPerCharLine}.Coords;
           angle1 = Node{i-nodesPerCharLine}.Theta - Node{i-nodesPerCharLine}.Mu;
           point2 = Node{i-1}.Coords;
           angle2 = Node{i-1}.Theta + Node{i-1}.Mu;

        % Assign Nodal Properties
        Node{i}.ID     = i;
        Node{i}.Theta  = Theta;
        Node{i}.Nu     = Nu;
        Node{i}.Mach   = Mach;
        Node{i}.Mu     = Mu;
        Node{i}.KPlus  = KPlus;
        Node{i}.KMinus = KMinus;
        Node{i}.Coords = computeIntersectionPoint(point1, angle1, point2, angle2);
%         Node{i}.R = R;
%         Node{i}.U = U;
%         Node{i}.V = V;
%         Node{i}.T = T;
%         Node{i}.P = P;

       end
    end

    % Update the centerline node ID and the nodes per line
    centerlineNodeID = i + 1;
    nodesPerCharLine = nodesPerCharLine - 1;

    % Check if the loop should terminate
    if nodesPerCharLine == 1
        keepComputing = false;
    end
end

%% Post Processing -- Group Nodes Into CMinus and CPlus Lines %%

% Group the CPlus lines (left running characteristics) first
centerlineNodeID = 1;
for i = 1:numCharLines
    for j = centerlineNodeID:nozzleWallNodeIDs(i)
        % Extract the proper node, make sure indexing begins at 1
        CPlus{i}{(j-centerlineNodeID)+1} = Node{j};
    end

    % Update centerline nodeID
    centerlineNodeID = nozzleWallNodeIDs(i) + 1;
end

% Group the CMinus lines (right running characteristics) second.
% Indexing gets tricky for CMinus lines, must continously change the number
% of nodes that are "skipped" to get to the next CMinus node
nodeStaggerMat = numCharLines:-1:0;
for i = 1:numCharLines
    % "starting" node of characteristic line
    nodeID = i;
    for j = 1:i
        CMinus{i}{j} = Node{nodeID};

        % Update the nodeID
        nodeID = nodeID + nodeStaggerMat(j);
    end
end

% Group the Nozzle wall nodes together
for i = 1:length(nozzleWallNodeIDs)
    NozzleNodes{i} = Node{nozzleWallNodeIDs(i)};
end

%% Post Processing -- Plot Characteristics and Nozzle Wall %%

% Close all figures
close all;

% Open and format a figure
Fig = figure;
set(Fig,'color','w');
Fig.Position = figPos;
daspect([1 1 1]);
hold on;

% Plot nodes
for i = 1:length(Node)
    coordsToPlot = Node{i}.Coords;
    plot(coordsToPlot(1),coordsToPlot(2),charNodeStyle,'MarkerSize',markerSize);
    drawnow;
end
clear coordsToPlot

% Plot CMinus Lines
for i = 1:length(CMinus)
    for j = 1:length(CMinus{i})
        coordsToPlot(j,:) = CMinus{i}{j}.Coords;
    end
    coordsToPlot = [0, AStar/2; coordsToPlot];
    plot(coordsToPlot(:,1),coordsToPlot(:,2),charLineStyle,'LineWidth',lineWidth-1);
    drawnow
    clear coordsToPlot
end

% Plot CPlus Lines
for i = 1:length(CPlus)
    for j = 1:length(CPlus{i})
        coordsToPlot(j,:) = CPlus{i}{j}.Coords;
        
    end
    plot(coordsToPlot(:,1),coordsToPlot(:,2),charLineStyle,'LineWidth',lineWidth-1);
    drawnow
    clear coordsToPlot
end

% Plot Nozzle Wall
coordsToPlot = [0, AStar/2];
for i = 1:length(NozzleNodes)
    coordsToPlot = [coordsToPlot; NozzleNodes{i}.Coords];
    plot(coordsToPlot(:,1),coordsToPlot(:,2),nozzleLineStyle,'LineWidth',lineWidth);
    drawnow
    clear coordsToPlot
    coordsToPlot = NozzleNodes{i}.Coords;
end
