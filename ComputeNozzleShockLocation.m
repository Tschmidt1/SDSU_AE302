%{ 
This script computes the x coordinate of a shock that appears in a supersonic
nozzle.

This is a 1D computation, so variables such as "AStar" are not
true area's, they are unit areas (a.k.a. lengths). 

The coordinate system is set such that (0,0) is located at the nozzle
throat, along the nozzle centerline. 

The necessary input variables are:
    AStar   -- "area" at the throat. For convienence, this is 1.0 by default
    Po      -- the chamber pressure in Pascals
    To      -- the chamber temperature in degrees Kelvin
    Pa      -- the ambient pressure
%}

%% Basic Setup %%

close all;
clear all;
clc;

% Add path to Matlab that are called in this script here
addpath HelperFunctions

% Set default interpreters
set(0,'defaulttextinterpreter','latex')
set(0,'defaulttextfontname','times')
set(0,'defaultaxesfontname','times')

warning('OFF','ALL');

%% Basic Variables %%

% Geometric variables
AStar = 1.0;    % "area" of the throat
xMin  = 0.0;    % x location of throut
xMax  = 144.0;  % x location of nozzle exit

% Flow variables
Me = 2.4;       % Mach number at exit plane
G  = 1.4;       % Ratio of specific heats
R  = 287;       % Gas constant
Po = 1500000;   % Chamber pressure     [Pascals]
To = 3000;      % Chamber temperature  [Kelvin]
Ro = Po/R/To;   % Chamber density      [Kg/m^3]
Pa = 101325;    % Ambient pressure     [Pascals] 


% Plotting variables
fontSize        = 18;
fontName        = 'times';
lineWidth       = 2;

% Figure setup variables
figPos = [680 338 1030 640];

% Miscellaneous variables
tol = 1e-6;   % Error tolerance used for iteration

%% Set Up -- Create Data Structures %% 


%% Set up -- Initiliaze Variables %%

% Create x and y coordinates for nozzle plotting later in the script
xNoz = linspace(xMin,xMax,1000);
yNoz = calculateNozzleContourFromEquation(AStar,xNoz);

% Compute the initial guess of the shock location as halfway though the nozzle
xShock = (xMax-xMin)/2;
yShock = calculateNozzleContourFromEquation(AStar,xShock);

%  Create lower and upper bounds of shock iteration scheme
xShockMin = xMin;
xShockMax = xMax;

% Create a variable to switch shock computations off once the location has been
% determined
computeShock = true;

% Compute exit "area"
Ae = 2*calculateNozzleContourFromEquation(AStar,xMax);

% Create the figure for plotting the nozzle and the shock.
Fig = figure;
set(Fig,'color','w');
Fig.Position = figPos;
hold on;
ax = gca;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$y$';
ax.FontSize = fontSize;
ax.FontName = fontName;

% Plot the nozzle contour and initial shock guess
plot(xNoz,yNoz,'k-','LineWidth',lineWidth);
plot(xNoz,-yNoz,'k-','LineWidth',lineWidth);
s = plot([xShock, xShock],[-yShock, yShock],'r-','LineWidth',lineWidth);
drawnow;
pause(0.005);


%% Set Up -- Check if Shock May Form %%

% Compute Mach at the exit plane, assuming SUPERSONIC isentropic flow up 
% to the exit plane
Me = MachFromAOverAStar(1,Ae/AStar,G);

% Compute corresponding state variables
[Pe,Te,Re] = calculateStateFromStag(Me,G,Po,To,Ro);

% If no shock will form, hurrah, no need to run this!
if Pe > Pa
    fprintf('Pe = %7.0f \t Pa = %7.0f\n',Pe,Pa);
    fprintf('Pe > Pa: this flow is under-exanded so no shock forms inside the nozzle!\n')

    computeShock = false;
end

%% Main -- Compute Remaining Characteristic Nodes %%
iterationCount = 0;
while computeShock
    % Increase iteration counter
    iterationCount = iterationCount + 1;
    
    % Compute area at guessed shock location
    A = 2*calculateNozzleContourFromEquation(AStar,xShock);

    % Compute Mach number at guessed shock location
    M1 = MachFromAOverAStar(1,A/AStar,G);

    % Compute upstream variables at guess shock location
    [P1,T1,R1] = calculateStateFromStag(M1,G,Po,To,Ro);

    % Compute Mach and state variables downstrean of normal shock
    [M2, R2, P2, T2] = calculateNormalShock_Downstream(M1, R1, P1, T1, G, R);

    % Compute stagnation properties across the shock
    [Po2, To2, Ro2] = calculateStagFromState(M2,G,P2,T2,R2);

    % Compute the "new" AStar. Remember, this is an imaginary state we just
    % use for reference!
    AOverAStar2  = AOverAStarFromMach(M2,G);
    AeOverAStar2 = (Ae/AStar)*(AStar/A)*(AOverAStar2);

    % Compute Mach at exit area, assuming the SUBSONIC solution
    Me = MachFromAOverAStar(0,AeOverAStar2,G);

    % Compute the state variables at the exit plane using the subsonic Mach
    % number and downstream stagnation values
    [Pe,Te,Re] = calculateStateFromStag(Me,G,Po2,To2,Ro2);


    % Check if the shock was "too strong" or "too weak". In the exact
    % solution Pe = Pa. Downstream of the normal shock the flow is
    % subsonic, and in subsonic flow no discontinuities can exist. Thus, if
    % Pe > Pa the shock was "too strong"; the pressure jump across the
    % shock was too large, which means the shock was assumed too far
    % downstream in the nozzle. Conversely, if Pe < Pa the shock was "too
    % weak", meaning the pressure jump across the shock was too small and
    % the shock was assumed to far upstream in the nozzle. 

    % Print out results:
    fprintf('Iteration: %i \t Guessed shock location: %.7f \t Pe/Pa: %.7f\n',...
             iterationCount, xShock, Pe/Pa);

    % Plot the guessed shock
    yShock = calculateNozzleContourFromEquation(AStar,xShock);
    delete(s);
    s = plot([xShock, xShock],[-yShock, yShock],'r-','LineWidth',lineWidth);
    drawnow;
    pause(0.005);

    if abs(Pe - Pa) <= tol
        break;
        computeShock = false;
    elseif (Pe < Pa)  % Shock was too strong, move it upstream
        xShockNew = mean([xShockMin,xShock]);
        xShockMax = xShock;
    else              % Shock was too weak, move it downstream
        xShockNew = mean([xShockMax,xShock]);
        xShockMin = xShock;
    end 

    % Check if the loop if xShock is not changing. If this is the case, the
    % shock does not occur in the nozzle. It's an oblique shock outside of
    % the nozzle
    if abs(xShock - xShockNew) <= 1e-12
        fprintf('No normal shock develops in the nozzle!\n');
        break
    end
    
    % update xShock guess
    xShock = xShockNew;
end


