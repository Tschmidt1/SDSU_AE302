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
xMax  = 4.0;    % x location of nozzle exit

% Flow variables
Me = 2.4;       % Mach number at exit plane
G  = 1.4;       % Ratio of specific heats
R  = 287;       % Gas constant
Po = 707392;    % Chamber pressure     [Pascals]
To = 520;       % Chamber temperature  [Kelvin]
Ro = Po/R/To;   % Chamber density      [Kg/m^3]

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


%% Set up -- Initiliaze Variables %%

% Create x coordinates for nozzle plotting later in the script
xNoz = linspace(xMin,xMax,1000);

% Compute the initial guess of the shock location as halfway though the nozzle
xShock = (xMax-xMin)/2;

% Create a variable to switch shock computations off once the location has been
% determined
computeShock = true;

% Compute exit "area"
Ae = 2*calculateNozzleContourFromEquation(AStar,xMax);

%% Set Up -- Check if Shock May Form %%

% Compute Mach at the exit plane, assuming isentropic flow up to the exit plane


%% Main -- Compute Remaining Characteristic Nodes %%



%% Post Processing -- Group Nodes Into CMinus and CPlus Lines %%



%% Post Processing -- Plot Characteristics and Nozzle Wall %%

% Close all figures
close all;

