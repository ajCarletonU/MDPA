%vertical stream of particles hitting another bunch of particles
%----------------------------------------------------------------


%PA2 - 20 Jan 2022
% Change to code 
% - added the call to AddEllipticalAtomicArray function 
% - configured the dimension to fit an elliptical shape
% - added two streams of particles at different locations (x-axis)
% - changed the number of particles that will interact with the atoms


doPlot = 1;
dt = 5e-15;
TStop = 3000 * dt;          %stop time
InitDist = 0.0;             %iniital disturbance any particles present\

%method of integration used
%--------------------------
Method = 'VE'; % VE -- verlot; FD -- Forward Difference


%sets the mass of two types of particles
%----------------------------
Mass0 = 14 * C.am; % Silicon
Mass1 = 100 * C.am; % Argon


%sets the atomic spacing size and parameters of the LJ potential
%----------------------------
AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing / 2^(1 / 6);
LJEpsilon = 1e-21;


PhiCutoff = 3 * AtomSpacing * 1.1;      %set the cutoff at which we stop calc forces

T = 30;         %set the temperatur of the system


% AddRectAtomicArray(10, 10, 0, 0, 0, 0, 0, T, 0);              % 10x10  at some T
% % vy0 = -sqrt(0.02*Ep/Mass1);
% % AddRectAtomicArray(4,4,0,12*AtomSpacing,0,vy0,0,T,1);
% AddParticleStream(5, 0.1, 10, -pi / 2, 1, Ep * C.q_0, 3);

%--------------------------------------------------------------
%calling Elliptical function and adding in an array of atoms
%parameter used will be length=11 and height=5  at some T  temp
%--------------------------------------------------------------
%(RADatomsA,RADatomsB,X0,Y0,VX0,VY0,InitDist,Temp,Type)
%the value of RADstomsA needs to be > RADatomsB
AddEllipticalAtomicArray(11, 5, 0, 0, 0, 0, 0, T, 0);   
Ep = 2;

%add two streams of particles at different locations and spacings 
%(num, x0, y0, PartAng, Type, Ep, Seper)
AddParticleStream(3, 1, 8, -pi / 2, 1, Ep * C.q_0, 4);
AddParticleStream(5,-4, 8, -pi / 2, 1, Ep * C.q_0, 2);


Size = 10*AtomSpacing;
Limits = [-Size +Size -Size +Size]; % square is good
PlDelt = 5 * dt;

PlotFile = 'BlockSt.gif';
PlotPosOnly = 1;
doPlotImage = 0;
PlotSize = [100, 100, 1049, 1049];

ScaleV = .02e-11;
ScaleF = 10;
