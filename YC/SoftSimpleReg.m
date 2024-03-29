winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt (c_mu_0/c_eps_0);


tSim = 200e-15 ;                   %Sets simulation time 
f = 230e12;                       %Sets frequency  
lambda = c_c/f;

xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;

epi{1} = ones(nx{1},ny{1})*c_eps_0;
epi{1}(60:80,55:95)= c_eps_0*11.3; %Line that adds the inclusion
%When commented out, the wave is unobstructed by a grating
epi{1}(20:40,55:95)= c_eps_0*11.3;
epi{1}(100:120,55:95)= c_eps_0*11.3;

epi{1}(60:80,15:45)= c_eps_0*11.3; 
epi{1}(20:40,15:45)= c_eps_0*11.3;
epi{1}(100:120,15:45)= c_eps_0*11.3;

epi{1}(60:80,105:135)= c_eps_0*11.3; 
epi{1}(20:40,105:135)= c_eps_0*11.3;
epi{1}(100:120,105:135)= c_eps_0*11.3;

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};                         %assigns numerical values to
dt = 0.25*dx/c_c;                           %dx and dt to make them diff-    
nSteps = round(tSim/dt*2);                  %erence eqns not differential
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx;

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;                       %Speed at which the simulation samples
Plot.MaxEz = 2;
Plot.MaxH = Plot.MaxEz/c_eta_0;     %Defines the scales of the plot
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];


bc{1}.NumS = 2;                     %Defines the number of plane wave sources
bc{1}.s(1).xpos = nx{1}/4+1;    %Specifies the starting position for the wave
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;      %Calles fct to create the plane wave for sim
% mag = -1/c_eta_0;                   

bc{1}.s(2).xpos = nx{1}/4+1;    %Specifies the starting position for the wave
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC; 

mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = 15e-15;
s = 0;
y0 = yMax*3/4;
y1 = yMax/4;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y1,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';  %Control the reflection or abosorbtion of the 
bc{1}.xp.type = 'e';  %edges of the simulation window from the left and right
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';



pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






