Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201;
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

%Nd = 1e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3

start_value = 1e16; % Define the start value
end_value = 20e16; % Define the end value

num_elements = nx; % Define the number of elements

% Generate the exponential gradient
% Convert start and end values to log10 scale
log10_start = log10(start_value);
log10_end = log10(end_value);

% Generate logarithmically spaced vector
Nd = logspace(log10_start, log10_end, num_elements);


NetDoping = ones(1,nx).*Nd; % doping

x0 = l/2;
nw = l/20;
%npDisturbance = 1e16*1e6*exp(-((x-x0)/nw).^2);
npDisturbance = 0;

LVbc = 0;
RVbc = 0;

TStop = 14200000*1e-18;
PlDelt = 100000*1e-18;

% PlotYAxis = {[-1e-15 2e-15] [-2e-9 2e-9] [-1.5e-12 1.5e-12]...
%     [1e22 2e22] [0 1e22] [0 20e43]...
%     [-20e33 15e33] [-2.5e34 2e34] [-1.1e8 1.1e8] ...
%     [-1e8 1e8] [-10e-3 10e-3] [0 2e22]};

doPlotImage = 0;
PlotFile = 'Gau2CarRC.gif';
