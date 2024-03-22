% Diode parameters
Is = 0.01e-12; % A
Ib = 0.1e-12;  % B (this will be fixed in the fit function)
Vb = 1.3;      % D (this will be fixed in the fit function)
Gp = 0.1;      % C (this is to be fitted)
VT = 0.025; % Thermal voltage (V)

% Voltage vector from -1.95 to 0.7 volts with 200 steps
V = linspace(-1.95, 0.7, 200);

% Ideal diode current
I_ideal = Is * (exp(V / VT) - 1) + Gp * V - Ib * (exp((V + Vb) / VT) - 1);

% Introduce 20% random noise
noise_factor = 0.20;
I_noisy = I_ideal .* (1 + noise_factor * (rand(size(V)) - 0.5));

% Polynomial fitting
% p4 = polyfit(V, I_noisy, 4); % 4th order polynomial fit
% p8 = polyfit(V, I_noisy, 8); % 8th order polynomial fit
% 
% % Evaluating the polynomial fits
% I_p4 = polyval(p4, V);
% I_p8 = polyval(p8, V);

% Plotting the I-V characteristic with noise and polynomial fits
% figure;
% subplot(1, 2, 1);
% plot(V, I_noisy, 'r.', V, I_p4, 'b-', V, I_p8, 'g-');
% title('Diode I-V Characteristics (Linear Scale)');
% xlabel('Voltage (V)');
% ylabel('Current (A)');
% legend('Noisy data', '4th order fit', '8th order fit');
% grid on;
% 
% subplot(1, 2, 2);
% semilogy(V, abs(I_noisy), 'r.', V, abs(I_p4), 'b-', V, abs(I_p8), 'g-'); % Using absolute value to handle negative currents in log scale
% title('Diode I-V Characteristics (Semi-log Scale)');
% xlabel('Voltage (V)');
% ylabel('Current (A)');
% legend('Noisy data', '4th order fit', '8th order fit');
% grid on;

% % Nonlinear curve fitting using fit function PART A
% fitType = fittype('A.*(exp(1.2*x/25e-3) - 1) + C.*x - B.*(exp(1.2*(x+D)/25e-3) - 1)', ...
%     'independent', 'x', ...
%     'dependent', 'y', ...
%     'coefficients', {'A', 'C'}, ...
%     'problem', {'B', 'D'});
% 
% % Start points for A and C
% startPoints = [Is, Gp];
% 
% % Using the values from equation 1 for B and D
% [B, D] = deal(Ib, Vb);
% 
% % Perform the fitting
% [fitResult, gof] = fit(V', I_noisy', fitType, 'problem', {B, D}, 'StartPoint', startPoints);
% 
% % Generate the fitted curve
% I_fit = feval(fitResult, V);

% % Nonlinear curve fitting using fit function PART B
% fitType = fittype('A.*(exp(1.2*x/25e-3) - 1) + C.*x - B.*(exp(1.2*(x+D)/25e-3) - 1)', ...
%     'independent', 'x', ...
%     'dependent', 'y', ...
%     'coefficients', {'A', 'B', 'C'}, ...
%     'problem', {'D'});
% 
% % Start points for A, B, and C
% startPoints = [Is, Ib, Gp];
% 
% % Perform the fitting
% [fitResult, gof] = fit(V', I_noisy', fitType, 'problem', D, 'StartPoint', startPoints);
% 
% % Generate the fitted curve
% I_fit_b = feval(fitResult, V);

% % Nonlinear curve fitting using fit function Part C
% fitType = fittype('A.*(exp(1.2*x/25e-3) - 1) + C.*x - B.*(exp(1.2*(x+D)/25e-3) - 1)', ...
%     'independent', 'x', ...
%     'dependent', 'y', ...
%     'coefficients', {'A', 'B', 'C', 'D'});
% 
% % Start points for A, B, C, and D
% startPoints = [Is, Ib, Gp, Vb];
% 
% % Perform the fitting
% [fitResult, gof] = fit(V', I_noisy', fitType, 'StartPoint', startPoints);
% 
% % Generate the fitted curve
% I_fit_all = feval(fitResult, V);

inputs = V.';
targets = I_noisy.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
view(net);
Inn = outputs;

% Plotting the I-V characteristic with noisy data and fitted curve
figure;

% Plot with noisy data and fitted curve using fit function
subplot(1, 2, 1);
plot(V, I_noisy, 'r.', V, Inn, 'b-');
title('Diode I-V Characteristics with Fit (Linear Scale)');
xlabel('Voltage (V)');
ylabel('Current (A)');
ylim([-6e21,0]);
legend('Noisy data', 'Fitted curve');
grid on;

% Fit on semi-log scale
subplot(1, 2, 2);
semilogy(V, abs(I_noisy), 'r.', V, abs(Inn), 'b-');
title('Diode I-V Characteristics with Fit (Semi-log Scale)');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Noisy data', 'Fitted curve');
grid on;

