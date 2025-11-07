close all
clear
clc

% Sompolinsky chaotic rate network integrated with ode45
%   dx_i/dt = -x_i + sum_j J_ij * tanh(x_j)
%   J_ij ~ N(0, g^2 / n), choose g slightly > 1 for weak chaos

rng(1); % reproducibility

% Parameters
numNeurons = 1000;      % n
g = 1.2;              % gain (slightly chaotic when > 1)
tspan = [0, 400];       % seconds

% Connectivity: Gaussian with mean 0 and variance g^2/n
J = (g / sqrt(numNeurons)) * randn(numNeurons, numNeurons);

% Initial condition (small random perturbation)
x0 = 0.1 * randn(numNeurons, 1);

% ODE right-hand side
rhs = @(t, x) -x + J * tanh(x);

% Integrate
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.05);
[t, X] = ode45(rhs, tspan, x0, opts);

% Plot: activity over time (all neurons)
figure('Color', 'w');
imagesc(t, 1:numNeurons, X.');
axis xy tight
colormap turbo
colorbar
xlabel('Time (s)')
ylabel('Neuron index')
title(sprintf('Sompolinsky network dynamics (n = %d, g = %.2f, phi = tanh)', numNeurons, g))

% Also show trajectories of a subset of neurons
figure('Color', 'w');
numToShow = min(10, numNeurons);
plot(t, X(:, 1:numToShow))
xlabel('Time (s)')
ylabel('x_i')
grid on
title(sprintf('Sample neuron trajectories (first %d)', numToShow))

