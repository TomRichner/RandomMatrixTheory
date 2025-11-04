% Threshold-linear (biased ReLU with saturation) plot
% phi(x) = clip(x + gamma, 0, phi_max)

clear; clc;

% --- Parameters (edit as desired) ---
gamma   = -1;      % bias (shifts the knee to x = -gamma)
phi_max = 2.0;      % saturation level

% --- Domain and definition ---
x  = linspace(-3, 3, 2001);
y  = min(max(x + gamma, 0), phi_max);   % clip(x+gamma, 0, phi_max)

% Breakpoints
x0 = -gamma;                % lower knee: where x+gamma crosses 0
x1 = phi_max - gamma;       % upper knee: where x+gamma reaches phi_max

% --- Plot ---
figure('Color','w'); hold on;
plot(x, y, 'LineWidth', 2);

% Guides at the knees and satu
