close all
clear
clc

% Sompolinsky 1988 chaotic rate network integrated with ode45
%   dx_i/dt = -x_i + sum_j J_ij * tanh(x_j)
%   J_ij ~ N(0, g^2 / n), choose g slightly > 1 for weak chaos

rng(4); % reproducibility

% Parameters
n = 20;      % n
g = 2.9;              % gain (slightly chaotic when > 1)
fs = 1000;          % sampling frequency (Hz)
dt = 1/fs;
tspan = [0, 400];       % seconds
t = linspace(tspan(1), tspan(2), (tspan(2)-tspan(1))*fs + 1)'; % time vector
nt = length(t)
mu_E = .25;
mu_I = -0.25;
M = [mu_E.*ones(n,n/2), mu_I.*ones(n,n/2)];
% Connectivity: Gaussian with mean 0 and variance g^2/n
J = (g / sqrt(n)) * (randn(n, n)+M);

J=J-mean(J,1)./n; % row sum to zero

% External input stimulus (n x time)
u = zeros(n, nt); % User can customize this as needed

amp = 0.25;
u(:,1:fix(nt/4)) = amp.*repmat(randn(n,1), 1, fix(nt/4));
u(:,fix(nt/4)+(1:fix(nt/4))) = amp.*repmat(randn(n,1), 1, fix(nt/4));
u(:,fix(nt/2)+(1:fix(nt/4))) = amp.*repmat(randn(n,1), 1, fix(nt/4));
u(:,fix(0.75*nt)+(1:fix(nt/4))-1) = amp.*repmat(randn(n,1), 1, fix(nt/4));

% Initial condition (small random perturbation)
x0 = 0.1 * randn(n, 1);

% ODE right-hand side (now with external stimulus)
rhs = @(t_current, x) rhs_somp(t_current, x, t, u, J);

% Integrate
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', dt);
[t_out, X] = ode45(rhs, t, x0, opts);

%% Compute largest Lyapunov exponent using Benettin's algorithm
fprintf('Computing largest Lyapunov exponent...\n');
params.n = n;
params.J = J;
d0 = 1e-3;
lya_dt = 5*dt;
lya_fs = 1/lya_dt;
T = tspan;
ode_solver = @ode45;
[LLE, local_lya, finite_lya, t_lya] = benettin_algorithm(X, t_out, dt, fs, d0, T, lya_dt, params, opts, @rhs_somp_wrapper, t, u, ode_solver);
fprintf('Largest Lyapunov Exponent: %.4f\n', LLE);

%% Plots
% Plot: activity over time (all neurons)
figure('Color', 'w');
set(gcf,'Position',[-2302         726        1574         543])
imagesc(t_out, 1:n, X.');
axis xy tight
colormap turbo
colorbar
xlabel('Time (s)')
ylabel('Neuron index')
title(sprintf('Sompolinsky network dynamics (n = %d, g = %.2f)', n, g))

% Also show trajectories of a subset of neurons and Lyapunov exponent
figure('Color', 'w');
set(gcf,'Position',[-2244          82        1340         800])

% Subplot 1: Sample neuron trajectories
subplot(2,1,1)
numToShow = min(30, n);
plot(t_out, X(:, 1:numToShow))
xlabel('Time (s)')
ylabel('x_i')
grid on
title(sprintf('Sample neuron trajectories (first %d)', numToShow))

% Subplot 2: Local Lyapunov exponent

[bL,aL] = butter(3,0.05/(lya_fs/2),'low');
local_lya_filt = filtfilt(bL,aL,local_lya);

subplot(2,1,2)
plot(t_lya, local_lya, 'LineWidth', 1.5)
hold on
plot(t_lya, local_lya_filt,'LineWidth',1.5)
yline(LLE, '--r', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Local Lyapunov Exponent')
grid on
title(sprintf('Local Lyapunov Exponent (LLE = %.4f)', LLE))
legend('Local', 'Mean LLE', 'Location', 'best')

%% Subfunction: ODE right-hand side with external stimulus
function dx_dt = rhs_somp(t_current, x, t_vec, u, J)
    % Persistent variables for efficient interpolation
    persistent u_interpolant t_vec_last;
    
    % Create or update griddedInterpolant if time vector has changed
    if isempty(u_interpolant) || isempty(t_vec_last) || ...
       numel(t_vec_last) ~= numel(t_vec) || t_vec_last(1) ~= t_vec(1) || t_vec_last(end) ~= t_vec(end)
        
        % Use 'linear' interpolation and 'none' for extrapolation
        % This returns NaN for out-of-bounds queries to catch errors
        u_interpolant = griddedInterpolant(t_vec, u', 'linear', 'none');
        t_vec_last = t_vec;
    end
    
    % Interpolate u at current time (u_interpolant returns 1 x n, transpose to n x 1)
    u_t = u_interpolant(t_current)';
    
    % Compute derivative: dx/dt = -x + J * relu(x) + u(t)
    % relu with saturation at 5
    dx_dt = -x + J * min(max(x, 0), 5) + u_t;
end

%% Wrapper function for benettin_algorithm
function dx_dt = rhs_somp_wrapper(tt, XX, t_ex, u_ex, params)
    dx_dt = rhs_somp(tt, XX, t_ex, u_ex, params.J);
end
