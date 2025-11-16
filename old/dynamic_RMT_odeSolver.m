close all
clear
clc

n = 10;
J = randn(n,n);
z_mask = randn(n,n)<0.5;
J(z_mask) = 0;
eigJ = eig(J);
absc = max(real(eigJ));
Jshifted = J-1.5*absc*eye(n,n);
eigJshifted = eig(Jshifted);


% Time and input setup
t_switch = 10;  % Time to switch inputs
t_end = 20;     % Final time

% Generate inputs (constant for each phase)
Input_1 = randn(n,1);
% Input_1(Input_1<0) = 0; % sparse positive input

Input_2 = randn(n,1);
% Input_2(Input_2<0) = 0;

% Initial condition
x0 = max(randn(n,1),0);

% Phase 1: Integrate from 0 to t_switch with Input_1
[t1, x1] = ode15s(@(t,x) odefun(t, x, Jshifted, Input_1), [0 t_switch], x0);

% Extract state at end of phase 1
x_fixedPoint_1 = x1(end,:)';

% Phase 2: Integrate from t_switch to t_end with Input_2
[t2, x2] = ode15s(@(t,x) odefun(t, x, Jshifted, Input_2), [t_switch t_end], x_fixedPoint_1);

% Extract state at end of phase 2
x_fixedPoint_2 = x2(end,:)';

% Combine results from both phases
t = [t1; t2(2:end)];  % Skip first point of phase 2 to avoid duplication
x = [x1; x2(2:end,:)];

% Jacobian analysis at fixed points
J_fixedPoint_1 = Jshifted.* double(x_fixedPoint_1'>0); % zero out columns of J based on fixed point
eig_J_x1 = eig(J_fixedPoint_1);
J_fixedPoint_2 = Jshifted.* double(x_fixedPoint_2'>0); % zero out columns of J based on fixed point
eig_J_x2 = eig(J_fixedPoint_2);

% Plotting
figure(1)
plot(t, x)
xlabel('Time (s)')
ylabel('State x')
title('State Evolution')

figure(2)
s(1) = subplot(1,3,1);
plot(real(eigJshifted),imag(eigJshifted),'ok')
title('Original Jshifted')
xlabel('Real')
ylabel('Imaginary')

s(2) = subplot(1,3,2);
plot(real(eig_J_x1),imag(eig_J_x1),'ok')
title('Jacobian at t=10s')
xlabel('Real')
ylabel('Imaginary')

s(3) = subplot(1,3,3);
plot(real(eig_J_x2),imag(eig_J_x2),'ok')
title('Jacobian at t=20s')
xlabel('Real')
ylabel('Imaginary')

linkaxes(s,'xy')
axis equal


%% Subfunction: ODE right-hand side
function dxdt = odefun(t, x, Jshifted, Input_vec)
    % Compute the derivative
    dxdt = Jshifted * x + Input_vec;
    
    % Enforce ReLU constraint: if x(i)<=0 and trying to go more negative, stop
    % This prevents x from becoming negative
    mask = (x <= 0) & (dxdt < 0);
    dxdt(mask) = 0;
end

