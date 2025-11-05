close all
clear
clc

n = 200;
J = randn(n,n);
z_mask = randn(n,n)<0.5;
J(z_mask) = 0;
eigJ = eig(J);
absc = max(real(eigJ));

% Jshifted = J;


fs = 1000;
dt = 1/fs;
t = dt:dt:20;
nt = length(t)

Input_1 = randn(n,1).*ones(n,nt/2);
% Input_1(Input_1<0) = 0; % sparse positive input

Input_2 = randn(n,1).*ones(n,nt/2);
% Input_2(Input_2<0) = 0;

Input = 1*[Input_1, Input_2];

x = zeros(n,nt);
r = zeros(n,nt);

x(:,1) = max(0.25*abs(randn(n,1)),0);

tau_d = 0.5;

for i_t = 2:nt

    r(:,i_t) = max(x(:,i_t-1),0);
    x(:,i_t) = x(:,i_t-1) + dt.*(-absc.*eye(n)*x(:,i_t-1) + J*r(:,i_t) + Input(:,i_t))./tau_d;

    if i_t == nt/2
        x_fixedPoint_1 = x(:,i_t);
    end

    if i_t == nt
        x_fixedPoint_2 = x(:,i_t);
    end

end

J_AllActive = ( -absc.*eye(n) + J )./tau_d;
eig_J_AllActive = eig(J_AllActive);

mask_1 = double(0 < x_fixedPoint_1');
J_fixedPoint_1 = ( -absc.*eye(n) + J .* mask_1 + diag(diag(J)) )./tau_d; % zero out off-diagonals of columns where x <= 0, keep diagonal
eig_J_x1 = eig(J_fixedPoint_1);
mask_2 = double(0 < x_fixedPoint_2');
J_fixedPoint_2 = ( -absc.*eye(n) + J .* mask_2 + diag(diag(J)) )./tau_d; % zero out off-diagonals of columns where x <= 0, keep diagonal
eig_J_x2 = eig(J_fixedPoint_2);

figure(1)
% plot(t',x','LineWidth', 3)
% hold on
plot(t',r', 'LineWidth', 1)
hold off


figure(2)
s(1) = subplot(1,3,1);
plot(real(eig_J_AllActive),imag(eig_J_AllActive),'ok')
s(2) = subplot(1,3,2);
plot(real(eig_J_x1),imag(eig_J_x1),'ok')
s(3) = subplot(1,3,3);
plot(real(eig_J_x2),imag(eig_J_x2),'ok')
linkaxes(s,'xy')

axis equal

