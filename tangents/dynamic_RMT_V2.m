close all
clear
clc

n = 16;
J = randn(n,n);
z_mask = randn(n,n)<0.5;
J(z_mask) = 0;
eigJ = eig(J);
orig_absc = 1.35*max(real(eigJ))
absc = orig_absc;

% Jshifted = J;


fs = 1000;
dt = 1/fs;
t = dt:dt:20;
nt = length(t)

Input_1 = randn(n,1).*ones(n,0.75*nt);
% Input_1(Input_1<0) = 0; % sparse positive input

Input_2 = randn(n,1).*ones(n,0.25*nt);
% Input_2(Input_2<0) = 0;

Input = 1*[Input_1, Input_2];

x = zeros(n,nt);
r = zeros(n,nt);

x(:,1) = max(0.25*abs(randn(n,1)),0);

tau_d = 0.5;

for i_t = 2:nt

    r(:,i_t) = max(x(:,i_t-1),0);
    x(:,i_t) = x(:,i_t-1) + dt.*(-absc.*eye(n)*x(:,i_t-1) + J*r(:,i_t) + Input(:,i_t))./tau_d;

    if i_t == nt/2 % correct negative feedback
        x_fp_1 = x(:,i_t);
        mask_1 = double(0 <= x_fp_1');
        A_fp_1 = (-absc.*eye(n) + J .* mask_1 + diag(diag(J)) )./tau_d;
        eig_A_fp_1 = eig(A_fp_1);
        % A_fp_1_unshift = ( 0 + J .* mask_1 + diag(diag(J)) )./tau_d; % undo shift get absc again, or just calc the absc_correction needed
        new_absc = absc + tau_d*max(real(eig(A_fp_1))) % not correct due to tau_d
        absc = 1.00*new_absc; % replace the absc for the stimulation
        A_fp_1_eos = ( -absc.*eye(n) + J .* mask_1 + diag(diag(J)) )./tau_d; % now on edge of stability with new 
        eig_A_fp_1_eos = eig(A_fp_1_eos);
    end

    if i_t == 0.75*nt-1

        mask_3 = double(0 <= r(:,i_t)');
        A_fp_3 = (-absc.*eye(n) + J .* mask_3 + diag(diag(J)) )./tau_d;
        eig_A_fp_3 = eig(A_fp_3);
        
    end

    if i_t == nt
        x_fp_2 = x(:,i_t);
        mask_2 = double(0 <= x_fp_2');
        A_fp_2 = ( -absc.*eye(n) + J .* mask_2 + diag(diag(J)) )./tau_d; 
        eig_A_fp_2 = eig(A_fp_2);
        end_absc = max(real(eig(A_fp_2)))
    end

end

J_AllActive = ( -absc.*eye(n) + J )./tau_d;
eig_J_AllActive = eig(J_AllActive);


figure(1)
% plot(t',x','LineWidth', 3)
% hold on
plot(t',r', 'LineWidth', 1)
hold off


figure(2)
s(1) = subplot(1,3,1);
plot(real(eig_J_AllActive),imag(eig_J_AllActive),'ok')
s(2) = subplot(1,3,2);
plot(real(eig_A_fp_1),imag(eig_A_fp_1),'or')
hold on
plot(real(eig_A_fp_1_eos),imag(eig_A_fp_1_eos),'*b')
hold on
plot(real(eig_A_fp_3),imag(eig_A_fp_3),'og')
hold off
s(3) = subplot(1,3,3);
plot(real(eig_A_fp_2),imag(eig_A_fp_2),'ok')
linkaxes(s,'xy')

axis equal

