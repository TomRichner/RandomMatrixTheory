close all
clear
clc

% rng(8)

fs = 1000;
dt = 1/fs;
t = dt:dt:100;
nt = length(t);
tau_d = 0.5;

n = 16;
J = 0.5*randn(n,n);
z_mask = randn(n,n)<0.5;
J(z_mask) = 0;
J(:,1:n/2) = abs(J(:,1:n/2));
J(:,n/2+1:end) = -abs(J(:,n/2+1:end));

% Make each row sum to zero by adjusting non-zero entries
for i = 1:n
    row_sum = sum(J(i,:));
    non_zero_mask = ~z_mask(i,:);  % Find non-zero entries in this row
    num_non_zero = sum(non_zero_mask);
    if num_non_zero > 0
        % Distribute the correction across non-zero entries
        J(i, non_zero_mask) = J(i, non_zero_mask) - row_sum / num_non_zero;
    end
end


eigJ = 0.5*eig(J./tau_d);
orig_absc = 1.05*max(real(eigJ))
absc = orig_absc;

Input_1 = 4*randn(n,1).*ones(n,0.75*nt);
Input_1(rand(n,1)<0.7) = 0;

Input_2 = 4*randn(n,1).*ones(n,0.25*nt);
Input_2(rand(n,1)<0.7) = 0;

Input = 1*[Input_1, Input_2];

x = zeros(n,nt);
r = zeros(n,nt);

x(:,1) = max(0.25*abs(randn(n,1)),0);

absc_vec = zeros(1,nt);

for i_t = 2:nt

    r(:,i_t) = max(x(:,i_t-1),0);
    x(:,i_t) = x(:,i_t-1) + dt.*(-absc.*eye(n)*x(:,i_t-1) + J*r(:,i_t) + Input(:,i_t))./tau_d;
    absc_vec(1,i_t) = absc;

    if  nt/2 < i_t && i_t <0.75*nt % adjusting abscissa toward zero
        x_fp_1 = x(:,i_t);
        mask_1 = double(0 <= x_fp_1');
        A_fp_1 = (-absc.*eye(n) + J .* mask_1 + diag(diag(J)) )./tau_d;
        eig_A_fp_1 = eig(A_fp_1);
        absc = absc + 1*dt*tau_d*max(real(eig(A_fp_1))); % adjusting abscissa with negative feedback
        A_fp_1_eos = ( -absc.*eye(n) + J .* mask_1 + diag(diag(J)) )./tau_d;
        eig_A_fp_1_eos = eig(A_fp_1_eos);
    end

    if i_t == 0.75*nt
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
plot(t',x','LineWidth', 3)
hold on
subplot(2,1,1)
plot(t',r', 'LineWidth', 1)
hold off
subplot(2,1,2)
plot(t',absc')


figure(2)
s(1) = subplot(1,3,1);
plot(real(eig_J_AllActive),imag(eig_J_AllActive),'ok')
s(2) = subplot(1,3,2);
% plot(real(eig_A_fp_1),imag(eig_A_fp_1),'or')
% hold on
% plot(real(eig_A_fp_1_eos),imag(eig_A_fp_1_eos),'*b')
% hold on
plot(real(eig_A_fp_3),imag(eig_A_fp_3),'og')
hold off
s(3) = subplot(1,3,3);
plot(real(eig_A_fp_2),imag(eig_A_fp_2),'ok')
linkaxes(s,'xy')

axis equal

