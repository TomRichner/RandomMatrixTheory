close all
clear
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
n = 500;
b = 1;
mu = 0;
f = 0.5;
mu_E = 0.5; % relative to b
mu_I = -f.*mu_E./(1-f);

mean_indegree = 100; % desired mean number of connections in and out
density = mean_indegree/(n);
dense_mask = rand(n,n) < density; % mask of kept connections

%% E and I logical vectors
E = zeros(n,1,'logical');
E(1:round(f*n),1) = true; 
I = not(E); 

%% Ginibre matrix
A = b*randn(n,n)+mu;

%% apply sparse mask
A(not(dense_mask)) = 0; 

%% compute row-zeroing adjustment
[~, Z_row_adj] = row_sum_to_zero(A, dense_mask);

%% Rajan matrix to shift the mean of E and I neurons
M = zeros(n,n);
M(:,E) = mu_E;
M(:,I) = mu_I;

%% shift along main diagonal
W = -b*sqrt(n)*sqrt(density).*eye(n,n);

%% Resulting matrix
% Ginibre A
A_1_1 = A;
r_1_1 = b*sqrt(n)*sqrt(density);

% Rajan
A_1_2 = A+M;
r_1_2 = b*sqrt(n)*sqrt(density);

% Ginibre with mu
A_2_1 = A+0.06;
r_2_1 = b*sqrt(n)*sqrt(density);

% Ginibre with shift to edge
A_2_2 = A+W;
r_2_2 = b*sqrt(n)*sqrt(density);
x_2_2 = -b*sqrt(n)*sqrt(density);

% Rajan with row zero sum
A_3_1 = A+M+Z_row_adj;
r_3_1 = b*sqrt(n)*sqrt(density);


%% Make figures of different scenarios

f1 = figure(1);
set(f1, 'Position', [-1715        -114         640        1060])
tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


ax(1) = nexttile;
draw_radius(r_1_1,0)
hold on
plot_eig_distribution(ax(1), A_1_1)
hold off

ax(2) = nexttile;
plot_eig_distribution(ax(2), A_1_2)
hold on
draw_radius(r_1_2,0)
hold off

ax(3) = nexttile;
plot_eig_distribution(ax(3), A_2_1)
hold on
draw_radius(r_2_1,0)
hold off

ax(4) = nexttile;
plot_eig_distribution(ax(4), A_2_2)
hold on
draw_radius(r_2_2,x_2_2)
hold off

ax(5) = nexttile;
plot_eig_distribution(ax(5), A_3_1)
hold on
draw_radius(r_3_1,0)
hold off

linkaxes(ax,'y')

function plot_eig_distribution(target_ax, B)

    B_eigs = eig(B);
    scatter(target_ax, real(B_eigs),imag(B_eigs), 5, 'MarkerEdgeColor',[0 0 0 ]);
    
    axis(target_ax, 'equal');
    xlabel(target_ax, 'Re($\lambda$)', 'Interpreter', 'latex');
    ylabel(target_ax, 'Im($\lambda$)', 'Interpreter', 'latex');
    
end

function draw_radius(r, x)
    % y_center is implicitly 0
    pos = [x - r, -r, 2*r, 2*r];
    rectangle('Position', pos, 'Curvature', [1,1], 'EdgeColor', 'k', 'LineWidth', 1);
end

function [Z_out, Z_adj] = row_sum_to_zero(Z_in, dense_mask)
    Z_out = Z_in;
    [n, ~] = size(Z_in);

    for i = 1:n
        row = Z_out(i, :);
        current_sum = sum(row);

        % Find non-zero elements in the row corresponding to the dense_mask
        adjustment_mask = dense_mask(i, :);
        num_to_adjust = sum(adjustment_mask);

        if num_to_adjust > 0
            adjustment = -current_sum / num_to_adjust;
            row(adjustment_mask) = row(adjustment_mask) + adjustment;
            Z_out(i, :) = row;
        end
    end
    Z_adj = Z_out - Z_in;
end