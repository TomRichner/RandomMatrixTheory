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
mean_indegree = 100; % desired mean number of connections in and out

%% Struct to hold RMT objects
G = struct();

%% Base Ginibre matrix
G(1).rmt = RMT(n, b, mu);
G(1).rmt.apply_sparsity(mean_indegree);
G(1).rmt.description = 'Ginibre';

%% Rajan
G(2).rmt = G(1).rmt.copy();
G(2).rmt.set_rajan_means(f, mu_E);
G(2).rmt.description = 'Rajan';

%% Ginibre with mu
G(3).rmt = G(1).rmt.copy();
G(3).rmt.add_constant(0.06);
G(3).rmt.description = 'Ginibre with mu';

%% Ginibre with shift to edge
G(4).rmt = G(1).rmt.copy();
shift = -G(4).rmt.b*sqrt(G(4).rmt.n)*sqrt(G(4).rmt.density);
G(4).rmt.shift_diagonal(shift);
G(4).rmt.description = 'Ginibre with shift';

%% Rajan with row zero sum
G(5).rmt = G(2).rmt.copy();
G(5).rmt.row_sum_to_zero();
G(5).rmt.description = 'Rajan row-zero';

%% Compute eigenvalues and set circles for all
for i=1:length(G)
    G(i).rmt.compute_eigenvalues();
    r = G(i).rmt.b*sqrt(G(i).rmt.n)*sqrt(G(i).rmt.density);
    x_center = 0;
    if i == 4 % Ginibre with shift
        x_center = -G(i).rmt.b*sqrt(G(i).rmt.n)*sqrt(G(i).rmt.density);
    end
    G(i).rmt.set_plot_circle(r, x_center);
end

%% Make figures of different scenarios
f1 = figure(1);
set(f1, 'Position', [-1715 -114 640 1060])
tiledlayout(length(G), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);
for i = 1:length(G)
    ax(i) = nexttile;
    G(i).rmt.plot_circle(ax(i));
    hold on
    G(i).rmt.plot_eigenvalue_distribution(ax(i));
    hold off
end

linkaxes(ax,'y')
