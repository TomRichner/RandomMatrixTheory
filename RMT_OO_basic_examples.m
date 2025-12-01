close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
make_imagesc_nans_transparent = false; % if false, NaNs are treated as zeros
n = 500;
scale_n = 1.5; % scale factor for larger matrices
assert(mod(scale_n * n, 1) == 0, 'scale_n * n must be an integer');
b = 1;
mu = 0;
f = 0.5;
mu_E = 1; % relative to b
mean_indegree = n; % desired mean number of connections in and out if sparse

% % % Figures 1-5 examples of changing n, changing density, changing mu, changing mu_E/mu_I and applying row zero
%% Struct to hold RMT objects
G = struct();

%% 1 Base Ginibre ensemble
g = 1;
G(g).rmt = RMT(n, b, mu);
G(g).rmt.description = ['n=' num2str(n) ', dense, \mu = ' num2str(mu)];

%% 2 Bigger Ginibre ensemble
g = g+1;
n_large = scale_n * n;
G(g).rmt = RMT(n_large, b, mu);
G(g).rmt.description = ['n=' num2str(n_large) ', dense, \mu = ' num2str(mu)];

%% 3 Sparse
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.apply_sparsity(mean_indegree);
G(g).rmt.description = ['n=' num2str(n_large) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu)];

%% 4 EI imbalanced
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
mu2 = 0.04;
G(g).rmt.add_constant(mu2);
G(g).rmt.description = ['n=' num2str(n_large) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu2)];

%% 5 Rajan
g = g+1;
G(g).rmt = G(1).rmt.copy();
G(g).rmt.set_rajan_means(f, mu_E);
G(g).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1'};

%% 6 Rajan zero sum
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.row_sum_to_zero();
G(g).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed'};

%% Ginibre with shift to edge
% g = g+1;
% G(g).rmt = G(1).rmt.copy();
% shift = -G(g).rmt.b*sqrt(G(g).rmt.n)*sqrt(G(g).rmt.density);
% G(g).rmt.shift_diagonal(shift);
% G(g).rmt.description = 'Ginibre with shift';



%% Compute eigenvalues and set circles for all
for i_g=1:length(G)
    G(i_g).rmt.compute_eigenvalues();
    r = G(i_g).rmt.b*sqrt(G(i_g).rmt.n)*sqrt(G(i_g).rmt.density);
    x_center = 0;
    G(i_g).rmt.set_plot_circle(r, x_center);
end

%% Make figures of different scenarios
f1 = figure(1);
set(f1, 'Position', [100 200 710 1060], 'Color', 'white')
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax(i_g) = nexttile;
    G(i_g).rmt.plot_circle(ax(i_g));
    hold on
    G(i_g).rmt.plot_eigenvalue_distribution(ax(i_g));
    hold off
end

linkaxes(ax,'xy')

% ensure left and right xlims are the same and +/- max(abs(xlim))
xl = xlim(ax(1));
new_lim = max(abs(xl));
xlim(ax, [-new_lim, new_lim]);

% Customize axes after linking
for i_g = 1:length(G)
    axes(ax(i_g));
    x_lim = xlim;
    y_lim = ylim;
    axis off;
    
    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    text(x_lim(1), y_lim(1), G(i_g).rmt.description, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end

% % % Figures 11-15 new set of figures, similar to figures 1-5 above in content.  Except this figure will have 1 x 4 subplots.  
% subplot 1 will be H(2), subplot 2 will be H(3), subplot 3 will be H(4), and subplot 4 will be H(5).  Note that H(1) will NOT be in the figures.
% first, we need a RMT with mu_E and mu_I similar to G(6) above, call it H(1)
% In this new set of figure, we want to simulate the sparse activation of a network
% Lets assume the neurons have an internal bias vector bias_b = randn(n,1) which is due to their own intrinsic excitability and average input from other neurons
% Neurons can't have negative firing rates, so if the bias is negative, we must zero out the corresponding columns of the A matrix to form H(2)
% for H(3) we compute the actual spectral radius of H(2) and shift H(3) such that it is on the edge of stability using RMT.shift_diagonal()
% this is to simulate homeostatic mechanisms which put it on the EOC
% now we add a random, sparse stimulus u to bias_b, which simulates a change in environmental stimulus. This instantaneously changes the activation of the 
% neurons.  Copy H(1) to make H(4), then zero out its columns for negative u+bias_b values.  Then shift it the same amount as we shifted H(2) to get H(3)
% this is because it takes a while for homeostatic mechanisms to change, so the shift can't change in time.
% due to the imbalance caused by u, we should now have outlier eigenvalues in H(4)
% for H(5) copy H(4), reduce the magnitude of all non-diagonal elements enough such that the new spectral radius is on the edge of stability.
% figure 11 are the distribution of eigenvalues

%% Struct to hold RMT objects for figures 11-15
H = struct();

%% H(1) - Base Rajan with row zero sum (similar to G(6))
h = 1;
H(h).rmt = G(1).rmt.copy();
H(h).rmt.set_rajan_means(f, mu_E);
H(h).rmt.row_sum_to_zero();
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'baseline'};

%% Create bias vector
bias_b = randn(n, 1);

%% H(2) - Zero columns where bias is negative
h = h+1;
H(h).rmt = H(1).rmt.copy();
negative_bias = bias_b < 0;
H(h).rmt.zero_columns(negative_bias);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed (neg bias)'};
% Get actual spectral radius of H(2) to use for all circles
H(2).rmt.compute_eigenvalues();
r_H2 = max(real(H(2).rmt.eigenvalues));
H(h).rmt.set_plot_circle(r_H2, 0);

%% H(3) - zero rows and Shift H(2) to edge of stability
h = h+1;
H(h).rmt = H(2).rmt.copy();
H(h).rmt.row_sum_to_zero();
H(h).rmt.compute_eigenvalues();
r_H3 = max(real(H(3).rmt.eigenvalues));
H(h).rmt.shift_diagonal(-r_H3);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed, shifted to EOC'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);

%% Create sparse stimulus u
u_density = 0.2;  % 10% of neurons receive stimulus
u = zeros(n, 1);
stimulated_neurons = rand(n, 1) < u_density;
u(stimulated_neurons) = randn(sum(stimulated_neurons), 1);

%% H(4) - Copy H(1), zero columns based on u+bias_b, apply same shift
h = h+1;
H(h).rmt = H(1).rmt.copy();
negative_activation = (u + bias_b) < 0;
H(h).rmt.zero_columns(negative_activation);
H(h).rmt.shift_diagonal(-r_H3);  % Same shift as H(3)
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed (stimulus), shifted'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);

%% H(5) - Scale non-diagonal elements to edge of stability
h = h+1;
H(h).rmt = H(4).rmt.copy();
H(h).rmt.scale_nondiagonal_to_spectral_radius(r_H2);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed, scaled to EOC'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);



%% Figure 11: Eigenvalue distributions
f11 = figure(11);
set(f11, 'Position', [100 200 1420 530], 'Color', 'white')
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax11 = gobjects(4, 1);
for i_h = 2:length(H)
    H(i_h).rmt.compute_eigenvalues();
    ax11(i_h-1) = nexttile;
    H(i_h).rmt.plot_circle(ax11(i_h-1));
    hold on
    H(i_h).rmt.plot_eigenvalue_distribution(ax11(i_h-1));
    hold off
end

linkaxes(ax11,'xy')

% ensure left and right xlims are the same and +/- max(abs(xlim))
xl = xlim(ax11(1));
new_lim = max(abs(xl));
xlim(ax11, [-new_lim, new_lim]);

% Customize axes after linking
for i_h = 2:length(H)
    axes(ax11(i_h-1));
    x_lim = xlim;
    y_lim = ylim;
    axis off;
    
    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    text(x_lim(1), y_lim(1), H(i_h).rmt.description, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end



