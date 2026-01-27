classdef RMT3 < handle
    properties
        N           % System size
        alpha       % Sparsity parameter (0 < alpha <= 1)
        f           % Fraction of excitatory neurons

        % Population statistics (can be scaled or unscaled, user responsibility)
        mu_e        % Mean of excitatory population
        mu_i        % Mean of inhibitory population
        sigma_e     % Std dev of excitatory population
        sigma_i     % Std dev of inhibitory population

        % Internal matrices
        A           % Base random matrix (Gaussian, mean 0, var 1)
        S           % Sparsity mask (logical)

        % Low-rank structure M = u * v'
        u           % Left vector(s) (N x k)
        v           % Right vector(s) (N x k)

        % Control flags
        zrs_mode    % 'none', 'ZRS', 'SZRS', 'Partial_SZRS'
        shift       % Scalar shift for eigenvalues

        % Dependent/Computed properties (for caching/display)
        E           % Logical index for Excitatory
        I           % Logical index for Inhibitory

        % Visualization
        eigenvalues
        description
        outlier_threshold  % Multiplier for R to determine outlier eigenvalues (default 1.03)
    end

    methods
        function obj = RMT3(N)
            % RMT3 Constructor
            obj.N = N;

            % Defaults
            obj.alpha = 1.0;
            obj.f = 0.5;
            obj.mu_e = 0;
            obj.mu_i = 0;
            obj.sigma_e = 1; % Default unit variance
            obj.sigma_i = 1;

            obj.zrs_mode = 'none';
            obj.shift = 0;
            obj.description = '';
            obj.outlier_threshold = 1.03;

            % Initialize random matrices
            obj.A = randn(N, N); % Mean 0, Var 1
            obj.update_sparsity();

            % Initialize E/I indices
            obj.update_populations();

            % Default u and v (Dale's law structure)
            obj.set_default_uv();
        end

        %% Parameter Setters
        function set.alpha(obj, val)
            obj.alpha = val;
            if ~isempty(obj.N)
                obj.update_sparsity();
            end
        end

        function set.f(obj, val)
            obj.f = val;
            if ~isempty(obj.N)
                obj.update_populations();
            end
        end

        function set_params(obj, mu_e, mu_i, sigma_e, sigma_i, f, alpha)
            if nargin > 1, obj.mu_e = mu_e; end
            if nargin > 2, obj.mu_i = mu_i; end
            if nargin > 3, obj.sigma_e = sigma_e; end
            if nargin > 4, obj.sigma_i = sigma_i; end
            if nargin > 5, obj.f = f; end
            if nargin > 6, obj.alpha = alpha; end

            % Always update default u/v based on current params
            obj.set_default_uv();
        end

        function set_alpha(obj, alpha)
            % Set alpha (sparsity) independently without affecting other parameters
            % This will trigger update_sparsity() via the property setter
            obj.alpha = alpha;
        end

        function set_zrs_mode(obj, mode)
            valid_modes = {'none', 'ZRS', 'SZRS', 'Partial_SZRS'};
            if ~ismember(mode, valid_modes)
                error('Invalid ZRS mode. decomposition choices: %s', strjoin(valid_modes, ', '));
            end
            obj.zrs_mode = mode;
        end

        %% Internal Updates
        function update_populations(obj)
            obj.E = false(obj.N, 1);
            obj.E(1:round(obj.f * obj.N)) = true;
            obj.I = ~obj.E;
        end

        function update_sparsity(obj)
            obj.S = rand(obj.N, obj.N) < obj.alpha;
        end

        function set_default_uv(obj)
            % Default u is all ones
            obj.u = ones(obj.N, 1);

            % Default v depends on population means
            obj.v = zeros(obj.N, 1);
            if ~isempty(obj.E)
                obj.v(obj.E) = obj.mu_e;
                obj.v(obj.I) = obj.mu_i;
            end
        end

        %% Core Computation
        function W = get_Jacobian(obj)
            % Construct Diagonal Variance Matrix D (Eq 11 in Harris 2023)
            D_vec = zeros(obj.N, 1);
            D_vec(obj.E) = obj.sigma_e;
            D_vec(obj.I) = obj.sigma_i;
            D = diag(D_vec);

            % Deterministic component M (Eq 6)
            M = obj.u * obj.v';

            % Calculate W based on mode
            % Note: Harris Eq 6: W = S .* (A*D + M)

            switch obj.zrs_mode
                case 'none'
                    % Standard construction
                    W_dense = (obj.A * D) + M;
                    W = obj.S .* W_dense;

                case 'ZRS'
                    % "Dense" ZRS using Projection Operator P (Eq 24, 25)
                    % Only valid/sensible if alpha = 1 (Dense)
                    % P = I - (u*u')/N  <-- Assumes u is ones vector.
                    % Generalized P projection onto nullspace of u:
                    % P = I - u * inv(u'*u) * u'

                    if obj.alpha < 1
                        warning('Using ''ZRS'' (projection) with sparse matrix. This will destroy sparsity. Consider ''SZRS''.');
                    end

                    P = eye(obj.N) - obj.u * pinv(obj.u' * obj.u) * obj.u';

                    % W = A*D*P + M
                    % Applying P to random component A*D ensures row sums of random part are 0.
                    W = (obj.A * D * P) + M;

                    if obj.alpha < 1
                        % If user insists on sparse + ZRS, we mask it, but warning above applies.
                        W = obj.S .* W;
                    end

                case 'SZRS'
                    % Sparse Zero Row Sum (Eq 30)
                    % W = S .* (A*D + M) - B

                    % Base sparse matrix
                    W_base = obj.S .* ((obj.A * D) + M);

                    % Correction matrix B
                    % Calculate row averages of non-zero elements

                    % Sum of weights per row
                    row_sums = sum(W_base, 2);
                    % Count of non-zero elements per row (sparsity pattern)
                    row_counts = sum(obj.S, 2);

                    % Avoid division by zero
                    row_counts(row_counts == 0) = 1;

                    Avg_Wi = row_sums ./ row_counts; % Vector of row averages

                    % Subtract average from existing connections only
                    % B_ij = S_ij * Avg_Wi
                    B = obj.S .* Avg_Wi;

                    W = W_base - B;

                case 'Partial_SZRS'
                    % Partial SZRS (Eq 32)
                    % Apply correction ONLY to random component J = S .* (A*D)
                    % Keep M component (S .* M) intact to preserve imbalance.

                    J_base = obj.S .* (obj.A * D);
                    M_base = obj.S .* M;

                    % Row averages of random component J
                    J_row_sums = sum(J_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;

                    Avg_Ji = J_row_sums ./ row_counts;

                    % Correction B_partial
                    B_partial = obj.S .* Avg_Ji;

                    % W = (J_base - B_partial) + M_base
                    W = (J_base - B_partial) + M_base;
            end

            % Apply shift
            W = W + obj.shift * eye(obj.N);
        end

        %% Theoretical Predictions
        function [lambda_O, R] = get_theoretical_stats(obj)
            % Equations 15, 17 and 18 from Harris 2023
            % Requires f, mu_e/i, sigma_e/i, alpha to be set.

            % Eq 15: mu_sk = alpha * mu_k (sparse effective mean)
            mu_se = obj.alpha * obj.mu_e;
            mu_si = obj.alpha * obj.mu_i;

            % Eq 17: lambda_O = N * [f * mu_se + (1-f) * mu_si]
            lambda_O = obj.N * (obj.f * mu_se + (1-obj.f) * mu_si);

            var_se = obj.alpha * (1 - obj.alpha) * obj.mu_e^2 + obj.alpha * obj.sigma_e^2;
            var_si = obj.alpha * (1 - obj.alpha) * obj.mu_i^2 + obj.alpha * obj.sigma_i^2;

            R = sqrt(obj.N * (obj.f * var_se + (1-obj.f) * var_si));
        end

        %% Helpers
        function display_parameters(obj)
            % Display measured statistics from W and compare to set property values
            W = obj.get_Jacobian();

            % Remove diagonal for statistics (often diagonal is shifted or specific)
            W_no_diag = W;
            W_no_diag(1:obj.N+1:end) = NaN;

            % Extract E and I columns (excluding diagonal)
            W_E = W_no_diag(:, obj.E);
            W_I = W_no_diag(:, obj.I);

            % For sparse matrices, only consider non-zero (connected) entries
            W_E_vals = W_E(~isnan(W_E) & (W_E ~= 0));
            W_I_vals = W_I(~isnan(W_I) & (W_I ~= 0));
            all_vals = [W_E_vals; W_I_vals];

            % Compute measured statistics
            measured_mu = mean(all_vals);
            measured_mu_E = mean(W_E_vals);
            measured_mu_I = mean(W_I_vals);
            measured_sigma_E = std(W_E_vals,1);
            measured_sigma_I = std(W_I_vals,1);

            [lambda_O, R] = obj.get_theoretical_stats();

            fprintf('\n========== RMT3 Parameter Summary ==========\n');
            if ~isempty(obj.description)
                fprintf('Description:        %s\n', obj.description);
            end
            fprintf('Mode:               %s\n', obj.zrs_mode);
            fprintf('Alpha:              %0.2f\n', obj.alpha);
            fprintf('                    Set Value    Measured (NZ)\n');
            fprintf('----------------------------------------------\n');
            fprintf('mu_E                %10.4f    %10.4f\n', obj.mu_e, measured_mu_E);
            fprintf('mu_I                %10.4f    %10.4f\n', obj.mu_i, measured_mu_I);
            fprintf('sigma_E             %10.4f    %10.4f\n', obj.sigma_e, measured_sigma_E);
            fprintf('sigma_I             %10.4f    %10.4f\n', obj.sigma_i, measured_sigma_I);
            fprintf('Theoretical Outlier: %10.4f\n', lambda_O);
            fprintf('Theoretical Radius:  %10.4f\n', R);
            fprintf('==============================================\n\n');
        end


        function compute_eigenvalues(obj)
            W = obj.get_Jacobian();
            obj.eigenvalues = eig(W);
        end

        function plot_spectrum(obj, ax)
            if nargin < 2
                figure; ax = gca;
            end

            if isempty(obj.eigenvalues)
                obj.compute_eigenvalues();
            end

            % Plot eigenvalues
            plot(ax, real(obj.eigenvalues), imag(obj.eigenvalues), 'ko', 'MarkerFaceColor', 'none');
            hold(ax, 'on');

            % Plot theoretical radius
            [lambda_O, R] = obj.get_theoretical_stats();
            theta = linspace(0, 2*pi, 100);

            xc = obj.shift;
            yc = 0;
            plot(ax, xc + R*cos(theta), yc + R*sin(theta), 'k-', 'LineWidth', 2);

            % Plot outlier eigenvalues (beyond outlier_threshold*R) as green filled circles
            threshold = obj.outlier_threshold * R;
            distances = abs(obj.eigenvalues - xc - 1i*yc);
            outlier_mask = distances > threshold;
            outlier_eigs = obj.eigenvalues(outlier_mask);
            if ~isempty(outlier_eigs)
                plot(ax, real(outlier_eigs), imag(outlier_eigs), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 .7 0], 'MarkerEdgeColor', [0 .7 0]);
            end

            xlabel(ax, 'Re(\lambda)');
            ylabel(ax, 'Im(\lambda)');
            grid(ax, 'on');
            axis(ax, 'equal');
            hold(ax, 'off');
        end
        function new_obj = copy(obj)
            % Create new object
            new_obj = RMT3(obj.N);

            % Copy all properties
            meta = metaclass(obj);
            props = meta.PropertyList;
            for i = 1:numel(props)
                prop = props(i);
                if ~prop.Dependent && ~prop.Constant && ~prop.Abstract
                    new_obj.(prop.Name) = obj.(prop.Name);
                end
            end
        end
    end
end
