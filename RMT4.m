classdef RMT4 < handle
    % RMT4 - Random Matrix Theory class following Harris et al. (2023)
    % Variable names match the paper's notation for equations 15-18, 24-25, 30-31

    properties
        N               % System size
        alpha           % Sparsity/connection probability (0 < alpha <= 1)
        f               % Fraction of excitatory neurons

        % Normalized population statistics (tilde notation from Harris 2023)
        % These are the pre-sparsity parameters: mu_tilde = mu/sqrt(N), sigma_tilde = sigma/sqrt(N)
        mu_tilde_e      % Normalized mean of excitatory population
        mu_tilde_i      % Normalized mean of inhibitory population
        sigma_tilde_e   % Normalized std dev of excitatory population
        sigma_tilde_i   % Normalized std dev of inhibitory population

        % Internal matrices
        A               % Base random matrix (Gaussian, mean 0, var 1)
        S               % Sparsity mask (logical)

        % Low-rank structure M = u * v' (Eq 12)
        u               % Left vector: ones(N,1)
        v               % Right vector: [mu_tilde_e repeated Nf times, mu_tilde_i repeated N(1-f) times]

        % Control flags
        zrs_mode        % 'none', 'ZRS', 'SZRS', 'Partial_SZRS'
        shift           % Scalar shift for eigenvalues (diagonal shift)

        % Population indices
        E               % Logical index for Excitatory neurons
        I               % Logical index for Inhibitory neurons

        % Visualization
        eigenvalues
        description
        outlier_threshold  % Multiplier for R to determine outlier eigenvalues (default 1.03)
    end

    methods
        function obj = RMT4(N)
            % RMT4 Constructor
            obj.N = N;

            % Defaults
            obj.alpha = 1.0;
            obj.f = 0.5;
            obj.mu_tilde_e = 0;
            obj.mu_tilde_i = 0;
            obj.sigma_tilde_e = 1/sqrt(N);  % Default: unit variance when scaled by sqrt(N)
            obj.sigma_tilde_i = 1/sqrt(N);

            obj.zrs_mode = 'none';
            obj.shift = 0;
            obj.description = '';
            obj.outlier_threshold = 1.03;

            % Initialize random matrices
            obj.A = randn(N, N);  % Mean 0, Var 1
            obj.update_sparsity();

            % Initialize E/I indices
            obj.update_populations();

            % Default u and v (Dale's law structure, Eq 12)
            obj.set_default_uv();
        end

        %% Property Setters with auto-updates
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

        %% Parameter Setters
        function set_params(obj, mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha)
            % Set all parameters in Harris 2023 notation
            % Arguments: mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha
            if nargin > 1, obj.mu_tilde_e = mu_tilde_e; end
            if nargin > 2, obj.mu_tilde_i = mu_tilde_i; end
            if nargin > 3, obj.sigma_tilde_e = sigma_tilde_e; end
            if nargin > 4, obj.sigma_tilde_i = sigma_tilde_i; end
            if nargin > 5, obj.f = f; end
            if nargin > 6, obj.alpha = alpha; end

            % Update u/v based on current params
            obj.set_default_uv();
        end

        function set_alpha(obj, alpha)
            % Set alpha (sparsity) independently
            obj.alpha = alpha;
        end

        function set_zrs_mode(obj, mode)
            % Set the zero row-sum mode
            valid_modes = {'none', 'ZRS', 'SZRS', 'Partial_SZRS'};
            if ~ismember(mode, valid_modes)
                error('Invalid ZRS mode. Valid choices: %s', strjoin(valid_modes, ', '));
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
            % Default u and v vectors (Eq 12 from Harris 2023)
            % u = (1, ..., 1)^T
            % v = (mu_tilde_e repeated Nf times, mu_tilde_i repeated N(1-f) times)^T
            obj.u = ones(obj.N, 1);

            obj.v = zeros(obj.N, 1);
            if ~isempty(obj.E)
                obj.v(obj.E) = obj.mu_tilde_e;
                obj.v(obj.I) = obj.mu_tilde_i;
            end
        end

        %% Sparse Statistics (Eq 15, 16)
        function [mu_se, mu_si] = get_sparse_means(obj)
            % Eq 15: mu_sk = alpha * mu_tilde_k
            mu_se = obj.alpha * obj.mu_tilde_e;
            mu_si = obj.alpha * obj.mu_tilde_i;
        end

        function [sigma_se_sq, sigma_si_sq] = get_sparse_variances(obj)
            % Eq 16: sigma_sk^2 = alpha*(1-alpha)*mu_tilde_k^2 + alpha*sigma_tilde_k^2
            sigma_se_sq = obj.alpha * (1 - obj.alpha) * obj.mu_tilde_e^2 + obj.alpha * obj.sigma_tilde_e^2;
            sigma_si_sq = obj.alpha * (1 - obj.alpha) * obj.mu_tilde_i^2 + obj.alpha * obj.sigma_tilde_i^2;
        end

        %% Theoretical Predictions (Eq 17, 18)
        function [lambda_O, R] = get_theoretical_stats(obj)
            % Eq 17: lambda_O = N * [f * mu_se + (1-f) * mu_si]
            % Eq 18: R = sqrt(N * [f * sigma_se^2 + (1-f) * sigma_si^2])

            [mu_se, mu_si] = obj.get_sparse_means();
            [sigma_se_sq, sigma_si_sq] = obj.get_sparse_variances();

            % Eq 17: Outlier eigenvalue
            lambda_O = obj.N * (obj.f * mu_se + (1 - obj.f) * mu_si);

            % Eq 18: Spectral radius
            R = sqrt(obj.N * (obj.f * sigma_se_sq + (1 - obj.f) * sigma_si_sq));
        end

        %% Core Computation: Build Jacobian/Weight Matrix
        function W = get_Jacobian(obj)
            % Construct weight matrix W based on Harris 2023 equations

            % Construct Diagonal Variance Matrix D (Eq 11)
            D_vec = zeros(obj.N, 1);
            D_vec(obj.E) = obj.sigma_tilde_e;
            D_vec(obj.I) = obj.sigma_tilde_i;
            D = diag(D_vec);

            % Deterministic component M = u * v' (low-rank structure from Eq 12)
            M = obj.u * obj.v';

            switch obj.zrs_mode
                case 'none'
                    % Standard construction: W = S .* (A*D + M)  (Eq 6)
                    W_dense = (obj.A * D) + M;
                    W = obj.S .* W_dense;

                case 'ZRS'
                    % Dense ZRS using Projection Operator P (Eq 24, 25)
                    % Eq 24: P = I_N - (u*u')/N
                    % Eq 25: W = A*D*P + u*v'

                    if obj.alpha < 1
                        warning('RMT4:SparsityWarning', 'Using ''ZRS'' (projection) with sparse matrix. This will destroy sparsity. Consider ''SZRS''.');
                    end

                    % Eq 24: Projection operator
                    P = eye(obj.N) - (obj.u * obj.u') / obj.N;

                    % Eq 25: W = A*D*P + M
                    W = (obj.A * D * P) + M;

                    if obj.alpha < 1
                        W = obj.S .* W;
                    end

                case 'SZRS'
                    % Sparse Zero Row Sum (Eq 30, 31)
                    % Eq 30: W = S .* (A*D + u*v') - B
                    % Eq 31: W_bar_i = sum_j W_ij / sum_j S_ij

                    % Base sparse matrix
                    W_base = obj.S .* ((obj.A * D) + M);

                    % Eq 31: Row averages of non-zero elements
                    row_sums = sum(W_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;  % Avoid division by zero
                    W_bar_i = row_sums ./ row_counts;

                    % Correction matrix B: B_ij = S_ij * W_bar_i
                    B = obj.S .* W_bar_i;

                    % Eq 30: Final matrix
                    W = W_base - B;

                case 'Partial_SZRS'
                    % Partial SZRS (Eq 32)
                    % Apply correction ONLY to random component J = S .* (A*D)
                    % Keep M component (S .* M) intact to preserve imbalance

                    % Random component
                    J_base = obj.S .* (obj.A * D);

                    % Mean structure component
                    M_base = obj.S .* M;

                    % Eq 32: Row averages of random component J only
                    J_row_sums = sum(J_base, 2);
                    row_counts = sum(obj.S, 2);
                    row_counts(row_counts == 0) = 1;
                    J_bar_i = J_row_sums ./ row_counts;

                    % Partial correction B
                    B_partial = obj.S .* J_bar_i;

                    % W = (J_base - B_partial) + M_base
                    W = (J_base - B_partial) + M_base;
            end

            % Apply diagonal shift
            W = W + obj.shift * eye(obj.N);
        end

        %% Display and Diagnostics
        function display_parameters(obj)
            % Display measured statistics from W and compare to theoretical predictions
            W = obj.get_Jacobian();

            % Remove diagonal for statistics
            W_no_diag = W;
            W_no_diag(1:obj.N+1:end) = NaN;

            % Extract E and I columns
            W_E = W_no_diag(:, obj.E);
            W_I = W_no_diag(:, obj.I);

            % For sparse matrices, only consider non-zero entries
            W_E_vals = W_E(~isnan(W_E) & (W_E ~= 0));
            W_I_vals = W_I(~isnan(W_I) & (W_I ~= 0));

            % Measured statistics (of non-zero entries)
            measured_mu_E = mean(W_E_vals);
            measured_mu_I = mean(W_I_vals);
            measured_sigma_E = std(W_E_vals, 1);
            measured_sigma_I = std(W_I_vals, 1);

            % Theoretical predictions
            [mu_se, mu_si] = obj.get_sparse_means();
            [sigma_se_sq, sigma_si_sq] = obj.get_sparse_variances();
            [lambda_O, R] = obj.get_theoretical_stats();

            fprintf('\n========== RMT4 Parameter Summary ==========\n');
            if ~isempty(obj.description)
                fprintf('Description:          %s\n', obj.description);
            end
            fprintf('Mode:                 %s\n', obj.zrs_mode);
            fprintf('N:                    %d\n', obj.N);
            fprintf('alpha:                %.4f\n', obj.alpha);
            fprintf('f:                    %.4f\n', obj.f);
            fprintf('\n--- Harris 2023 Notation ---\n');
            fprintf('                      Set Value    Sparse Eff.  Measured(NZ)\n');
            fprintf('--------------------------------------------------------------\n');
            fprintf('mu_tilde_e            %10.4f   mu_se=%.4f   %.4f\n', obj.mu_tilde_e, mu_se, measured_mu_E);
            fprintf('mu_tilde_i            %10.4f   mu_si=%.4f   %.4f\n', obj.mu_tilde_i, mu_si, measured_mu_I);
            fprintf('sigma_tilde_e         %10.4f   sigma_se=%.4f  %.4f\n', obj.sigma_tilde_e, sqrt(sigma_se_sq), measured_sigma_E);
            fprintf('sigma_tilde_i         %10.4f   sigma_si=%.4f  %.4f\n', obj.sigma_tilde_i, sqrt(sigma_si_sq), measured_sigma_I);
            fprintf('\n--- Theoretical Predictions (Eq 17, 18) ---\n');
            fprintf('lambda_O (outlier):   %.4f\n', lambda_O);
            fprintf('R (radius):           %.4f\n', R);
            fprintf('==============================================\n\n');
        end

        %% Eigenvalue Computation
        function compute_eigenvalues(obj)
            W = obj.get_Jacobian();
            obj.eigenvalues = eig(W);
        end

        %% Plotting
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

            % Plot theoretical radius (Eq 18)
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

        %% Deep Copy
        function new_obj = copy(obj)
            % Create new object with same N
            new_obj = RMT4(obj.N);

            % Copy all non-dependent, non-constant properties
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
