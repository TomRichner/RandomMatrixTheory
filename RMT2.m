classdef RMT2 < handle
    properties
        n
        mu % Base mean scalar
        A % Base random matrix (Gaussian, mean 0, var 1)
        f
        E
        I
        b_E
        b_I
        mu_E
        mu_I
        density
        S
        shift
        D
        M
        eigenvalues
        r
        x
        description
        zero_pop_mean
        row_sum_zero_enabled
        post_sparsification_zrs_enabled
        color_outliers
        outlier_factor
    end

    methods
        function obj = RMT2(n)
            % constructor
            obj.n = n;
            obj.A = randn(n,n);

            % defaults
            obj.mu = 0;
            obj.f = 0.5;
            obj.density = 1;
            obj.shift = 0;
            obj.b_E = 1;
            obj.b_I = 1;
            obj.mu_E = 0;
            obj.color_outliers = false;
            obj.outlier_factor = 1.02;

            % Initialize f, E, I and dependent matrices (mu_I, M, D)
            obj.set_f(0.5);

            obj.S = true(n,n);
            obj.D = eye(n);
            obj.M = zeros(n,n);
            obj.row_sum_zero_enabled = false;
            obj.post_sparsification_zrs_enabled = false;
            obj.zero_pop_mean = false;
        end

        function set_mu(obj, mu)
            obj.mu = mu;
        end

        function set_f(obj, f)
            obj.f = f;
            obj.E = false(obj.n, 1);
            obj.E(1:round(f*obj.n)) = true;
            obj.I = ~obj.E;

            % Update dependent matrices if parameters are already set
            obj.update_D();
            obj.update_M();
        end

        function set_mean_parameters(obj, mu_E)
            obj.mu_E = mu_E;
            if isempty(obj.f)
                % If f is not set, we can't calculate mu_I yet.
                % But constructor sets default f, so this branch might be moot
                % unless f was cleared. Assuming f is always valid.
                warning('RMT2:fNotSet', 'f is not set properly, using default or previous value.');
            end

            obj.update_M();
        end

        function set_stdev_parameters(obj, b_E, b_I)
            obj.b_E = b_E;
            obj.b_I = b_I;
            obj.update_D();
        end

        function update_D(obj)
            if isempty(obj.E) || isempty(obj.I)
                % Can happen if set_f hasn't been called, but constructor handles defaults
                % If needed, re-call set_f logic or just wait
                return;
            end

            diag_vals = zeros(obj.n, 1);
            diag_vals(obj.E) = obj.b_E;
            diag_vals(obj.I) = obj.b_I;
            obj.D = diag(diag_vals);
        end

        function update_M(obj)
            if isempty(obj.f) || isempty(obj.E)
                return;
            end

            obj.mu_I = -obj.f * obj.mu_E / (1-obj.f);

            obj.M = zeros(obj.n, obj.n);
            obj.M(:, obj.E) = obj.mu_E;
            obj.M(:, obj.I) = obj.mu_I;
        end

        function set_density(obj, density)
            obj.density = density;
            obj.S = rand(obj.n, obj.n) < obj.density;
        end

        function set_zero_pop_mean(obj, enable)
            obj.zero_pop_mean = enable;
        end

        function set_row_sum_zero(obj, enable)
            obj.row_sum_zero_enabled = enable;
        end

        function set_post_sparsification_zrs(obj, enable)
            obj.post_sparsification_zrs_enabled = enable;
        end

        function set_color_outliers(obj, enable)
            obj.color_outliers = enable;
        end

        function set_outlier_factor(obj, factor)
            obj.outlier_factor = factor;
        end

        function is_dense = check_S_is_dense(obj)
            % Check if S is fully dense (all ones)
            is_dense = all(obj.S(:));
        end

        function is_valid = check_S_is_binary(obj)
            % Check if S contains only 0s and 1s
            is_valid = all(obj.S(:) == 0 | obj.S(:) == 1);
        end

        function W = get_W(obj)
            % Validate mutual exclusion
            if obj.row_sum_zero_enabled && obj.post_sparsification_zrs_enabled
                error('RMT2:InvalidZRS', 'Cannot enable both row_sum_zero and post_sparsification_zrs simultaneously.');
            end

            if obj.zero_pop_mean && obj.post_sparsification_zrs_enabled
                error('RMT2:InvalidZRS', 'Cannot enable zero_pop_mean when post_sparsification_zrs is enabled.');
            end

            if obj.row_sum_zero_enabled
                % Dense ZRS: zeros A*D+M together, then adds mu
                if ~obj.check_S_is_dense()
                    error('RMT2:InvalidS', 'row_sum_zero_enabled requires S to be fully dense (all ones).');
                end

                A_D_M = obj.A * obj.D + obj.M;
                row_means = mean(A_D_M, 2);         % n×1 vector
                if obj.zero_pop_mean
                    % If enabled, mu is effectively included in the ZRS (so it is removed)
                    % Matrix has row sums exactly 0
                    W = A_D_M - row_means;
                else
                    % Standard ZRS: Center A*D+M, then add mu
                    % Matrix has row sums equal to mu
                    W = A_D_M + obj.mu - row_means;     % add mu, subtract row mean
                end

            elseif obj.post_sparsification_zrs_enabled
                % Sparse ZRS: zeros A*D at sparse locations, then adds mu+M
                if ~obj.check_S_is_binary()
                    error('RMT2:InvalidS', 'post_sparsification_zrs requires S to contain only 0s and 1s.');
                end

                % Step 1: Compute sparsified A*D
                A_D = obj.A * obj.D;
                A_D_sparse = obj.S .* A_D;

                % Step 2: Zero each row of the sparse A*D
                row_sums = sum(A_D_sparse, 2);
                row_counts = sum(obj.S, 2);
                correction = row_sums ./ max(row_counts, 1);  % avoid div by 0
                A_D_zrs = A_D_sparse - obj.S .* correction;   % only subtract at S locations

                % Step 3: Add mu and M AFTER zeroing (preserves low-rank)
                W = A_D_zrs + obj.S .* (obj.mu + obj.M);

            else
                % No ZRS: standard computation
                term = obj.A * obj.D + obj.mu + obj.M;
                W = obj.S .* term;
            end

            W = W + obj.shift * eye(obj.n);
        end


        function compute_eigenvalues(obj)
            W = obj.get_W();
            obj.eigenvalues = eig(W);
        end

        function plot_eigenvalue_distribution(obj, target_ax)
            % Plot style flag: 1 = semitransparent grey filled, 2 = black empty circles
            plot_style_flag = 2;

            if obj.color_outliers
                R = obj.compute_expected_radius();
                outlier_threshold = obj.outlier_factor * R;
                is_outlier = abs(obj.eigenvalues - obj.shift) > outlier_threshold;

                if plot_style_flag == 1
                    % Option 1: Semitransparent grey circles, filled, no edge
                    scatter(target_ax, real(obj.eigenvalues(~is_outlier)), imag(obj.eigenvalues(~is_outlier)), 30, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
                else
                    % Option 2: Black circles with no fill
                    scatter(target_ax, real(obj.eigenvalues(~is_outlier)), imag(obj.eigenvalues(~is_outlier)), 25, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none');
                end

                wasHeld = ishold(target_ax);
                hold(target_ax, 'on');

                % Plot outliers (green filled)
                scatter(target_ax, real(obj.eigenvalues(is_outlier)), imag(obj.eigenvalues(is_outlier)), 35, 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', 'none');

                if ~wasHeld
                    hold(target_ax, 'off');
                end
            else
                if plot_style_flag == 1
                    % Option 1: Semitransparent grey circles, filled, no edge
                    scatter(target_ax, real(obj.eigenvalues), imag(obj.eigenvalues), 30, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
                else
                    % Option 2: Black circles with no fill
                    scatter(target_ax, real(obj.eigenvalues), imag(obj.eigenvalues), 25, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none');
                end
            end

            axis(target_ax, 'equal');
            xlabel(target_ax, 'Re($\lambda$)', 'Interpreter', 'latex');
            ylabel(target_ax, 'Im($\lambda$)', 'Interpreter', 'latex');
        end

        function plot_discrete_eigenvalue_distribution(obj, target_ax, dt)
            % Plot eigenvalues of the equivalent discrete-time system
            % Uses matrix exponential: A_discrete = expm(W * dt)
            % For discrete systems, stability requires |lambda| < 1
            %
            % Arguments:
            %   target_ax - axes handle to plot on
            %   dt - time step (default: 0.01)

            if nargin < 3
                dt = 0.01;
            end

            % Compute discrete system matrix via matrix exponential
            W = obj.get_W();
            A_discrete = expm(W * dt);

            % Compute eigenvalues of discrete system
            discrete_eigenvalues = eig(A_discrete);

            % Plot style flag: 1 = semitransparent grey filled, 2 = black empty circles
            plot_style_flag = 2;

            if plot_style_flag == 1
                % Option 1: Semitransparent grey circles, filled, no edge
                scatter(target_ax, real(discrete_eigenvalues), imag(discrete_eigenvalues), 30, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            else
                % Option 2: Black circles with no fill
                scatter(target_ax, real(discrete_eigenvalues), imag(discrete_eigenvalues), 25, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none');
            end

            % Draw unit circle (stability boundary for discrete systems)
            wasHeld = ishold(target_ax);
            hold(target_ax, 'on');
            theta = linspace(0, 2*pi, 100);
            plot(target_ax, cos(theta), sin(theta), 'k--', 'LineWidth', 1);
            if ~wasHeld
                hold(target_ax, 'off');
            end

            axis(target_ax, 'equal');
            xlabel(target_ax, 'Re($\lambda$)', 'Interpreter', 'latex');
            ylabel(target_ax, 'Im($\lambda$)', 'Interpreter', 'latex');
            title(target_ax, sprintf('Discrete System Eigenvalues (dt = %.4f)', dt));
        end

        function max_real = get_max_real_eig(obj)
            if isempty(obj.eigenvalues)
                obj.compute_eigenvalues();
            end
            max_real = max(real(obj.eigenvalues));
        end

        function [R, total_variance] = compute_expected_radius(obj)
            % Implements Equation 18 from Harris et al. (2023)
            % R = sqrt(N * [f*sigma_se^2 + (1-f)*sigma_si^2])
            % Where sigma_sk^2 (Eq 16) accounts for both variance and mean due to sparsity.

            alpha = obj.density;

            % The mean of the entries before sparsification (dense mean)
            % For Excitatory cols: mu + mu_E
            % For Inhibitory cols: mu + mu_I
            mu_dense_E = obj.mu + obj.mu_E;
            mu_dense_I = obj.mu + obj.mu_I;

            % The standard deviation of the entries before sparsification (dense sigma)
            % For Excitatory cols: b_E
            % For Inhibitory cols: b_I
            sigma_dense_E = obj.b_E;
            sigma_dense_I = obj.b_I;

            % Variance of the sparse entries (Eq. 16/B4/B5 in paper)
            % Var = alpha * sigma_dense^2 + alpha * (1 - alpha) * mu_dense^2
            var_sparse_E = alpha * sigma_dense_E^2 + alpha * (1 - alpha) * mu_dense_E^2;
            var_sparse_I = alpha * sigma_dense_I^2 + alpha * (1 - alpha) * mu_dense_I^2;

            % Total Variance weighted by population fraction f
            total_variance = obj.f * var_sparse_E + (1 - obj.f) * var_sparse_I;

            R = sqrt(obj.n * total_variance);
        end

        function sigma_I = compute_sigma_I_from_eff(obj, sigma_eff, sigma_E)
            term = (sigma_eff^2 - obj.f * sigma_E^2) / (1 - obj.f);
            if term < 0
                error('RMT2:InvalidSigma', 'Invalid parameters: Effective sigma cannot be achieved with given sigma_E and f.');
            end
            sigma_I = sqrt(term);
        end

        function set_plot_circle(obj, r, x)
            obj.r = r;
            obj.x = x;
        end

        function plot_circle(obj, target_ax)
            pos = [obj.x - obj.r, -obj.r, 2*obj.r, 2*obj.r];
            rectangle(target_ax, 'Position', pos, 'Curvature', [1,1], 'EdgeColor', 'k', 'LineWidth', 1);
        end

        function plot_weight_histogram(obj, target_ax, num_bins)
            % Plots overlapping histograms of E and I weights
            % E weights: 50% opacity red
            % I weights: 50% opacity blue
            % No box, no bar outlines

            if nargin < 3
                num_bins = 50;
            end

            W = obj.get_W();

            % Remove diagonal (set to 0, which are removed below)
            W(1:obj.n+1:end) = 0;

            % Extract weights from E and I columns
            W_E = W(:, obj.E);
            W_I = W(:, obj.I);

            % Flatten to vectors
            weights_E = W_E(:);
            weights_I = W_I(:);

            % Remove zero weights
            weights_E = weights_E(weights_E ~= 0);
            weights_I = weights_I(weights_I ~= 0);

            % Determine common bin edges
            all_weights = [weights_E; weights_I];
            bin_edges = linspace(min(all_weights), max(all_weights), num_bins + 1);

            % Plot histograms
            hold(target_ax, 'on');

            histogram(target_ax, weights_E, bin_edges, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            histogram(target_ax, weights_I, bin_edges, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

            % Remove box
            box(target_ax, 'off');

            % Labels
            xlabel(target_ax, 'Weight');
            ylabel(target_ax, 'Count');
            legend(target_ax, {'E', 'I'}, 'Box', 'off');

            hold(target_ax, 'off');
        end

        function display_parameters(obj)
            % Display measured statistics from W and compare to set property values
            W = obj.get_W();

            % Remove diagonal for statistics
            W_no_diag = W;
            W_no_diag(1:obj.n+1:end) = NaN;

            % Extract E and I columns (excluding diagonal)
            W_E = W_no_diag(:, obj.E);
            W_I = W_no_diag(:, obj.I);

            % For sparse matrices, only consider non-zero (connected) entries
            if obj.density < 1
                mask_E = obj.S(:, obj.E);
                mask_I = obj.S(:, obj.I);
                weights_E = W_E(mask_E & ~isnan(W_E));
                weights_I = W_I(mask_I & ~isnan(W_I));
                all_weights = [weights_E; weights_I];
            else
                weights_E = W_E(~isnan(W_E));
                weights_I = W_I(~isnan(W_I));
                all_weights = W_no_diag(~isnan(W_no_diag));
            end

            % Compute measured statistics
            measured_mu = mean(all_weights);
            measured_mu_E = mean(weights_E);
            measured_mu_I = mean(weights_I);
            measured_sigma_E = std(weights_E);
            measured_sigma_I = std(weights_I);
            % Compute measured_sigma_pop as weighted average of variances (matching expected formula)
            measured_var_pop = obj.f * measured_sigma_E^2 + (1 - obj.f) * measured_sigma_I^2;
            measured_sigma_pop = sqrt(measured_var_pop);

            % Compute expected population sigma from total_variance
            [~, total_variance] = obj.compute_expected_radius();
            expected_sigma_pop = sqrt(total_variance);

            % Display results
            fprintf('\n========== RMT2 Parameter Summary ==========\n');
            if ~isempty(obj.description)
                fprintf('Description:        %s\n', obj.description);
            end
            fprintf('                    Set Value    Measured Value\n');
            fprintf('----------------------------------------------\n');
            fprintf('mu (overall mean)   %10.4f    %10.4f\n', obj.mu, measured_mu);
            fprintf('mu_E                %10.4f    %10.4f\n', obj.mu_E, measured_mu_E);
            fprintf('mu_I                %10.4f    %10.4f\n', obj.mu_I, measured_mu_I);
            fprintf('sigma_E (b_E)       %10.4f    %10.4f\n', obj.b_E, measured_sigma_E);
            fprintf('sigma_I (b_I)       %10.4f    %10.4f\n', obj.b_I, measured_sigma_I);
            fprintf('sigma_pop           %10.4f    %10.4f\n', expected_sigma_pop, measured_sigma_pop);
            fprintf('==============================================\n\n');
        end

        function new_obj = copy(obj)
            new_obj = RMT2(obj.n);
            props = properties(obj);
            for i = 1:length(props)
                if ~strcmp(props{i}, 'A') % A is already initialized in constructor with randn
                    new_obj.(props{i}) = obj.(props{i});
                end
            end
            new_obj.A = obj.A; % now copy A to ensure identical noise
        end
    end
end