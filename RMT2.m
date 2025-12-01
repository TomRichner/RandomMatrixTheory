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
        AD_ZRS
        eigenvalues
        r
        x
        description
        row_sum_zero_enabled
        post_sparsification_zrs_enabled
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
            
            obj.S = true(n,n);
            obj.D = eye(n);
            obj.M = zeros(n,n);
            obj.AD_ZRS = zeros(n,n);
            obj.row_sum_zero_enabled = false;
            obj.post_sparsification_zrs_enabled = false;
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

        function set_sparsity(obj, density)
            obj.density = density;
            obj.S = rand(obj.n, obj.n) < obj.density;
        end
        
        function set_row_sum_zero(obj, enable)
            obj.row_sum_zero_enabled = enable;
        end

        function set_post_sparsification_zrs(obj, enable)
            obj.post_sparsification_zrs_enabled = enable;
        end

        function compute_AD_ZRS(obj)
            Z_in = obj.A * obj.D + obj.mu;
            vec = mean(Z_in, 2);
            obj.AD_ZRS = repmat(vec, 1, obj.n);
        end
        
        function W = get_W(obj)
            if obj.row_sum_zero_enabled
                obj.compute_AD_ZRS();
                term = (obj.A * obj.D + obj.mu - obj.AD_ZRS + obj.M);
            else
                term = (obj.A * obj.D + obj.mu + obj.M);
            end
            
            W_temp = obj.S .* term;
            
            if obj.post_sparsification_zrs_enabled
                 % Calculate row sums of the sparsified matrix
                 row_sums = sum(W_temp, 2);
                 
                 % Calculate number of non-zero elements per row in S
                 row_counts = sum(obj.S ~= 0, 2);
                 
                 % Compute correction per non-zero element
                 correction_vals = zeros(obj.n, 1);
                 mask_nz = row_counts > 0;
                 correction_vals(mask_nz) = row_sums(mask_nz) ./ row_counts(mask_nz);
                 
                 % Create W_ZRS matrix
                 W_ZRS = obj.S .* repmat(correction_vals, 1, obj.n);
                 
                 W_temp = W_temp - W_ZRS;
            end
            
            W = W_temp + obj.shift * eye(obj.n);
        end
        
        function compute_eigenvalues(obj)
            W = obj.get_W();
            obj.eigenvalues = eig(W);
        end
        
        function plot_eigenvalue_distribution(obj, target_ax)
            % Option 1: Semitransparent grey circles, filled, no edge
            scatter(target_ax, real(obj.eigenvalues), imag(obj.eigenvalues), 36, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
            
            % Option 2: Black circles with no fill
            % scatter(target_ax, real(obj.eigenvalues), imag(obj.eigenvalues), 36, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none');
            
            axis(target_ax, 'equal');
            xlabel(target_ax, 'Re($\lambda$)', 'Interpreter', 'latex');
            ylabel(target_ax, 'Im($\lambda$)', 'Interpreter', 'latex');
        end
        
        function max_real = get_max_real_eig(obj)
            if isempty(obj.eigenvalues)
                obj.compute_eigenvalues();
            end
            max_real = max(real(obj.eigenvalues));
        end

        function R = compute_expected_radius(obj)
            R = sqrt(obj.density * obj.n * (obj.f * obj.b_E^2 + (1 - obj.f) * obj.b_I^2));
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
