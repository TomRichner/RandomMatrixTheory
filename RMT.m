classdef RMT < handle
    properties
        n
        b
        mu
        A
        dense_mask
        density
        E
        I
        f
        mu_E
        mu_I
        M
        eigenvalues
        r
        x
        description
    end
    
    methods
        function obj = RMT(n, b, mu)
            % constructor
            obj.n = n;
            obj.b = b;
            obj.mu = mu;
            obj.A = b*randn(n,n)+mu;
            obj.density = 1; % default to 1
            obj.dense_mask = true(n,n);
        end
        
        function apply_sparsity(obj, mean_indegree)
            obj.density = mean_indegree / obj.n;
            obj.dense_mask = rand(obj.n, obj.n) < obj.density;
            obj.A(~obj.dense_mask) = 0;
        end
        
        function set_rajan_means(obj, f, mu_E)
            obj.f = f;
            obj.mu_E = mu_E;
            
            obj.E = false(obj.n, 1);
            obj.E(1:round(f*obj.n)) = true;
            obj.I = ~obj.E;
            
            obj.mu_I = -f * mu_E / (1-f);
            
            obj.M = zeros(obj.n, obj.n);
            obj.M(:, obj.E) = obj.mu_E;
            obj.M(:, obj.I) = obj.mu_I;
            
            obj.A = obj.A + obj.M;
        end
        
        function row_sum_to_zero(obj)
            Z_in = obj.A;
            Z_out = Z_in;
            for i = 1:obj.n
                row = Z_out(i, :);
                current_sum = sum(row);
                
                adjustment_mask = obj.dense_mask(i, :);
                num_to_adjust = sum(adjustment_mask);
                
                if num_to_adjust > 0
                    adjustment = -current_sum / num_to_adjust;
                    row(adjustment_mask) = row(adjustment_mask) + adjustment;
                    Z_out(i, :) = row;
                end
            end
            obj.A = Z_out;
        end
        
        function shift_diagonal(obj, shift)
            W = shift .* eye(obj.n, obj.n);
            obj.A = obj.A + W;
        end
        
        function compute_eigenvalues(obj)
            obj.eigenvalues = eig(obj.A);
        end
        
        function plot_eigenvalue_distribution(obj, target_ax)
            scatter(target_ax, real(obj.eigenvalues), imag(obj.eigenvalues), 18, 'MarkerEdgeColor',[0 0 0]);
            axis(target_ax, 'equal');
            xlabel(target_ax, 'Re($\lambda$)', 'Interpreter', 'latex');
            ylabel(target_ax, 'Im($\lambda$)', 'Interpreter', 'latex');
            title(target_ax, obj.description, 'FontWeight', 'normal');
        end
        
        function max_real = get_max_real_eig(obj)
            if isempty(obj.eigenvalues)
                obj.compute_eigenvalues();
            end
            max_real = max(real(obj.eigenvalues));
        end
        
        function set_plot_circle(obj, r, x)
            obj.r = r;
            obj.x = x;
        end
        
        function plot_circle(obj, target_ax)
             pos = [obj.x - obj.r, -obj.r, 2*obj.r, 2*obj.r];
             rectangle(target_ax, 'Position', pos, 'Curvature', [1,1], 'EdgeColor', 'k', 'LineWidth', 1);
        end

        function add_constant(obj, const)
            obj.A = obj.A + const;
        end

        function new_obj = copy(obj)
            new_obj = RMT(obj.n, obj.b, obj.mu);
            props = properties(obj);
            for i = 1:length(props)
                if ~strcmp(props{i}, 'A') % A is already initialized in constructor with randn
                    new_obj.(props{i}) = obj.(props{i});
                end
            end
            new_obj.A = obj.A; % now copy A
        end
    end
    
end
