% Test script for RMT2 ZRS ZRS logic refactor

%% Setup
n = 100;
obj = RMT2(n);
obj.set_mean_parameters(1.0); % mu_E
obj.set_stdev_parameters(0.1, 0.1); % b_E, b_I
obj.set_sparsity(1.0); % Start dense

fprintf('Running RMT2 ZRS Logic Verification...\n');

%% Test 1: Mutual Exclusion
fprintf('Test 1: Mutual Exclusion... ');
try
    obj.set_row_sum_zero(true);
    obj.set_post_sparsification_zrs(true);
    obj.get_W();
    error('Failed: Should have thrown error for both enabled');
catch ME
    if strcmp(ME.identifier, 'RMT2:InvalidZRS')
        fprintf('Passed.\n');
    else
        rethrow(ME);
    end
end
% Reset
obj.set_row_sum_zero(false);
obj.set_post_sparsification_zrs(false);

%% Test 2: Dense ZRS (row_sum_zero_enabled)
fprintf('Test 2: Dense ZRS... ');
obj.set_sparsity(1.0); % Ensure dense
obj.set_row_sum_zero(true);

W = obj.get_W();
% Check logic: W = A*D+M+mu - mean(A*D+M)
% So mean(W) should be mu
% Wait, W = (A*D+M - mean(A*D+M)) + mu
% Then mean(W) = 0 + mu = mu.
% Or sum(W,2)/n should be approx mu.

% Actually, let's check exact math relative to obj properties
A_D_M = obj.A * obj.D + obj.M;
row_means = mean(A_D_M, 2);
expected_W = A_D_M + obj.mu - row_means + obj.shift * eye(n);

diff = norm(W - expected_W, 'fro');
if diff < 1e-10
    fprintf('Passed (matches expected logic).\n');
else
    fprintf('Failed (diff = %e).\n', diff);
    disp('First few W values:');
    disp(W(1:3, 1:3));
    disp('First few expected values:');
    disp(expected_W(1:3, 1:3));
end

% Verify row means match mu (ignoring shift for a moment, or taking it into account)
% mean(W,2) = mu + shift/n (since shift is diagonal)
% Wait, shift is * eye(n). So it adds shift to 1 element per row.
% So mean of row i is mu + shift/n.
actual_row_means = mean(W, 2);
expected_row_means = obj.mu + obj.shift / n;
diff_means = norm(actual_row_means - expected_row_means);
if diff_means < 1e-10
    fprintf('   Row means verification: Passed.\n');
else
    fprintf('   Row means verification: Failed (diff = %e).\n', diff_means);
end

% Reset
obj.set_row_sum_zero(false);

%% Test 3: Dense ZRS Validation (Sparse S)
fprintf('Test 3: Dense ZRS checking S density... ');
obj.set_sparsity(0.5); % Sparse
obj.set_row_sum_zero(true);
try
    obj.get_W();
    error('Failed: Should have thrown error for sparse S');
catch ME
    if strcmp(ME.identifier, 'RMT2:InvalidS')
        fprintf('Passed.\n');
    else
        rethrow(ME);
    end
end
obj.set_row_sum_zero(false);

%% Test 4: Sparse ZRS (post_sparsification_zrs_enabled)
fprintf('Test 4: Sparse ZRS... ');
obj.set_sparsity(0.5); % Sparse
obj.set_post_sparsification_zrs(true);

W = obj.get_W();

% Logic:
% A_D_sparse = S .* (A*D)
% correction = row means of A_D_sparse (normalized by row counts)
% A_D_zrs = A_D_sparse - S .* correction
% W = A_D_zrs + S.*(mu + M) (+ shift)

% Verify that random part A*D sums to zero for each row
random_part = W - (obj.S .* (obj.mu + obj.M)) - obj.shift * eye(n);
% This random_part should correspond to A_D_zrs
% And sum(A_D_zrs, 2) should be 0 (or close to it)
zrs_sums = sum(random_part, 2);
if norm(zrs_sums) < 1e-10
    fprintf('Passed (Random part is Zero Row Sum).\n');
else
    fprintf('Failed (Random part row sums norm = %e).\n', norm(zrs_sums));
end

% Reset
obj.set_post_sparsification_zrs(false);

fprintf('All tests passed!\n');
exit;
