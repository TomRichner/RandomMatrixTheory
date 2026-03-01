# Comparison: RMT.m (ConnectivityAdaptation) vs RMT4.m (RandomMatrixTheory)

> **Assumption**: `RMT.m` was branched from `RMT4.m` and then further developed. This comparison lists all differences. Changes labeled **RMT only** represent development done after branching.

## Summary

`RMT.m` is a **cleaned-up version** of `RMT4.m`. The substantive changes are:
- Removal of all legacy / backward-compatibility methods (`get_sparse_means`, `get_sparse_variances`, `get_theoretical_stats`, `get_Jacobian`)
- Removal of transitional comments (e.g., "This replaces the old get_Jacobian() method", "for backward compatibility")
- One minor comment improvement (clarifying *why* the diagonal is removed in `display_parameters`)

All properties, dependent property getters, setters, parameter logic, eigenvalue caching, plotting, and the `copy` method are **functionally identical** (modulo the class name rename).

## Detailed Differences

| #   | Location                                                | RMT4.m (original)                                                                                                   | RMT.m (branched, further developed)                             | Type                                             |
| --- | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- | ------------------------------------------------ |
| 1   | Class name                                              | `classdef RMT4 < handle`                                                                                            | `classdef RMT < handle`                                         | Rename                                           |
| 2   | Constructor                                             | `function obj = RMT4(N)` / `% RMT4 Constructor`                                                                     | `function obj = RMT(N)` / `% RMT Constructor`                   | Rename                                           |
| 3   | `get.W` comment (line ~120)                             | `% This replaces the old get_Jacobian() method`                                                                     | *(removed)*                                                     | Cleanup — legacy reference removed               |
| 4   | Warning ID in `get.W` ZRS case                          | `'RMT4:SparsityWarning'`                                                                                            | `'RMT:SparsityWarning'`                                         | Rename                                           |
| 5   | Legacy methods (lines 329–352 in RMT4)                  | `get_sparse_means()`, `get_sparse_variances()`, `get_theoretical_stats()`, `get_Jacobian()` — all marked DEPRECATED | *(all removed)*                                                 | **Cleanup — removed 4 legacy methods**           |
| 6   | `display_parameters` comment                            | `% Remove diagonal for statistics`                                                                                  | `% Remove diagonal for statistics since they may contain shift` | Improvement — added clarifying reason            |
| 7   | `display_parameters` header                             | `'RMT4 Parameter Summary'`                                                                                          | `'RMT Parameter Summary'`                                       | Rename                                           |
| 8   | Error ID in `compute_sigma_tilde_i_for_target_variance` | `'RMT4:SparseCaseNotSupported'`                                                                                     | `'RMT:SparseCaseNotSupported'`                                  | Rename                                           |
| 9   | Error ID in `compute_sigma_tilde_i_for_target_variance` | `'RMT4:InvalidTargetVariance'`                                                                                      | `'RMT:InvalidTargetVariance'`                                   | Rename                                           |
| 10  | `compute_eigenvalues` comment                           | `% Force recomputation and cache update (for backward compatibility)`                                               | `% Force recomputation and cache update`                        | Cleanup — removed "(for backward compatibility)" |
| 11  | `copy()` method                                         | `new_obj = RMT4(obj.N)`                                                                                             | `new_obj = RMT(obj.N)`                                          | Rename                                           |

## Categorized Changes

### Class Rename (trivial)
Rows 1, 2, 4, 7, 8, 9, 11 — all instances of `RMT4` to `RMT` in class name, constructor, error/warning IDs, display header, and `copy()`.

### Legacy Code Removed (substantive cleanup)
- **Row 5**: Four deprecated methods removed entirely:
  - `get_sparse_means()` — use `obj.mu_se`, `obj.mu_si` directly
  - `get_sparse_variances()` — use `obj.sigma_se_sq`, `obj.sigma_si_sq` directly
  - `get_theoretical_stats()` — use `obj.lambda_O`, `obj.R` directly
  - `get_Jacobian()` — use `obj.W` directly

### Comment Cleanup
- **Row 3**: Removed `"% This replaces the old get_Jacobian() method"` from `get.W`
- **Row 10**: Removed `"(for backward compatibility)"` from `compute_eigenvalues`

### Comment Improvement
- **Row 6**: Added `"since they may contain shift"` to explain *why* the diagonal is removed in `display_parameters`

## Conclusion

Your intuition is correct — **RMT.m had further development after branching from RMT4.m**. The development was focused on **cleaning up legacy cruft**: removing deprecated wrapper methods and transitional comments that referenced the old API. There are no new features or behavioral changes; the cleanup makes RMT.m a leaner, more forward-looking version of the same class.

> [!TIP]
> If `RMT4.m` no longer needs to support callers using the old API (`get_Jacobian`, `get_sparse_means`, etc.), you could apply the same cleanup to bring it in line with `RMT.m`.
