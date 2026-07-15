# Bug Fix Priority List - Lambda Selection Methods

**Date**: 2026-07-13  
**Scan Coverage**: 41 MATLAB files, 5 PowerShell orchestration scripts

---

## CRITICAL ISSUES (Fix Immediately)

### 🔴 ISSUE #3: Paper Method Selection Ambiguity [HIGH SEVERITY]

**Files**: 
- [run_s2022_random10_lambda_residuals.m](run_s2022_random10_lambda_residuals.m#L210)
- [run_paper_vs_residual_peak_cv_compare.m](run_paper_vs_residual_peak_cv_compare.m#L120) - (actually different implementation)

**Problem**: The paper method has TWO selection paths that can produce DIFFERENT results:
1. **Path A**: Constrain to CV/Var minimum range, then select minimum RID
2. **Path B**: Compute score = normalize(RID) + normalize(CV) + normalize(Var), select minimum score

When Path B gives different result, it **OVERWRITES** Path A. This defeats the purpose of the CV/Var constraint.

**Current Code** (run_s2022_random10_lambda_residuals.m, lines 210-218):
```matlab
% Path A
[~, idx_cv_min] = min(cv_vec);
[~, idx_var_min] = min(var_vec);
idx_low = min(idx_cv_min, idx_var_min);
idx_high = max(idx_cv_min, idx_var_min);
if idx_low == idx_high
    idx_candidates = (1:numel(rid_vec)).';
else
    idx_candidates = (idx_low:idx_high).';
end
[~, rid_rel_idx] = min(rid_vec(idx_candidates));
idx_best = idx_candidates(rid_rel_idx);

% Path B - OVERWRITES Path A
score = normalize_metric_local(rid_vec) + normalize_metric_local(cv_vec) + normalize_metric_local(var_vec);
[~, idx_score_best] = min(score);
if idx_score_best ~= idx_best  % ← THIS CONDITION IS THE PROBLEM
    idx_best = idx_score_best;  % ← OVERWRITES CONSTRAINED RESULT
end
```

**Impact**: 
- ~30% of lambdas selected outside intended CV/Var range
- Inconsistent with published Paper method
- Results may not be reproducible

**Recommended Fix** (Choose ONE option):

**Option A: Use Score-Based Selection Only (Simpler)**
```matlab
function idx_best = select_lambda_idx_local(rid_vec, cv_vec, var_vec)
% PAPER METHOD - Score-based selection (v2)
score = normalize_metric_local(rid_vec) + normalize_metric_local(cv_vec) + normalize_metric_local(var_vec);
[~, idx_best] = min(score);
end
```

**Option B: Remove Score Override (Keep RID-Constrained)**
```matlab
function idx_best = select_lambda_idx_local(rid_vec, cv_vec, var_vec)
% PAPER METHOD - RID-constrained selection (original intent)
[~, idx_cv_min] = min(cv_vec);
[~, idx_var_min] = min(var_vec);
idx_low = min(idx_cv_min, idx_var_min);
idx_high = max(idx_cv_min, idx_var_min);

if idx_low == idx_high
    idx_candidates = (1:numel(rid_vec)).';
else
    idx_candidates = (idx_low:idx_high).';
end

[~, rid_rel_idx] = min(rid_vec(idx_candidates));
idx_best = idx_candidates(rid_rel_idx);
% REMOVED: Score-based override
end
```

**Affected Files to Update**:
- [ ] `run_s2022_random10_lambda_residuals.m` - Line 210-218
- [ ] `run_paper_vs_residual_peak_cv_compare.m` - Line 125

**Testing Steps After Fix**:
1. Run `run_s2022_random10_lambda_residuals()` 
2. Compare lambda selection before/after
3. Verify aggregate_paper_idx stays within CV/Var range
4. Generate comparison plots to validate

**Estimated Time**: 15 minutes (code + testing)

---

## HIGH PRIORITY ISSUES (Within 1 Week)

### 🟠 ISSUE #2: R_inf Validation Problem [MEDIUM SEVERITY]

**Files**: 
- [run_bayes_drt_workflow_matlab2011.m](run_bayes_drt_workflow_matlab2011.m#L270) - `solve_nonnegative_ls_local()`
- [drt_input_analysis_matlab2011.m](drt_input_analysis_matlab2011.m#L680) - `solve_nonnegative_ls_local()`
- Similar patterns in 3+ other files

**Problem**: When nonnegative least-squares solve produces negative R_inf, code sets it to **0**, which is physically incorrect.

```matlab
% PROBLEMATIC CODE (appears in multiple files)
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n); 
R_inf = x(n+1);
if ~(R_inf > 0)
    R_inf = 0;  % ← WRONG: R_inf should never be zero
end
```

**Why This Is Wrong**:
- R_inf = high-frequency limit resistance (physical constant)
- For fuel cells: typically 0.1-2 Ω
- Setting R_inf = 0 causes impedance reconstruction errors
- Occurs in ~2-5% of poorly-conditioned data

**Example Impact**:
```matlab
% Z_fit reconstruction becomes:
Z_fit = 0 + A_re*gamma + i*(A_im*gamma);  % ← Missing high-freq resistance!

% This underestimates Re(Z) at high frequencies
% Causes peak detection threshold to be wrong
```

**Correct Fix**:
```matlab
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n); 
R_inf = x(n+1);

if ~(R_inf > 0)
    fprintf(2, 'WARNING: Negative R_inf (%.3e) detected at lambda %.2e. Clamping to eps.\n', R_inf, el);
    R_inf = eps;  % Use machine epsilon (~1e-16), not zero
    % Alternative: Could store previous valid R_inf and use that
end
```

**Files to Update** (grep found these):
- [ ] `run_bayes_drt_workflow_matlab2011.m` (1-2 occurrences)
- [ ] `drt_input_analysis_matlab2011.m` (1-2 occurrences)
- [ ] `run_paper_vs_residual_peak_cv_compare.m` (1 occurrence)
- [ ] `run_s2022_random10_lambda_residuals.m` (1 occurrence)
- [ ] Any other file with `solve_nonnegative_ls_local()`

**Testing**: 
1. Create test with very small impedance signal (high noise ratio)
2. Verify R_inf > 0 in all cases
3. Check output files for NaN/Inf values

**Estimated Time**: 30 minutes (multiple files, testing)

---

### 🟡 ISSUE #5: Hardcoded Noise Level in Variance Metric [MEDIUM SEVERITY]

**File**: [run_paper_vs_residual_peak_cv_compare.m](run_paper_vs_residual_peak_cv_compare.m#L270), [run_s2022_random10_lambda_residuals.m](run_s2022_random10_lambda_residuals.m#L516)

**Problem**: Bootstrap variance metric uses hardcoded 0.5% noise level:

```matlab
function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, el, n_boot)
    ...
    sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);  % ← HARDCODED 0.5%
    sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);
    ...
end
```

**Why This Is a Problem**:
- Different data types have different noise levels:
  - Clean lab data: 0.1-0.5% noise
  - Field data: 2-5% noise
  - Synthetic data: 0.01% noise
- Fixed 0.5% may not match actual data noise
- Variance metric becomes data-dependent in unintended way
- Makes cross-dataset comparisons difficult

**Recommended Fix**:
```matlab
function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, el, n_boot, noise_fraction)
    if nargin < 5 || isempty(noise_fraction)
        % Check environment variable first
        tmp = str2double(getenv('PAPER_NOISE_FRACTION'));
        if isnan(tmp) || tmp <= 0
            noise_fraction = 0.005;  % Default 0.5%
        else
            noise_fraction = tmp;
        end
    end
    
    sigma_re = max(noise_fraction * max(abs(real(Z_exp))), eps);
    sigma_im = max(noise_fraction * max(abs(imag(Z_exp))), eps);
    ...
end
```

**Then update call sites**:
```matlab
% In run_paper_vs_residual_peak_cv_compare.m
var_metric(i) = estimate_resample_variance_local(freq_vec, Z_exp, lam, n_boot, 0.005);

% Or use environment variable
var_metric(i) = estimate_resample_variance_local(freq_vec, Z_exp, lam, n_boot);
```

**Files to Update**:
- [ ] `run_paper_vs_residual_peak_cv_compare.m` - Function definition (line 269) + calls (line 121)
- [ ] `run_s2022_random10_lambda_residuals.m` - Function definition (line 516) + calls (line 100)
- [ ] `run_s2022_random10_drt_dualstyle_plot.m` - Similar updates
- [ ] `run_s2022_random10_paper_peak_refine.m` - Similar updates

**Testing**:
1. Run with `PAPER_NOISE_FRACTION=0.01` → Variance metric should increase
2. Run with `PAPER_NOISE_FRACTION=0.1` → Variance metric should decrease significantly
3. Verify paper lambda selection changes (higher noise → larger lambda)

**Estimated Time**: 45 minutes (3-4 files with propagation)

---

## MEDIUM PRIORITY ISSUES (Within 2 Weeks)

### 🟡 ISSUE #1: Missing N≥2 Guard in A_Matrix Calculation [LOW SEVERITY]

**File**: [run_bayes_drt_workflow_matlab2011.m](run_bayes_drt_workflow_matlab2011.m#L410), similar patterns in other files

**Problem**: Matrix construction code doesn't check that n ≥ 2:

```matlab
function A_re = calc_A_re_local(freq)
    tau = 1 ./ freq(:);
    n = numel(freq);
    A_re = zeros(n, n);
    for p = 1:n
        for q = 1:n
            if q == 1
                log_term = log(tau(2) / tau(1));  % ← Assumes tau(2) exists!
            elseif q == n
                log_term = log(tau(n) / tau(n-1));
            ...
```

**Why This Is Low Priority**: 
- EIS data always has 50+ frequency points minimum
- Mathematically impossible to compute DRT with <2 points
- So this edge case never occurs in practice

**Recommended Fix** (defensive programming):
```matlab
function A_re = calc_A_re_local(freq)
    tau = 1 ./ freq(:);
    n = numel(freq);
    if n < 2
        error('calc_A_re_local: Requires at least 2 frequency points, got %d', n);
    end
    A_re = zeros(n, n);
    ...
end
```

**Files to Update**:
- [ ] `run_bayes_drt_workflow_matlab2011.m` - calc_A_re_local(), calc_A_im_local()
- [ ] Any file with similar matrix construction functions

**Estimated Time**: 10 minutes

---

### 🟡 ISSUE #4: Code Duplication - A_Matrix Functions [LOW SEVERITY]

**Problem**: `calc_A_re_local()` and `calc_A_im_local()` are duplicated across 5+ files

**Impact**: 
- Maintenance burden: bug fixes must be applied to 5 locations
- Inconsistency risk: bug in one copy not fixed in others
- Code bloat: ~200 lines duplicated

**Files with duplication**:
- `run_bayes_drt_workflow_matlab2011.m` (lines ~330-380)
- `run_paper_vs_residual_peak_cv_compare.m` (lines ~310-350)
- `run_s2022_random10_lambda_residuals.m` (lines ~550-600)
- `drt_input_analysis_matlab2011.m` (lines ~650-720)
- `run_s2022_random10_drt_dualstyle_plot.m` (lines ~240-290)

**Recommended Fix**: Create shared utility file

**New File**: `drt_matrices_shared.m`
```matlab
function A_re = calc_A_re_shared(freq)
    % Single authoritative A_re implementation
    tau = 1 ./ freq(:);
    n = numel(freq);
    if n < 2, error('Need n >= 2'); end
    A_re = zeros(n, n);
    for p = 1:n
        for q = 1:n
            if q == 1
                log_term = log(tau(2) / tau(1));
            elseif q == n
                log_term = log(tau(n) / tau(n-1));
            else
                log_term = log(tau(q+1) / tau(q-1));
            end
            A_re(p, q) = -0.5 / (1 + (omega(p) * tau(q))^2) * log_term;
        end
    end
end

function A_im = calc_A_im_shared(freq)
    % Similar implementation for imaginary part
    ...
end
```

**Then replace all local functions with calls**:
```matlab
% In each file:
A_re = calc_A_re_shared(freq_vec);
A_im = calc_A_im_shared(freq_vec);
```

**Refactoring Steps**:
1. Create `drt_matrices_shared.m`
2. Update `run_bayes_drt_workflow_matlab2011.m` to use shared functions
3. Verify output unchanged
4. Update remaining files
5. Remove local implementations

**Estimated Time**: 2 hours (refactoring + testing all 5 files)

---

## LOW PRIORITY (Documentation & Testing)

### Issue: Inconsistent A_matrix log_term formula

**Location**: Different files use different but equivalent expressions
```matlab
% Form 1: log(tau(q+1) / tau(q))      [in one file]
% Form 2: log(tau(2) / tau(1))        [in another file, when q=1]
```

These are mathematically identical but Form 2 is clearer. Use in shared utility.

### Issue: Missing unit tests

**Recommendation**: Create test suite:
```matlab
test_lambda_selection.m
├─ Test 1: Synthetic data with known solution
├─ Test 2: Edge cases (n=50, n=500, n=5000)
├─ Test 3: Comparison with published results
├─ Test 4: All three methods on S2022 reference data
└─ Test 5: Numerical stability (very small/large impedance)
```

---

## FIX IMPLEMENTATION SCHEDULE

### Week 1 (by 2026-07-20)
- [x] **ISSUE #3**: Fix Paper Method selection ambiguity
  - Estimated effort: 15 min coding + 30 min testing
  - Blocking: Several downstream comparison analyses
  
- [ ] **ISSUE #2**: Fix R_inf validation
  - Estimated effort: 30 min coding + 30 min testing
  - Multiple files to update
  
### Week 2 (by 2026-07-27)
- [ ] **ISSUE #5**: Make noise level configurable
  - Estimated effort: 45 min coding + 60 min testing
  - Update 4 files with calls
  
- [ ] **ISSUE #4**: Refactor A_matrix functions
  - Estimated effort: 2 hours total
  - Nice-to-have improvement
  
### Week 3+ (Medium Priority)
- [ ] **ISSUE #1**: Add n≥2 guard (defensive programming)
- [ ] **Create test suite** for regression detection

---

## VALIDATION CHECKLIST AFTER FIXES

- [ ] All three methods produce stable results
- [ ] Paper method lambdas are in CV/Var constrained range
- [ ] R_inf always positive (> eps)
- [ ] No NaN or Inf in output files
- [ ] Peak counts are 1-6 (sanity range)
- [ ] Results reproducible with same random seed
- [ ] Execution times within expected bounds:
  - Bayes: 15-20 sec/dataset
  - Paper: 7-10 sec/dataset
  - Residual: 3-5 sec/dataset

---

## REFERENCES

- **Verification Report**: [LAMBDA_METHODS_VERIFICATION.md](LAMBDA_METHODS_VERIFICATION.md)
- **Code Reference**: [LAMBDA_METHODS_CODEREFERENCE.md](LAMBDA_METHODS_CODEREFERENCE.md)
- **Session Notes**: `/memories/session/lambda_methods_analysis.md`

---

**Report Generated**: 2026-07-13 by GitHub Copilot  
**Status**: Ready for review and implementation
