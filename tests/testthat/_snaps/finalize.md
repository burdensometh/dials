# finalize validates inputs

    Code
      finalize(penalty(), mtcars, force = "yes")
    Condition
      Error in `finalize()`:
      ! `force` must be `TRUE` or `FALSE`, not the string "yes".

---

    Code
      get_p(penalty(), mtcars, log_vals = "yes")
    Condition
      Error in `get_p()`:
      ! `log_vals` must be `TRUE` or `FALSE`, not the string "yes".

---

    Code
      get_n_frac(mtry(), mtcars, frac = "half")
    Condition
      Error in `get_n_frac()`:
      ! `frac` must be a number, not the string "half".

---

    Code
      get_n_frac(mtry(), mtcars, frac = 1.5)
    Condition
      Error in `get_n_frac()`:
      ! `frac` must be a number between 0 and 1, not the number 1.5.

---

    Code
      get_n_frac_range(mtry(), mtcars, frac = 0.5)
    Condition
      Error in `get_n_frac_range()`:
      ! `frac` must be a numeric vector of length 2 with values between 0 and 1, not the number 0.5.

---

    Code
      get_n_frac_range(mtry(), mtcars, frac = c(0.1, 1.5))
    Condition
      Error in `get_n_frac_range()`:
      ! `frac` must be a numeric vector of length 2 with values between 0 and 1, not a double vector.

---

    Code
      get_rbf_range("not a param", mtcars)
    Condition
      Error in `get_rbf_range()`:
      ! `object` must be a <param> object, not the string "not a param".

# estimate columns

    Code
      get_p(1:10)
    Condition
      Error in `get_p()`:
      ! `object` must be a <param> object, not an integer vector.

---

    Code
      get_p(1:10, 1:10)
    Condition
      Error in `get_p()`:
      ! `object` must be a <param> object, not an integer vector.

---

    Code
      get_p(mtry(), 1:10)
    Condition
      Error in `get_p()`:
      ! Cannot determine number of columns. Is `x` a 2D data object?

# estimate rows

    Code
      get_n(1:10)
    Condition
      Error in `get_n()`:
      ! `object` must be a <param> object, not an integer vector.

---

    Code
      get_n(1:10, 1:10)
    Condition
      Error in `get_n()`:
      ! `object` must be a <param> object, not an integer vector.

---

    Code
      get_n(mtry(), 1:10)
    Condition
      Error in `get_n_frac()`:
      ! Cannot determine number of columns. Is `x` a 2D data object?

# `get_batch_size() is deprecated

    Code
      bsizes <- get_batch_sizes(batch_size(), iris, frac = c(0.3, 0.7))
    Condition
      Warning:
      `get_batch_sizes()` was deprecated in dials 1.4.2.

# estimate sigma

    Code
      get_rbf_range(rbf_sigma(), iris)
    Condition
      Error in `get_rbf_range()`:
      ! The matrix version of the initialization data is not numeric.

# sigma sample and memory budgets are respected

    Code
      get_rbf_range(rbf_sigma(), x, n_pairs = 0)
    Condition
      Error in `get_rbf_range()`:
      ! `n_pairs` must be a whole number larger than or equal to 1, not the number 0.

---

    Code
      get_rbf_range(rbf_sigma(), x, chunk_size = 0)
    Condition
      Error in `get_rbf_range()`:
      ! `chunk_size` must be a whole number larger than or equal to 1, not the number 0.

---

    Code
      get_rbf_range(rbf_sigma(), x, frac = 0)
    Condition
      Error in `get_rbf_range()`:
      ! `frac` = 0 samples fewer than one pair from 400 rows.
      i Increase `frac` or supply more rows.

# sigma estimation has explicit rules for degenerate inputs

    Code
      get_rbf_range(rbf_sigma(), na_x, seed = 1L)
    Condition
      Error in `get_rbf_range()`:
      ! The initialization data must contain only finite numeric values.
      i Remove or impute missing (`NA`/`NaN`) or infinite values before finalizing `rbf_sigma`.

---

    Code
      get_rbf_range(rbf_sigma(), constant, seed = 1L)
    Condition
      Warning in `get_rbf_range()`:
      1 variable has zero variance; scaling is disabled for all columns.
      i This matches the `kernlab::sigest()` fallback behavior.
    Message
      Radial Basis Function sigma (quantitative)
      Transformer: log-10 [1e-100, Inf]
      Range (transformed scale): [-1.28, -0.873]

---

    Code
      get_rbf_range(rbf_sigma(), matrix(1, nrow = 20, ncol = 3), seed = 1L)
    Condition
      Warning in `get_rbf_range()`:
      3 variables have zero variance; scaling is disabled for all columns.
      i This matches the `kernlab::sigest()` fallback behavior.
      Error in `get_rbf_range()`:
      ! All sampled pairs have zero distance, so the sigma quantiles cannot be estimated. Do the data contain any variation between rows?

---

    Code
      get_rbf_range(rbf_sigma(), iris)
    Condition
      Error in `get_rbf_range()`:
      ! The matrix version of the initialization data is not numeric.

---

    Code
      get_rbf_range(rbf_sigma(), 1:10)
    Condition
      Error in `get_rbf_range()`:
      ! Cannot determine number of columns. Is `x` a 2D data object?

# sigma argument validation is explicit

    Code
      get_rbf_range(rbf_sigma(), x, scaled = c(TRUE, FALSE))
    Condition
      Error in `get_rbf_range()`:
      ! `scaled` must have length 1 or one value per column.
      i The data has 5 columns; `scaled` has length 2.

---

    Code
      get_rbf_range(rbf_sigma(), x, scaled = NA)
    Condition
      Error in `get_rbf_range()`:
      ! `scaled` must be a non-empty logical vector without missing values, not `NA`.

---

    Code
      get_rbf_range(rbf_sigma(), x, bogus = 1)
    Condition
      Error in `get_rbf_range()`:
      ! `...` must be empty.
      x Problematic argument:
      * bogus = 1

# finalize interfaces

    Code
      finalize("threshold", mtcars)
    Condition
      Error in `finalize()`:
      ! Cannot finalize a string.

---

    Code
      finalize(list(TRUE, 3), mtcars)
    Condition
      Error in `map()`:
      i In index: 2.
      Caused by error in `.f()`:
      ! Cannot finalize a number.

