#' Functions to finalize data-specific parameter ranges
#'
#' These functions take a parameter object and modify the unknown parts of
#' `ranges` based on a data set and simple heuristics.
#'
#' @param object A `param` object or a list of `param` objects.
#'
#' @param x The predictor data. In some cases (see below) this should only
#' include numeric data.
#'
#' @param force A single logical that indicates that even if the parameter
#' object is complete, should it update the ranges anyway?
#'
#' @param log_vals A logical: should the ranges be set on the log10 scale?
#'
#' @param ... Other arguments to pass to the underlying parameter
#' finalizer functions. These must be empty for `get_rbf_range()`; its
#' arguments are named explicitly (`frac`, `scaled`, `n_pairs`,
#' `chunk_size`, `seed`).
#'
#' @param frac A double for the fraction of the data to be used for the upper
#' bound. For `get_n_frac_range()` and `get_batch_sizes()`, a vector of two
#' fractional values are required. For `get_rbf_range()`, it is the fraction of
#' rows that contribute distance pairs (as in [kernlab::sigest()]); the number
#' of pairs is capped at `n_pairs`.
#'
#' @param seed An integer to control the randomness of the calculations.
#'
#' @param scaled A single logical or a logical vector with one element per
#' column of `x` indicating which columns should be scaled to unit variance
#' before distances are computed. If any scaled column has zero variance,
#' scaling is disabled for every column and a warning is issued, matching the
#' `kernlab::sigest()` heuristic.
#'
#' @param n_pairs A positive integer: the sample budget, i.e. the maximum
#' number of randomly sampled pairs of rows whose squared distances are used
#' to estimate the quantiles. Capping the budget keeps both the random draws
#' and the resulting distance vector bounded regardless of the size of `x`.
#'
#' @param chunk_size A positive integer: the memory budget for distance
#' calculations. Distances are computed in blocks of at most this many pairs;
#' only the rows involved in one block are materialized. For sparse `x` the
#' full data set is never densified.
#'
#' @return
#'
#' An updated `param` object or a list of updated `param` objects depending
#' on what is provided in `object`.
#'
#' @details
#'
#' `finalize()` runs the embedded finalizer function contained in the `param`
#' object (`object$finalize`) and returns the updated version. The finalization
#' function is one of the `get_*()` helpers.
#'
#' The `get_*()` helper functions are designed to be used with the pipe
#' and update the parameter object in-place.
#'
#' `get_p()` and `get_log_p()` set the upper value of the range to be
#' the number of columns in the data (on the natural and
#' log10 scale, respectively).
#'
#' `get_n()` and `get_n_frac()` set the upper value to be the number of
#' rows in the data or a fraction of the total number of rows.
#'
#' `get_rbf_range()` sets both bounds based on the heuristic defined in
#' [kernlab::sigest()]. It requires that all columns in `x` be numeric. The
#' same 10th, 50th and 90th quantiles of pairwise squared distances are
#' estimated, but the computation is memory bounded: rows are sampled through
#' a local random stream, distances are computed in blocks no larger than
#' `chunk_size` pairs, and sparse matrices keep their sparse structure (only
#' the rows participating in a block are materialized). The number of sampled
#' pairs is capped at `n_pairs`, so neither the random draw nor the distance
#' vector grows with the size of `x`.
#'
#' Inputs follow explicit rules:
#'
#' * all columns must be numeric;
#' * all values must be finite, so missing (`NA`), `NaN` and infinite values
#'   raise an error instead of being silently dropped;
#' * columns with zero variance (constant columns) disable scaling for all
#'   columns, with a warning;
#' * data whose sampled rows are all identical leave no positive distances,
#'   which raises an error because the sigma quantiles cannot be estimated.
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' library(dplyr)
#' car_pred <- select(mtcars, -mpg)
#'
#' # Needs an upper bound
#' mtry()
#' finalize(mtry(), car_pred)
#'
#' # Nothing to do here since no unknowns
#' penalty()
#' finalize(penalty(), car_pred)
#'
#' library(kernlab)
#' library(tibble)
#' library(purrr)
#'
#' params <-
#'   tribble(
#'     ~parameter, ~object,
#'     "mtry", mtry(),
#'     "num_terms", num_terms(),
#'     "rbf_sigma", rbf_sigma()
#'   )
#' params
#'
#' # Note that `rbf_sigma()` has a default range that does not need to be
#' # finalized but will be changed if used in the function:
#' complete_params <-
#'   params |>
#'   mutate(object = map(object, finalize, car_pred))
#' complete_params
#'
#' params |>
#'   dplyr::filter(parameter == "rbf_sigma") |>
#'   pull(object)
#' complete_params |>
#'   dplyr::filter(parameter == "rbf_sigma") |>
#'   pull(object)
#'
#' @export
finalize <- function(object, ...) {
  UseMethod("finalize")
}

#' @export
#' @rdname finalize
finalize.list <- function(object, x, force = TRUE, ...) {
  check_bool(force)
  map(object, finalize, x, force, ...)
}

#' @export
#' @rdname finalize
finalize.param <- function(object, x, force = TRUE, ...) {
  check_bool(force)
  if (is.null(object$finalize)) {
    return(object)
  }
  if (!has_unknowns(object) & !force) {
    return(object)
  }
  object$finalize(object, x = x, ...)
}


safe_finalize <- function(object, x, force = TRUE, ...) {
  if (all(is.na(object))) {
    res <- NA
  } else {
    res <- finalize(object, x = x, force = force, ...)
  }
  res
}

#' @export
#' @rdname finalize
finalize.parameters <- function(object, x, force = TRUE, ...) {
  object$object <- map(object$object, safe_finalize, x = x, force = force, ...)
  object
}

# These two finalize methods are for cases when a tuning parameter has no
# parameter object or isn't listed in the tunable method.

#' @export
#' @rdname finalize
finalize.logical <- function(object, x, force = TRUE, ...) {
  object
}

#' @export
#' @rdname finalize
finalize.default <- function(object, x, force = TRUE, ...) {
  if (all(is.na(object))) {
    return(object)
  } else {
    cli::cli_abort(
      "Cannot finalize {.obj_type_friendly {object}}."
    )
  }
  object
}


#' @export
#' @rdname finalize
get_p <- function(object, x, log_vals = FALSE, ...) {
  check_inherits(object, "param")
  check_bool(log_vals)
  check_dots_empty()

  rngs <- range_get(object, original = FALSE)
  if (!is_unknown(rngs$upper)) {
    return(object)
  }

  x_dims <- dim(x)
  if (is.null(x_dims)) {
    cli::cli_abort(
      "Cannot determine number of columns. Is {.arg x} a 2D data object?"
    )
  }

  if (log_vals) {
    rngs[2] <- log10(x_dims[2])
  } else {
    rngs[2] <- x_dims[2]
  }

  if (object$type == "integer" & is.null(object$trans)) {
    rngs <- as.integer(rngs)
  }

  range_set(object, rngs)
}

#' @export
#' @rdname finalize
get_log_p <- function(object, x, ...) {
  check_inherits(object, "param")
  check_dots_empty()
  get_p(object, x, log_vals = TRUE)
}

#' @export
#' @rdname finalize
get_n_frac <- function(object, x, log_vals = FALSE, frac = 1 / 3, ...) {
  check_inherits(object, "param")
  check_bool(log_vals)
  check_number_decimal(frac, min = 0, max = 1)
  check_dots_empty()

  rngs <- range_get(object, original = FALSE)
  if (!is_unknown(rngs$upper)) {
    return(object)
  }

  x_dims <- dim(x)
  if (is.null(x_dims)) {
    cli::cli_abort(
      "Cannot determine number of columns. Is {.arg x} a 2D data object?"
    )
  }

  n_frac <- floor(x_dims[1] * frac)

  if (log_vals) {
    rngs[2] <- log10(n_frac)
  } else {
    rngs[2] <- n_frac
  }
  if (object$type == "integer" & is.null(object$trans) & !log_vals) {
    rngs <- as.integer(rngs)
  }
  range_set(object, rngs)
}

#' @export
#' @rdname finalize
get_n_frac_range <- function(
  object,
  x,
  log_vals = FALSE,
  frac = c(1 / 10, 5 / 10),
  ...
) {
  check_inherits(object, "param")
  check_bool(log_vals)
  check_frac_range(frac)
  check_dots_empty()

  rngs <- range_get(object, original = FALSE)
  if (!is_unknown(rngs$upper)) {
    return(object)
  }

  x_dims <- dim(x)
  if (is.null(x_dims)) {
    cli::cli_abort(
      "Cannot determine number of columns. Is {.arg x} a 2D data object?"
    )
  }

  n_frac <- sort(floor(x_dims[1] * frac))

  if (log_vals) {
    rngs <- log10(n_frac)
  } else {
    rngs <- n_frac
  }

  if (object$type == "integer" & is.null(object$trans) & !log_vals) {
    rngs <- as.integer(rngs)
  }

  range_set(object, rngs)
}

#' @export
#' @rdname finalize
get_n <- function(object, x, log_vals = FALSE, ...) {
  check_inherits(object, "param")
  check_dots_empty()
  get_n_frac(object, x, log_vals, frac = 1)
}

#' @export
#' @rdname finalize
get_rbf_range <- function(
  object,
  x,
  seed = sample.int(10^5, 1),
  frac = 0.5,
  scaled = TRUE,
  n_pairs = 1000L,
  chunk_size = 250L,
  ...
) {
  check_inherits(object, "param")
  check_number_decimal(frac, min = 0, max = 1)
  check_logical_vec(scaled)
  check_number_whole(n_pairs, min = 1)
  check_number_whole(chunk_size, min = 1)
  check_dots_empty()

  x <- coerce_rbf_input(x)

  with_seed(
    seed,
    rng <- estimate_rbf_sigma(
      x,
      frac = frac,
      scaled = scaled,
      n_pairs = as.integer(n_pairs),
      chunk_size = as.integer(chunk_size)
    )[-2L]
  )
  rng <- log10(rng)
  range_set(object, rng)
}

# Coerce and validate the `x` argument of `get_rbf_range()` without ever
# densifying a whole sparse matrix.
coerce_rbf_input <- function(x, call = caller_env()) {
  if (is.data.frame(x)) {
    non_numeric <- names(x)[!vapply(x, is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      cli::cli_abort(
        "The matrix version of the initialization data is not numeric.",
        call = call
      )
    }
    return(as.matrix(x))
  }

  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      cli::cli_abort(
        "The matrix version of the initialization data is not numeric.",
        call = call
      )
    }
    return(x)
  }

  if (is_sparse_matrix(x)) {
    if (!is_numeric_sparse(x)) {
      cli::cli_abort(
        "The matrix version of the initialization data is not numeric.",
        call = call
      )
    }
    if (!inherits(x, "CsparseMatrix")) {
      x <- rbf_as_csparse(x)
    }
    return(x)
  }

  cli::cli_abort(
    "Cannot determine number of columns. Is {.arg x} a 2D data object?",
    call = call
  )
}

is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix") && length(dim(x)) == 2L
}

# Internal seam for coercing alternative sparse representations (row-based,
# triplet) to the column-compressed form used by the estimator.
rbf_as_csparse <- function(x) {
  methods::as(x, "CsparseMatrix")
}

# Pattern (`nsparseMatrix`) and logical (`lsparseMatrix`) sparse matrices do
# not store numbers; only numeric sparse matrices (double or integer) are
# usable.
is_numeric_sparse <- function(x) {
  inherits(x, "dsparseMatrix") || inherits(x, "isparseMatrix")
}

check_logical_vec <- function(x, ..., arg = caller_arg(x),
                              call = caller_env()) {
  if (!is.logical(x) || anyNA(x) || length(x) < 1L) {
    stop_input_type(
      x,
      "a non-empty logical vector without missing values",
      arg = arg,
      call = call
    )
  }
  invisible(x)
}

# Pure-R, memory-bounded re-implementation of the `kernlab::sigest()`
# heuristic. The full data set is never densified: two index vectors of at
# most `n_pairs` rows are drawn from a local random stream, and squared
# distances are evaluated block by block, materializing only the rows of one
# block at a time.
estimate_rbf_sigma <- function(x, frac, scaled, n_pairs, chunk_size,
                               call = caller_env()) {
  m <- nrow(x)
  p <- ncol(x)
  if (is.null(m) || is.null(p) || m < 1L || p < 1L) {
    cli::cli_abort(
      "The initialization data must have at least one row and one column.",
      call = call
    )
  }

  check_finite_values(x, call = call)

  if (length(scaled) == 1L) {
    scaled <- rep(scaled, p)
  } else if (length(scaled) != p) {
    cli::cli_abort(
      c(
        "{.arg scaled} must have length 1 or one value per column.",
        i = "The data has {p} column{?s}; {.arg scaled} has length
             {length(scaled)}."
      ),
      call = call
    )
  }

  # Column standard deviations used to scale distances below. `sds` stays
  # length `p`, so no copy of `x` is made for scaling.
  sds <- rep(1, p)
  if (any(scaled)) {
    vars <- rbf_col_vars(x[, scaled, drop = FALSE])
    if (any(!is.finite(vars)) || any(vars == 0)) {
      n_const <- sum(!is.finite(vars) | vars == 0)
      cli::cli_warn(
        c(
          "{n_const} variable{?s} {?has/have} zero variance; scaling is
           disabled for all columns.",
          i = "This matches the {.fun kernlab::sigest} fallback behavior."
        ),
        call = call
      )
      scaled[] <- FALSE
    } else {
      sds[scaled] <- sqrt(vars)
    }
  }

  n <- floor(frac * m)
  if (n < 1L) {
    cli::cli_abort(
      c(
        "{.arg frac} = {frac} samples fewer than one pair from {m} row{?s}.",
        i = "Increase {.arg frac} or supply more rows."
      ),
      call = call
    )
  }
  n <- min(n, n_pairs)

  # Two independent samples with replacement, as in `kernlab::sigest()`.
  index <- sample(seq_len(m), n, replace = TRUE)
  index2 <- sample(seq_len(m), n, replace = TRUE)

  dist <- numeric(n)
  for (start in seq(1L, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n)
    rows_1 <- index[start:end]
    rows_2 <- index2[start:end]
    block <- x[rows_1, , drop = FALSE] - x[rows_2, , drop = FALSE]
    # Distance on scaled columns: divide each column by its standard
    # deviation. `block` is the only materialized, possibly dense, object.
    block <- rbf_scale_block(block, sds)
    dist[start:end] <- as.numeric(Matrix::rowSums(block^2))
  }

  nonzero <- dist[dist != 0]
  if (length(nonzero) == 0L) {
    cli::cli_abort(
      "All sampled pairs have zero distance, so the sigma quantiles cannot be
       estimated. Do the data contain any variation between rows?",
      call = call
    )
  }

  1 / stats::quantile(nonzero, probs = c(0.9, 0.5, 0.1))
}

check_finite_values <- function(x, call = caller_env()) {
  if (is_sparse_matrix(x)) {
    # Structural zeros are finite by definition; only explicit values need
    # checking, and the full matrix is not materialized.
    bad <- !is.finite(x@x)
  } else {
    bad <- !is.finite(x)
  }
  if (any(bad)) {
    cli::cli_abort(
      c(
        "The initialization data must contain only finite numeric values.",
        i = "Remove or impute missing ({.code NA}/{.code NaN}) or infinite
             values before finalizing {.arg rbf_sigma}."
      ),
      call = call
    )
  }
  invisible(NULL)
}

# Column variances computed without densifying sparse input.
rbf_col_vars <- function(x) {
  if (is_sparse_matrix(x)) {
    n <- nrow(x)
    p <- ncol(x)
    col_sums <- Matrix::colSums(x)
    nnz <- diff(x@p)
    # Two-pass formula on explicit (nonzero) entries; structural zeros
    # contribute their squared mean through the `(n - nnz)` term.
    full <- numeric(p)
    if (length(x@x) > 0L) {
      grp <- rep.int(seq_len(p), nnz)
      sums_sq <- rowsum(
        (x@x - col_sums[grp] / n)^2,
        grp,
        reorder = FALSE
      )
      full[as.integer(rownames(sums_sq))] <- sums_sq
    }
    (full + (n - nnz) * (col_sums / n)^2) / (n - 1)
  } else {
    apply(x, 2L, stats::var)
  }
}

# Divide every column of a block (dense matrix or sparse Matrix) by the
# column standard deviation, keeping block-level sparsity when possible.
rbf_scale_block <- function(block, sds) {
  if (is_sparse_matrix(block)) {
    # A sparse diagonal multiplier keeps the block sparse and costs O(p)
    # memory rather than the O(p^2) of the base `diag()` matrix.
    block %*% Matrix::Diagonal(x = 1 / sds)
  } else {
    block / rep(sds, each = nrow(block))
  }
}

#' Get batch sizes
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @inheritParams finalize
#' @keywords internal
#' @export
get_batch_sizes <- function(object, x, frac = c(1 / 10, 1 / 3), ...) {
  lifecycle::deprecate_warn("1.4.2", "get_batch_sizes()")
  check_dots_empty()

  rngs <- range_get(object, original = FALSE)
  if (!is_unknown(rngs$lower) & !is_unknown(rngs$upper)) {
    return(object)
  }

  x_dims <- dim(x)
  if (is.null(x_dims)) {
    cli::cli_abort(
      "Cannot determine number of columns. Is {.arg x} a 2D data object?"
    )
  }

  n_frac <- sort(floor(x_dims[1] * frac))
  n_frac <- log2(n_frac)

  if (object$type == "integer" & is.null(object$trans)) {
    n_frac <- as.integer(n_frac)
  }

  range_set(object, n_frac)
}
