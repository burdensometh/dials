test_that("finalize validates inputs", {
  expect_snapshot(error = TRUE, finalize(penalty(), mtcars, force = "yes"))
  expect_snapshot(error = TRUE, get_p(penalty(), mtcars, log_vals = "yes"))
  expect_snapshot(error = TRUE, get_n_frac(mtry(), mtcars, frac = "half"))
  expect_snapshot(error = TRUE, get_n_frac(mtry(), mtcars, frac = 1.5))
  expect_snapshot(error = TRUE, get_n_frac_range(mtry(), mtcars, frac = 0.5))
  expect_snapshot(
    error = TRUE,
    get_n_frac_range(mtry(), mtcars, frac = c(0.1, 1.5))
  )
  expect_snapshot(error = TRUE, get_rbf_range("not a param", mtcars))
})

test_that("estimate columns", {
  expect_snapshot(error = TRUE, get_p(1:10))
  expect_snapshot(error = TRUE, get_p(1:10, 1:10))
  expect_snapshot(error = TRUE, get_p(mtry(), 1:10))

  expect_equal(
    range_get(get_p(mtry(), mtcars)),
    list(lower = 1, upper = ncol(mtcars))
  )
  expect_equal(
    range_get(get_log_p(mtry_long(), mtcars), original = FALSE),
    list(lower = 0, upper = log10(ncol(mtcars)))
  )
})


test_that("estimate rows", {
  expect_snapshot(error = TRUE, get_n(1:10))
  expect_snapshot(error = TRUE, get_n(1:10, 1:10))
  expect_snapshot(error = TRUE, get_n(mtry(), 1:10))

  expect_equal(
    range_get(
      get_n_frac(mtry_long(), mtcars, log_vals = TRUE),
      original = FALSE
    ),
    list(lower = 0, upper = 1)
  )
  expect_equal(
    get_n_frac_range(
      sample_size(c(unknown(), unknown())),
      warpbreaks,
      frac = c(.3, .7)
    )$range,
    list(lower = 16, upper = 37)
  )
  expect_equal(
    get_n_frac_range(
      sample_size(c(unknown(), unknown()), trans = transform_log10()),
      x = warpbreaks,
      log_vals = TRUE,
      frac = c(.3, .7)
    )$range,
    list(lower = log10(16), upper = log10(37))
  )

  expect_equal(
    get_n_frac_range(
      sample_size(c(16L, unknown())),
      x = warpbreaks,
      log_vals = FALSE,
      frac = c(.3, .7)
    )$range,
    list(lower = 16, upper = 37)
  )

  expect_equal(get_n_frac(mtry(c(1, 2)), mtcars), mtry(c(1, 2)))
})

test_that("`get_batch_size() is deprecated", {
  expect_snapshot(
    bsizes <- get_batch_sizes(batch_size(), iris, frac = c(0.3, 0.7))
  )
})

test_that("`get_batch_size() works", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mock_batch_size_with_unknown <- new_quant_param(
    type = "integer",
    range = c(unknown(), unknown()),
    inclusive = c(TRUE, TRUE),
    trans = transform_log2(),
    label = c(batch_size = "Batch Size"),
    finalize = get_batch_sizes
  )
  expect_equal(
    get_batch_sizes(
      mock_batch_size_with_unknown,
      iris,
      frac = c(.3, .7)
    )$range,
    list(lower = log2(45), upper = log2(105))
  )
})


test_that("estimate sigma", {
  skip_if_not_installed("kernlab")
  suppressMessages(library(kernlab))

  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), iris))

  run_1 <- range_get(get_rbf_range(rbf_sigma(), mtcars, seed = 5624))
  run_2 <- range_get(get_rbf_range(rbf_sigma(), mtcars, seed = 5624))

  expect_equal(run_1, run_2)
})

test_that("sigma estimate matches the kernlab::sigest() heuristic", {
  skip_if_not_installed("kernlab")
  suppressMessages(library(kernlab))

  # The memory-bounded implementation estimates the same 10/90 percentiles of
  # pairwise squared distance as sigest() when the sample budget is not hit.
  for (seed in c(5624L, 1L, 42L, 7L)) {
    new <- range_get(
      get_rbf_range(rbf_sigma(), mtcars, seed = seed, n_pairs = 1000L),
      original = FALSE
    )
    old <- log10(
      withr::with_seed(
        seed,
        kernlab::sigest(as.matrix(mtcars))
      )[-2]
    )
    expect_equal(c(new$lower, new$upper), unname(old), ignore_attr = TRUE)
  }
})

test_that("sigma estimate works on sparse matrices without densifying", {
  skip_if_not_installed("Matrix")

  set.seed(2)
  dense <- matrix(rnorm(4000), nrow = 400)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)

  dense_range <- range_get(
    get_rbf_range(rbf_sigma(), dense, seed = 11L),
    original = FALSE
  )
  sparse_range <- range_get(
    get_rbf_range(rbf_sigma(), sparse, seed = 11L),
    original = FALSE
  )
  expect_equal(sparse_range, dense_range, tolerance = 1e-8)

  # Other sparse representations are coerced before estimation.
  r_sparse <- as(sparse, "RsparseMatrix")
  r_range <- range_get(
    get_rbf_range(rbf_sigma(), r_sparse, seed = 11L),
    original = FALSE
  )
  expect_equal(r_range, sparse_range, tolerance = 1e-8)

  # The full data set must never be densified. Two guards:
  # (a) a Csparse input never triggers the representation-coercion seam;
  big <- Matrix::rsparsematrix(50000, 8000, density = 0.0005)
  local_mocked_bindings(
    rbf_as_csparse = function(x) {
      stop("column-compressed input must not be coerced again")
    }
  )
  expect_no_error(
    get_rbf_range(
      rbf_sigma(), big,
      seed = 3L, n_pairs = 400L, chunk_size = 100L
    )
  )

  # (b) peak V-cell growth stays a small fraction of the dense size
  # (400 million doubles here), since only chunk-sized row slices exist.
  rm(big)
  big <- Matrix::rsparsematrix(50000, 8000, density = 0.0005)
  vcells <- function() gc()["Vcells", 2L]
  gc(full = TRUE)
  before <- vcells()
  get_rbf_range(
    rbf_sigma(), big,
    seed = 3L, n_pairs = 400L, chunk_size = 100L
  )
  growth_mb <- vcells() - before
  expect_lt(growth_mb, 50)
})

test_that("sigma sample and memory budgets are respected", {
  set.seed(4)
  x <- matrix(rnorm(4000), nrow = 400)

  # The pair count is min(floor(frac * nrow(x)), n_pairs).
  capped <- get_rbf_range(
    rbf_sigma(), x, seed = 5L, n_pairs = 25L, chunk_size = 10L
  )
  expect_s3_class(capped, "quant_param")

  # Different budgets are allowed to give different ranges, but both are
  # finite and ordered.
  for (budget in c(10L, 100L, 1000L)) {
    r <- range_get(
      get_rbf_range(rbf_sigma(), x, seed = 5L, n_pairs = budget),
      original = FALSE
    )
    expect_true(is.finite(r$lower) && is.finite(r$upper))
    expect_lt(r$lower, r$upper)
  }

  # Chunk sizes smaller than the pair budget still use every sampled pair.
  full <- range_get(
    get_rbf_range(rbf_sigma(), x, seed = 6L, chunk_size = 7L),
    original = FALSE
  )
  one <- range_get(
    get_rbf_range(rbf_sigma(), x, seed = 6L, chunk_size = 1000L),
    original = FALSE
  )
  expect_equal(full, one)

  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), x, n_pairs = 0))
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), x, chunk_size = 0))
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), x, frac = 0))
})

test_that("sigma estimation has explicit rules for degenerate inputs", {
  skip_if_not_installed("Matrix")

  set.seed(8)
  x <- matrix(rnorm(100), nrow = 20)

  # Non-finite values are rejected (densely and sparsely) rather than dropped.
  na_x <- x
  na_x[1, 1] <- NA
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), na_x, seed = 1L))

  inf_x <- x
  inf_x[2, 2] <- Inf
  expect_error(get_rbf_range(rbf_sigma(), inf_x, seed = 1L), "finite")

  sparse <- Matrix::Matrix(x, sparse = TRUE)
  na_sparse <- sparse
  na_sparse@x[1] <- NA
  expect_error(
    get_rbf_range(rbf_sigma(), na_sparse, seed = 1L),
    "finite"
  )

  # Constant columns disable scaling with a warning, matching sigest().
  constant <- cbind(x, const = 3)
  expect_snapshot(
    get_rbf_range(rbf_sigma(), constant, seed = 1L)
  )
  dense_c <- suppressWarnings(
    range_get(get_rbf_range(rbf_sigma(), constant, seed = 1L), original = FALSE)
  )
  sparse_c <- suppressWarnings(
    range_get(
      get_rbf_range(rbf_sigma(), Matrix::Matrix(constant, sparse = TRUE),
                    seed = 1L),
      original = FALSE
    )
  )
  expect_equal(sparse_c, dense_c, tolerance = 1e-8)

  # Data with no variation between rows leaves no positive distances.
  expect_snapshot(
    error = TRUE,
    get_rbf_range(rbf_sigma(), matrix(1, nrow = 20, ncol = 3), seed = 1L)
  )

  # Non-numeric and non-2D inputs follow the same rules as before.
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), iris))
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), 1:10))
  pattern <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L), j = c(1L, 2L, 3L), dims = c(4L, 4L)
  )
  expect_error(get_rbf_range(rbf_sigma(), pattern, seed = 1L), "numeric")
})

test_that("sigma argument validation is explicit", {
  x <- matrix(rnorm(100), nrow = 20)

  expect_snapshot(
    error = TRUE,
    get_rbf_range(rbf_sigma(), x, scaled = c(TRUE, FALSE))
  )
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), x, scaled = NA))
  expect_snapshot(error = TRUE, get_rbf_range(rbf_sigma(), x, bogus = 1))
})

test_that("estimated sigma range enters the native grid", {
  finalized <- get_rbf_range(rbf_sigma(), mtcars, seed = 5624)
  bounds <- range_get(finalized)
  grid <- grid_regular(finalized, levels = 5)

  expect_true(all(grid$rbf_sigma >= bounds$lower))
  expect_true(all(grid$rbf_sigma <= bounds$upper))

  grid_random <- grid_random(finalized, size = 7)
  expect_true(all(grid_random$rbf_sigma >= bounds$lower))
  expect_true(all(grid_random$rbf_sigma <= bounds$upper))
})


test_that("force", {
  rbf_sigma_final <- finalize(rbf_sigma(), mtcars)
  rbf_sigma_same <- finalize(rbf_sigma(), mtcars, force = FALSE)

  expect_false(rbf_sigma_final$range$lower == rbf_sigma()$range$lower)
  expect_false(rbf_sigma_final$range$upper == rbf_sigma()$range$upper)
  expect_true(rbf_sigma_same$range$lower == rbf_sigma()$range$lower)
  expect_true(rbf_sigma_same$range$upper == rbf_sigma()$range$upper)
})

test_that("finalize interfaces", {
  ## list interface
  set.seed(1)
  rbf_sigma_final <- finalize(rbf_sigma(), mtcars)
  set.seed(1)
  rbf_sigma_list <- finalize(list(rbf = rbf_sigma()), mtcars)
  expect_true(is.list(rbf_sigma_list))
  expect_equal(rbf_sigma_list$rbf, rbf_sigma_final)

  ## parameters interface
  set.seed(1)
  param <- parameters(list(rbf = rbf_sigma(), thresh = threshold()))
  param_final <- finalize(param, mtcars)
  expect_equal(
    param_final$object[[which(param_final$id == "rbf")]],
    rbf_sigma_final
  )
  expect_equal(
    param_final$object[[which(param_final$id == "thresh")]],
    threshold()
  )

  ## no parameter object
  param_na <- param
  param_na$object[1] <- NA
  expect_no_error(param_na_final <- finalize(param_na, mtcars))
  expect_true(is.na(param_na_final$object[[1]]))

  ## nothing occurs
  expect_equal(finalize(threshold(), mtcars), threshold())

  ## logical
  expect_true(finalize(TRUE, mtcars))
  expect_true(is.na(finalize(NA, mtcars))) # <- a logical

  ## NA and other objects
  expect_true(is.na(finalize(NA_character_, mtcars)))
  expect_snapshot(finalize("threshold", mtcars), error = TRUE)
  expect_snapshot(finalize(list(TRUE, 3), mtcars), error = TRUE)
})
