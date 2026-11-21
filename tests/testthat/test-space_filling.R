test_that("`grid_max_entropy()` is deprecated", {
  set.seed(1)
  expect_snapshot(
    grid_max_entropy(mixture(), trees(), size = 2)
  )
})

test_that("max entropy designs", {
  withr::local_options(lifecycle_verbosity = "quiet")

  grid_1 <- grid_max_entropy(
    cost(),
    mixture(),
    size = 11,
    original = FALSE
  )
  expect_equal(nrow(grid_1), 11L)
  expect_true(all(grid_1$mixture > 0 & grid_1$mixture < 1))
  expect_true(all(grid_1$cost > -10 & grid_1$cost < 5))

  grid_2 <- grid_max_entropy(
    cost(),
    mixture(),
    size = 11,
    original = TRUE
  )
  expect_true(all(grid_2$cost > 2^-10 & grid_2$cost < 2^5))

  grid_3 <- grid_max_entropy(
    cost(),
    size = 11,
    original = FALSE
  )
  expect_equal(ncol(grid_3), 1L)

  grid_4 <- grid_max_entropy(
    list(cost = cost(), mix = mixture()),
    size = 11,
    original = FALSE
  )
  expect_equal(nrow(grid_4), 11L)
  expect_true(all(grid_4$mix > 0 & grid_4$mix < 1))
  expect_true(all(grid_4$cost > -10 & grid_4$cost < 5))

  expect_snapshot(
    error = TRUE,
    grid_max_entropy(
      mtry(),
      size = 11,
      original = FALSE
    )
  )

  expect_snapshot(
    error = TRUE,
    grid_max_entropy(
      mixture(),
      levels = 11
    )
  )
})

test_that("`grid_latin_hypercube()` is deprecated", {
  set.seed(1)
  expect_snapshot(
    grid_latin_hypercube(mixture(), trees(), size = 2)
  )
})

test_that("latin square designs", {
  withr::local_options(lifecycle_verbosity = "quiet")

  grid_1 <- grid_latin_hypercube(
    cost(),
    mixture(),
    size = 11,
    original = FALSE
  )
  expect_equal(nrow(grid_1), 11L)
  expect_true(all(grid_1$mixture > 0 & grid_1$mixture < 1))
  expect_true(all(grid_1$cost > -10 & grid_1$cost < 5))

  grid_2 <- grid_latin_hypercube(
    cost(),
    mixture(),
    size = 11,
    original = TRUE
  )
  expect_true(all(grid_2$cost > 2^-10 & grid_2$cost < 2^5))

  grid_3 <- grid_latin_hypercube(
    cost(),
    size = 11,
    original = FALSE
  )
  expect_equal(ncol(grid_3), 1L)

  grid_4 <- grid_latin_hypercube(
    list(cost = cost(), mix = mixture()),
    size = 11,
    original = FALSE
  )
  expect_equal(nrow(grid_4), 11L)
  expect_true(all(grid_4$mix > 0 & grid_4$mix < 1))
  expect_true(all(grid_4$cost > -10 & grid_4$cost < 5))

  expect_lt(
    nrow(grid_latin_hypercube(prod_degree(), prune_method(), size = 20)),
    20
  )

  expect_snapshot(
    error = TRUE,
    grid_latin_hypercube(
      mtry(),
      size = 11,
      original = FALSE
    )
  )

  expect_snapshot(
    error = TRUE,
    grid_latin_hypercube(
      mixture(),
      levels = 11
    )
  )
})


test_that("sfd package designs - default", {
  size <- 11
  prm <- parameters(mom = momentum(), mixture())
  vls <- prm$object |> purrr::map(\(.x) dials::value_seq(.x, size))

  dials_2_any <- grid_space_filling(prm, size = size, type = "any")
  sfd_2_any <- sfd::get_design(2, num_points = size, type = "any")
  names(sfd_2_any) <- prm$id
  for (i in 1:2) {
    sfd_2_any[[i]] <- vls[[i]][sfd_2_any[[i]]]
  }

  expect_equal(dials_2_any, sfd_2_any)
})

test_that("sfd package designs AE", {
  size <- 11
  prm <- parameters(mom = momentum(), mixture())
  vls <- prm$object |> purrr::map(\(.x) dials::value_seq(.x, size))

  dials_2_any <- grid_space_filling(prm, size = size, type = "any")
  sfd_2_any <- sfd::get_design(2, num_points = size, type = "any")

  dials_2_ae <- grid_space_filling(prm, size = size, type = "audze_eglais")
  sfd_2_ae <- sfd::get_design(2, num_points = size, type = "audze_eglais")
  names(sfd_2_ae) <- prm$id
  for (i in 1:2) {
    sfd_2_ae[[i]] <- vls[[i]][sfd_2_ae[[i]]]
  }

  expect_equal(dials_2_ae, sfd_2_ae)
  expect_equal(dials_2_any, sfd_2_ae)
})

test_that("sfd package designs - MaxMin L1", {
  size <- 11
  prm <- parameters(mom = momentum(), mixture())
  vls <- prm$object |> purrr::map(\(.x) dials::value_seq(.x, size))

  dials_2_mml1 <- grid_space_filling(prm, size = size, type = "max_min_l1")
  sfd_2_mml1 <- sfd::get_design(2, num_points = size, type = "max_min_l1")
  names(sfd_2_mml1) <- prm$id
  for (i in 1:2) {
    sfd_2_mml1[[i]] <- vls[[i]][sfd_2_mml1[[i]]]
  }

  expect_equal(dials_2_mml1, sfd_2_mml1)
})

test_that("sfd package designs - MaxMin L2", {
  size <- 11
  prm <- parameters(mom = momentum(), mixture())
  vls <- prm$object |> purrr::map(\(.x) dials::value_seq(.x, size))

  dials_2_mml2 <- grid_space_filling(prm, size = size, type = "max_min_l2")
  sfd_2_mml2 <- sfd::get_design(2, num_points = size, type = "max_min_l2")
  names(sfd_2_mml2) <- prm$id
  for (i in 1:2) {
    sfd_2_mml2[[i]] <- vls[[i]][sfd_2_mml2[[i]]]
  }

  expect_equal(dials_2_mml2, sfd_2_mml2)
})

test_that("sfd package designs - uniform", {
  size <- 11
  prm <- parameters(mom = momentum(), mixture())
  vls <- prm$object |> purrr::map(\(.x) dials::value_seq(.x, size))

  dials_2_unif <- grid_space_filling(prm, size = size, type = "uniform")
  sfd_2_unif <- sfd::get_design(2, num_points = size, type = "uniform")
  names(sfd_2_unif) <- prm$id
  for (i in 1:2) {
    sfd_2_unif[[i]] <- vls[[i]][sfd_2_unif[[i]]]
  }

  expect_equal(dials_2_unif, sfd_2_unif)
})

test_that("DiceDesign package designs - max entropy", {
  withr::local_options(lifecycle_verbosity = "quiet")

  size <- 11
  prm <- parameters(mom = momentum(), mixture())

  set.seed(1)
  dials_2_maxent <- grid_space_filling(prm, size = size, type = "max_entropy")

  set.seed(1)
  dials_2_exp <- grid_max_entropy(prm, size = size)

  expect_equal(dials_2_maxent, dials_2_exp)

  ###

  set.seed(1)
  dials_2_maxent <- grid_space_filling(
    prm,
    size = size,
    type = "max_entropy",
    iter = 1
  )

  set.seed(1)
  dials_2_exp <- grid_max_entropy(prm, size = size, iter = 1)

  expect_equal(dials_2_maxent, dials_2_exp)

  ###

  set.seed(1)
  dials_2_maxent <-
    grid_space_filling(
      prm,
      size = size,
      type = "max_entropy",
      variogram_range = 0.1
    )

  set.seed(1)
  dials_2_exp <- grid_max_entropy(prm, size = size, variogram_range = 0.1)

  expect_equal(dials_2_maxent, dials_2_exp)
})

test_that("DiceDesign package designs - latin hypercube", {
  withr::local_options(lifecycle_verbosity = "quiet")

  size <- 11
  prm <- parameters(mom = momentum(), mixture())

  set.seed(1)
  dials_2_lh <- grid_space_filling(prm, size = size, type = "latin_hypercube")

  set.seed(1)
  dials_2_exp <- grid_latin_hypercube(prm, size = size)

  expect_equal(dials_2_lh, dials_2_exp)
})


test_that("no pre-made design", {
  size <- 501
  prm <- parameters(mom = momentum(), mixture())

  set.seed(1)
  dials_2_big <- grid_space_filling(prm, size = size, type = "any", iter = 2)

  set.seed(1)
  dials_2_exp <- grid_space_filling(
    prm,
    size = size,
    type = "max_entropy",
    iter = 2
  )

  expect_equal(dials_2_big, dials_2_exp)
})


test_that("S3 methods for space-filling", {
  size <- 12
  prm <- parameters(mixture(), mom = momentum(), activation(c("relu", "tanh")))

  design_paramset <- grid_space_filling(prm, size = size, type = "uniform")
  design_dots <-
    grid_space_filling(
      mixture(),
      mom = momentum(),
      activation(c("relu", "tanh")),
      size = size,
      type = "uniform"
    )
  expect_equal(design_paramset, design_dots)

  ###

  design_list <-
    grid_space_filling(
      list(
        mixture(),
        mom = momentum(),
        activation(c("relu", "tanh"))
      ),
      size = size,
      type = "uniform"
    )
  expect_equal(design_paramset, design_list)

  ## also:
  expect_snapshot(error = TRUE, {
    grid_space_filling(prm, levels = size, type = "uniform")
  })
})

test_that("1-point grid", {
  expect_silent({
    set.seed(1)
    grid <- grid_space_filling(parameters(neighbors()), size = 1)
  })
  expect_equal(nrow(grid), 1L)
})

test_that("grid_space_filling validates inputs", {
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), size = "five"))
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), size = -1))
  expect_snapshot(
    error = TRUE,
    grid_space_filling(penalty(), variogram_range = -1)
  )
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), iter = "many"))
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), original = "yes"))
})

test_that("grid_space_filling() errors with non-param inputs", {
  # default method
  expect_snapshot(error = TRUE, grid_space_filling())
  expect_snapshot(error = TRUE, grid_space_filling("not a param"))

  # param method
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), "min_n"))
  expect_snapshot(error = TRUE, grid_space_filling(mtry(), "min_n"))

  # list method
  expect_snapshot(error = TRUE, grid_space_filling(list()))
  expect_snapshot(error = TRUE, grid_space_filling(list(penalty(), "min_n")))
  expect_snapshot(error = TRUE, grid_space_filling(list(mtry(), "min_n")))
})

test_that("grid_space_filling.parameters() checks for NA", {
  p <- parameters(penalty())
  p <- update(p, penalty = NA)
  expect_snapshot(error = TRUE, grid_space_filling(p))
})

test_that("grid_space_filling() errors with params containing unknowns", {
  # parameters method
  expect_snapshot(error = TRUE, grid_space_filling(parameters(mtry())))

  # param method
  expect_snapshot(error = TRUE, grid_space_filling(mtry()))
  expect_snapshot(error = TRUE, grid_space_filling(mtry(), sample_size()))

  # list method
  expect_snapshot(error = TRUE, grid_space_filling(list(mtry())))
  expect_snapshot(
    error = TRUE,
    grid_space_filling(list(mtry_custom_name = mtry()))
  )
  expect_snapshot(error = TRUE, grid_space_filling(list(mtry(), sample_size())))
})

test_that("grid_space_filling() errors with duplicate parameter ids", {
  # param method
  expect_snapshot(error = TRUE, grid_space_filling(penalty(), penalty()))

  # list method
  expect_snapshot(
    error = TRUE,
    grid_space_filling(list(a = penalty(), a = mtry()))
  )
})

test_that("pre-made designs respect the 'original argument", {
  # See issue #409

  prms <- parameters(penalty(), mixture())
  types <- c(
    "audze_eglais",
    "max_min_l1",
    "max_min_l2",
    "uniform",
    "max_entropy",
    "latin_hypercube"
  )

  for (i in types) {
    grd_orig <- grid_space_filling(prms, type = i, original = TRUE)
    grd_trans <- grid_space_filling(prms, type = i, original = FALSE)

    expect_true(all(grd_orig$penalty > 0))
    expect_true(all(grd_trans$penalty <= 0))
  }
})

# anchors ----------------------------------------------------------------------

test_that("anchors are kept, in order, and new rows are added", {
  prm <- parameters(penalty(), mixture())
  anchors <- tibble::tibble(
    penalty = c(1e-4, 1e-2),
    mixture = c(0.1, 0.9)
  )

  set.seed(1)
  grid <- grid_space_filling(prm, size = 6, anchors = anchors)

  expect_equal(nrow(grid), 6L)
  # The anchors are returned first, unchanged and in their original order.
  expect_identical(grid[1:2, ], anchors)
  # New rows are unique and are not anchors.
  expect_equal(anyDuplicated(as.data.frame(grid)), 0L)
  new_rows <- grid[3:6, ]
  anchor_keys <- paste(anchors$penalty, anchors$mixture)
  expect_false(any(paste(new_rows$penalty, new_rows$mixture) %in% anchor_keys))
})

test_that("anchors are added to the `param` and `list` methods", {
  anchors <- tibble::tibble(
    penalty = c(1e-4, 1e-2),
    mixture = c(0.1, 0.9)
  )

  from_param <- grid_space_filling(
    penalty(),
    mixture(),
    size = 5,
    anchors = anchors
  )
  from_list <- grid_space_filling(
    list(penalty(), mixture()),
    size = 5,
    anchors = anchors
  )

  expect_identical(from_param[1:2, ], anchors)
  expect_identical(from_list[1:2, ], anchors)
  expect_identical(from_param, from_list)
})

test_that("new anchors do not depend on anchor order", {
  prm <- parameters(penalty(), mixture())
  anchors <- tibble::tibble(
    penalty = c(1e-4, 1e-2, 1),
    mixture = c(0.05, 0.5, 0.95)
  )
  shuffled <- anchors[c(3, 1, 2), ]

  grid_1 <- grid_space_filling(prm, size = 8, anchors = anchors)
  grid_2 <- grid_space_filling(prm, size = 8, anchors = shuffled)

  expect_identical(grid_1[1:3, ], anchors)
  expect_identical(grid_2[1:3, ], shuffled)

  new_1 <- grid_1[4:8, ]
  new_2 <- grid_2[4:8, ]
  new_1 <- new_1[order(new_1$penalty, new_1$mixture), ]
  new_2 <- new_2[order(new_2$penalty, new_2$mixture), ]
  expect_identical(new_1, new_2)
})

test_that("anchor augmentation is deterministic", {
  prm <- parameters(penalty(), mixture())
  anchors <- tibble::tibble(penalty = c(1e-4, 1e-2), mixture = c(0.1, 0.9))

  grid_1 <- grid_space_filling(prm, size = 6, anchors = anchors)
  grid_2 <- grid_space_filling(prm, size = 6, anchors = anchors)
  expect_identical(grid_1, grid_2)
})

test_that("anchors are preserved exactly", {
  prm <- parameters(penalty(), mixture())
  penalty_values <- 10^seq(-9, 0, length.out = 4)
  mixture_values <- seq(0.01, 0.99, length.out = 4)
  anchors <- tibble::tibble(
    penalty = penalty_values,
    mixture = mixture_values
  )

  grid <- grid_space_filling(prm, size = 8, anchors = anchors)
  expect_identical(grid$penalty[1:4], penalty_values)
  expect_identical(grid$mixture[1:4], mixture_values)
})

test_that("integer anchors return tuning-ready integer values", {
  prm <- parameters(trees(), min_n())
  anchors <- tibble::tibble(trees = c(1L, 2000L), min_n = c(2L, 40L))

  grid <- grid_space_filling(prm, size = 8, anchors = anchors)

  expect_identical(grid[1:2, ], anchors)
  expect_type(grid$trees, "integer")
  expect_type(grid$min_n, "integer")
  expect_true(all(grid$trees >= 1 & grid$trees <= 2000))
  expect_true(all(grid$min_n >= 2 & grid$min_n <= 40))
  expect_equal(anyDuplicated(as.data.frame(grid)), 0L)
})

test_that("the maximin rule maximizes the minimum distance", {
  prm <- parameters(mixture())
  anchors <- tibble::tibble(mixture = c(0, 1))

  grid <- grid_space_filling(prm, size = 5, anchors = anchors)
  # The new points split the remaining space into equal gaps.
  values <- sort(grid$mixture)
  expect_equal(values, c(0, 0.25, 0.5, 0.75, 1), tolerance = 1e-6)
  expect_true(min(diff(values)) > 0.2)

  # Adding one more point makes the minimum gap strictly smaller.
  bigger <- grid_space_filling(prm, size = 6, anchors = anchors)
  expect_lt(min(diff(sort(bigger$mixture))), min(diff(values)))
})

test_that("fixed parameters are held constant", {
  prm <- parameters(penalty(), mixture() |> value_set(0.5))
  anchors <- tibble::tibble(penalty = c(1e-4, 1e-2), mixture = 0.5)

  expect_silent(withr::with_options(
    list(warn = 2),
    grid <- grid_space_filling(prm, size = 5, anchors = anchors)
  ))

  expect_equal(nrow(grid), 5L)
  expect_true(all(grid$mixture == 0.5))
  expect_identical(grid[1:2, ], anchors)
})

test_that("fixed parameter spaces report that nothing can be added", {
  prm <- parameters(mixture() |> value_set(0.5))
  anchors <- tibble::tibble(mixture = 0.5)

  expect_snapshot(
    grid <- grid_space_filling(prm, size = 3, anchors = anchors)
  )
  expect_equal(grid, anchors)
})

test_that("exhausted discrete spaces return only the remaining rows", {
  prm <- parameters(
    penalty() |> value_set(-4:-1),
    mixture() |> value_set(c(0.1, 0.5))
  )
  anchors <- tibble::tibble(penalty = 1e-4, mixture = 0.1)

  expect_snapshot(
    grid <- grid_space_filling(prm, size = 10, anchors = anchors)
  )

  # 4 * 2 feasible combinations in total, one of which is an anchor.
  expect_equal(nrow(grid), 8L)
  expect_equal(anyDuplicated(as.data.frame(grid)), 0L)
  expect_identical(grid[1, ], anchors)
})

test_that("zero-row anchors fall back to a standard design", {
  prm <- parameters(penalty(), mixture())
  empty <- tibble::tibble(penalty = numeric(0), mixture = numeric(0))

  set.seed(1)
  with_empty <- grid_space_filling(prm, size = 4, type = "uniform", anchors = empty)
  standard <- grid_space_filling(prm, size = 4, type = "uniform")

  expect_identical(with_empty, standard)
})

test_that("`size` no larger than the anchors returns only the anchors", {
  prm <- parameters(penalty(), mixture())
  anchors <- tibble::tibble(penalty = c(1e-4, 1e-2), mixture = c(0.1, 0.9))

  expect_snapshot(
    grid <- grid_space_filling(prm, size = 2, anchors = anchors)
  )
  expect_identical(grid, anchors)

  expect_warning(
    sz_one <- grid_space_filling(prm, size = 1, anchors = anchors),
    "smaller than the number of anchors"
  )
  expect_identical(sz_one, anchors)
})

test_that("anchors work with transformations", {
  prm <- parameters(batch_size())
  anchors <- tibble::tibble(batch_size = c(4L, 128L))

  grid <- grid_space_filling(prm, size = 5, anchors = anchors)
  expect_identical(grid[1:2, ], anchors)
  expect_type(grid$batch_size, "integer")

  # With `original = FALSE`, the anchors are on the transformed scale.
  trans_anchors <- tibble::tibble(batch_size = c(2, 7))
  trans_grid <- grid_space_filling(
    prm,
    size = 4,
    original = FALSE,
    anchors = trans_anchors
  )
  expect_identical(trans_grid[1:2, ], trans_anchors)
})

test_that("anchors work with qualitative parameters", {
  prm <- parameters(penalty(), activation(c("relu", "tanh", "elu")))
  anchors <- tibble::tibble(
    penalty = c(1e-4, 1e-2),
    activation = c("relu", "elu")
  )

  grid <- grid_space_filling(prm, size = 5, anchors = anchors)

  expect_identical(grid[1:2, ], anchors)
  expect_equal(nrow(grid), 5L)
  expect_true(all(grid$activation %in% c("relu", "tanh", "elu")))
})

test_that("`anchors` validates its input", {
  prm <- parameters(penalty(), mixture())

  expect_snapshot(
    error = TRUE,
    grid_space_filling(prm, anchors = list(penalty = 0.1, mixture = 0.5))
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(prm, anchors = tibble::tibble(penalty = 0.1))
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      prm,
      anchors = tibble::tibble(penalty = 0.1, mixture = 0.5, other = 1)
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      prm,
      anchors = tibble::tibble(penalty = 0.1, mixture = NA_real_)
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      prm,
      anchors = tibble::tibble(penalty = c(0.1, 0.1), mixture = c(0.5, 0.5))
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      prm,
      anchors = tibble::tibble(penalty = 0.1, mixture = 1.5)
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      parameters(trees()),
      anchors = tibble::tibble(trees = 2.5)
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      parameters(penalty(), activation(c("relu", "tanh"))),
      anchors = tibble::tibble(penalty = 0.1, activation = "foo")
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      parameters(penalty()),
      anchors = tibble::tibble(penalty = "a")
    )
  )
  expect_snapshot(
    error = TRUE,
    grid_space_filling(
      parameters(activation(c("relu", "tanh"))),
      anchors = tibble::tibble(activation = 1)
    )
  )
})
