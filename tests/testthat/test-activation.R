# ------------------------------------------------------------------------------
# Activation conditions
# ------------------------------------------------------------------------------

test_pset <- function() {
  parameters(list(
    prune = prune_method(c("backward", "none", "cost")),
    penalty = penalty(),
    mixture = mixture()
  ))
}

test_that("activate() stores conditions in an activation list-column", {
  p <- test_pset()
  pa <- activate(p, penalty = prune == "cost")

  expect_s3_class(pa, "parameters")
  expect_true("activation" %in% names(pa))
  expect_true(is_quosure(pa$activation[[which(pa$id == "penalty")]]))
  expect_null(pa$activation[[which(pa$id == "prune")]])
  expect_true(has_activation(pa))
  expect_false(has_activation(p))
})

test_that("parameters() and parameters_constr() include an empty activation column", {
  p <- test_pset()
  expect_true("activation" %in% names(p))
  expect_true(all(map_lgl(p$activation, is.null)))

  pc <- parameters_constr(
    name = "a",
    id = "a",
    source = "list",
    component = "unknown",
    component_id = "unknown",
    object = list(penalty())
  )
  expect_true("activation" %in% names(pc))
  expect_null(pc$activation[[1]])
})

test_that("activate() replaces an existing condition for the same parameter", {
  p <- activate(test_pset(), penalty = prune == "cost")
  p2 <- activate(p, penalty = prune == "none")
  quo <- p2$activation[[which(p2$id == "penalty")]]
  expect_identical(quo_get_expr(quo), expr(prune == "none"))
})

test_that("activate() rejects unknown target identifiers", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), not_a_param = prune == "cost")
  })
})

test_that("activate() rejects unknown references", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), penalty = not_a_param == "cost")
  })
})

test_that("activate() rejects self references", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), penalty = penalty > 0.1)
  })
})

test_that("activate() rejects constant conditions", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), penalty = TRUE)
  })
})

test_that("activate() rejects unnamed and duplicated condition names", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), prune == "cost")
  })
  expect_snapshot(error = TRUE, {
    activate(
      test_pset(),
      penalty = prune == "cost",
      penalty = prune == "none"
    )
  })
})

test_that("activate() rejects dependency cycles", {
  p <- activate(test_pset(), penalty = mixture > 0.5, .strict = FALSE)
  expect_snapshot(error = TRUE, {
    activate(p, mixture = penalty > 0.5)
  })
})

test_that("activate() rejects conflicting (never satisfiable) conditions", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), penalty = prune == "not_a_value")
  })
  # strict = FALSE keeps the condition, the plan is still built
  p <- activate(
    test_pset(),
    penalty = prune == "not_a_value",
    .strict = FALSE
  )
  expect_true(has_activation(p))
})

test_that("activate() rejects non-logical conditions", {
  expect_snapshot(error = TRUE, {
    activate(test_pset(), penalty = prune)
  })
})

test_that("activate() works with unknown parameters not referenced by a condition", {
  p <- parameters(list(
    prune = prune_method(c("backward", "none", "cost")),
    mtry = mtry()
  ))
  expect_no_condition(activate(p, mtry = prune == "cost"))
})

test_that("activate() handles conditions referencing several parameters", {
  p <- parameters(list(
    a = prune_method(c("backward", "none", "cost")),
    b = activation(c("relu", "tanh", "elu")),
    c = penalty()
  ))
  expect_no_condition(
    p <- activate(p, c = a %in% c("cost", "none") & b == "relu")
  )
  g <- grid_regular(p, levels = 3)
  active <- g$a %in% c("cost", "none") & g$b == "relu"
  expect_true(all(is.na(g$c) == !active))
})

test_that("activate() accepts conditions on quantitative parameters", {
  p <- parameters(list(t = trees(c(1, 50)), pen = penalty()))
  p <- activate(p, pen = t >= 25)
  g <- grid_regular(p, levels = 4)
  expect_true(all(is.na(g$pen) == !(g$t >= 25)))
})

test_that("strict probing does not drop qualitative levels", {
  p <- parameters(list(
    a = prune_method(c("backward", "none", "cost")),
    b = activation(c("relu", "tanh", "elu")),
    c = penalty()
  ))
  # satisfiable only through a single level combination
  expect_no_condition(
    activate(p, c = a == "cost" & b == "elu")
  )
  # a genuinely impossible combination is still rejected
  expect_snapshot(error = TRUE, {
    activate(p, c = a == "cost" & b == "no_such_activation")
  })
})

test_that("deactivate() removes conditions", {
  p <- activate(
    test_pset(),
    penalty = prune == "cost",
    mixture = prune != "none"
  )

  p1 <- deactivate(p, penalty)
  expect_null(p1$activation[[which(p1$id == "penalty")]])
  expect_true(!is.null(p1$activation[[which(p1$id == "mixture")]]))

  p2 <- deactivate(p, c(penalty, mixture))
  expect_true(all(map_lgl(p2$activation, is.null)))

  p3 <- deactivate(p)
  expect_true(all(map_lgl(p3$activation, is.null)))
})

test_that("deactivate() rejects unknown identifiers", {
  expect_snapshot(error = TRUE, {
    deactivate(activate(test_pset(), penalty = prune == "cost"), not_a_param)
  })
})

# ------------------------------------------------------------------------------
# Regular grid
# ------------------------------------------------------------------------------

test_that("grid_regular() honors activation while generating candidates", {
  p <- activate(test_pset(), penalty = prune == "cost")
  g <- grid_regular(p, levels = 3)

  # inactive cells carry the explicit, type-consistent NA sentinel
  expect_true(is.double(g$penalty))
  expect_true(all(is.na(g$penalty[g$prune != "cost"])))
  # active rows always carry real values
  expect_true(all(!is.na(g$penalty[g$prune == "cost"])))
  # the inactive parameter never takes part in the Cartesian expansion:
  # "cost" rows expand over penalty x mixture, other rows only over mixture
  n_prune <- length(unique(g$prune))
  expect_equal(nrow(g), 3 * 3 + (n_prune - 1) * 3)
  # inactive placeholders do not create distinct candidates
  expect_equal(nrow(unique(g)), nrow(g))
})

test_that("grid_regular() expands active rows over the full value sequence", {
  p <- activate(test_pset(), penalty = prune == "cost")
  g <- grid_regular(p, levels = 3)
  cost_rows <- g[g$prune == "cost", ]
  expect_equal(length(unique(cost_rows$penalty)), 3)
  expect_equal(length(unique(cost_rows$mixture)), 3)
})

test_that("grid_regular() handles chained dependencies in order", {
  p <- parameters(list(
    a = prune_method(c("backward", "none", "cost")),
    b = penalty(),
    c = mixture()
  ))
  p <- activate(p, b = a == "cost", c = b > 0.1)
  g <- grid_regular(p, levels = 4)

  expect_true(all(is.na(g$b) == (g$a != "cost")))
  expect_true(all(is.na(g$c) == (is.na(g$b) | g$b <= 0.1)))
})

test_that("grid_regular() keeps filter for post-generation row filtering", {
  p <- activate(test_pset(), penalty = prune == "cost")
  g <- grid_regular(p, levels = 3, filter = mixture >= 0.5)
  expect_true(all(g$mixture >= 0.5))
  expect_true(all(is.na(g$penalty[g$prune != "cost"])))
})

test_that("grid_regular() without activation is unchanged", {
  p <- test_pset()
  expected <- expand.grid(
    prune = value_seq(p$object[[1]], 3),
    penalty = value_seq(p$object[[2]], 3),
    mixture = value_seq(p$object[[3]], 3),
    stringsAsFactors = FALSE
  )
  attr(expected, "out.attrs") <- NULL
  expect_equal(grid_regular(p, levels = 3), tibble::as_tibble(expected))
})

# ------------------------------------------------------------------------------
# Random grid
# ------------------------------------------------------------------------------

test_that("grid_random() honors activation while sampling candidates", {
  p <- activate(test_pset(), penalty = prune == "cost")
  set.seed(1)
  g <- grid_random(p, size = 50)

  expect_equal(nrow(g), 50)
  expect_true(all(is.na(g$penalty[g$prune != "cost"])))
  expect_true(all(!is.na(g$penalty[g$prune == "cost"])))
  # sentinels do not make inactive rows distinct candidates
  expect_equal(nrow(unique(g)), nrow(g))
})

test_that("grid_random() handles chained dependencies", {
  p <- parameters(list(
    a = prune_method(c("backward", "none", "cost")),
    b = penalty(),
    c = mixture()
  ))
  p <- activate(p, b = a == "cost", c = b > 0.1)
  set.seed(2)
  g <- grid_random(p, size = 40)

  expect_true(all(is.na(g$b) == (g$a != "cost")))
  expect_true(all(is.na(g$c) == (is.na(g$b) | g$b <= 0.1)))
})

test_that("grid_random() supports filter together with activation", {
  p <- activate(test_pset(), penalty = prune == "cost")
  g <- with_seed(
    11L,
    grid_random(p, size = 50, filter = mixture >= 0.5)
  )
  expect_true(all(g$mixture >= 0.5))
  expect_true(all(is.na(g$penalty[g$prune != "cost"])))
})

# ------------------------------------------------------------------------------
# Hybrid grid
# ------------------------------------------------------------------------------

test_that("grid_hybrid() applies the same activation plan", {
  p <- activate(test_pset(), penalty = prune == "cost")
  h <- grid_hybrid(
    p,
    reg_parameters = prune,
    sfd_size = 8,
    reg_levels = 3
  )

  expect_true(all(is.na(h$penalty[h$prune != "cost"])))
  expect_true(all(!is.na(h$penalty[h$prune == "cost"])))
  expect_equal(nrow(unique(h)), nrow(h))
})

test_that("grid_hybrid() honors activation when all parameters are space-filling", {
  p <- activate(test_pset(), penalty = prune == "cost")
  h <- grid_hybrid(
    p,
    reg_parameters = character(0),
    sfd_size = 12
  )
  expect_true(all(is.na(h$penalty[h$prune != "cost"])))
  expect_true(all(!is.na(h$penalty[h$prune == "cost"])))
})

test_that("grid_hybrid() honors activation when all parameters are regular", {
  p <- activate(test_pset(), penalty = prune == "cost")
  h <- grid_hybrid(p, reg_parameters = everything(), reg_levels = 2)
  expect_true(all(is.na(h$penalty[h$prune != "cost"])))
})

test_that("grid_hybrid() handles conditions crossing the sfd/regular split", {
  # penalty is a regular parameter but its condition references prune, which is
  # not part of the regular subset.
  p <- activate(test_pset(), penalty = prune == "cost")
  expect_no_condition(
    h <- grid_hybrid(
      p,
      reg_parameters = penalty,
      sfd_size = 6,
      reg_levels = 3
    )
  )
  expect_true(all(is.na(h$penalty[h$prune != "cost"])))
})

# ------------------------------------------------------------------------------
# Shared plan and sentinel semantics
# ------------------------------------------------------------------------------

test_that("all three grid functions agree on inactive cells", {
  p <- activate(test_pset(), penalty = prune == "cost")

  g <- grid_regular(p, levels = 3)
  set.seed(7)
  r <- grid_random(p, size = 60)
  h <- grid_hybrid(p, reg_parameters = prune, sfd_size = 6, reg_levels = 2)

  inactive <- function(d) is.na(d$penalty) == (d$prune != "cost")
  expect_true(all(inactive(g)))
  expect_true(all(inactive(r)))
  expect_true(all(inactive(h)))
})

test_that("inactive sentinels match the parameter type", {
  p <- parameters(list(
    prune = prune_method(c("backward", "none", "cost")),
    balance = balance_probabilities(),
    trees = trees()
  ))
  p <- activate(
    p,
    balance = prune == "cost",
    trees = prune != "none"
  )
  g <- grid_regular(p, levels = 2)

  expect_type(g$balance, "logical")
  expect_type(g$trees, "integer")
  expect_type(g$prune, "character")
  expect_true(all(is.na(g$balance[g$prune != "cost"])))
  expect_true(all(is.na(g$trees[g$prune == "none"])))
})

test_that("activation survives parameters row subsetting", {
  p <- activate(test_pset(), penalty = prune == "cost")
  sub <- p[1:2, ]
  expect_s3_class(sub, "parameters")
  expect_true("activation" %in% names(sub))
})

test_that("deactivated sets fall back to full Cartesian expansion", {
  p <- activate(test_pset(), penalty = prune == "cost")
  pd <- deactivate(p)
  expect_equal(
    nrow(grid_regular(pd, levels = 3)),
    3 * 3 * 3
  )
})

test_that("activation is generation, not post-generation row filtering", {
  p <- activate(test_pset(), penalty = prune == "cost")
  g <- grid_regular(p, levels = 3)

  # A post-generation filter on the full Cartesian product could only remove
  # rows, never mark a cell inactive: every inactive row survives with an NA
  # cell, which cannot be produced by filtering complete combinations.
  inactive_rows <- g[g$prune != "cost", ]
  expect_gt(nrow(inactive_rows), 0L)
  expect_true(all(is.na(inactive_rows$penalty)))
  # and the inactive parameter contributes no expansion factor: each prune x
  # mixture combination appears once, with cost rows additionally expanded
  # over the 3 penalty values
  base <- nrow(expand.grid(
    prune = value_seq(p$object[[1]], 3),
    mixture = value_seq(p$object[[3]], 3),
    stringsAsFactors = FALSE
  ))
  expect_equal(nrow(g), base + 3 * 3 - 3)
  expect_equal(nrow(g), 15L)
})

test_that("grid_random() draws the requested number of candidates", {
  p <- activate(test_pset(), penalty = prune == "cost")
  set.seed(4)
  g <- grid_random(p, size = 25)
  expect_equal(nrow(g), 25L)
})
