# Augmenting a space-filling design with fixed anchor combinations.
#
# Supplied anchors are kept exactly as-is and are never moved. New combinations
# are added from a deterministic candidate set using a greedy maximin rule:
# each new point maximizes the minimum squared Euclidean distance to the anchors
# and to the points already added, measured on the [0, 1] scale of the
# non-fixed parameters. Ties are broken by the canonical candidate order (which
# does not depend on the anchor order), so reordering `anchors` does not change
# which new points are selected.

max_anchor_enumeration <- 20000L
max_anchor_candidates <- 20000L
max_integer_discrete <- 1e5

make_anchored_sfd <- function(
  params,
  size,
  original,
  anchors,
  call = caller_env()
) {
  info <- map2(params, names(params), function(param, nm) {
    param_augment_info(param, original = original, nm = nm, call = call)
  })
  names(info) <- names(params)

  anchors_tbl <- check_anchor_grid(
    anchors,
    params,
    info,
    original = original,
    call = call
  )
  # A data frame with no rows means "no anchors"; the caller falls back to a
  # standard design.
  if (nrow(anchors_tbl) == 0L) {
    return(NULL)
  }

  active <- map_lgl(info, function(.x) .x$active)

  anchors_tbl <- cast_anchor_grid(anchors_tbl, info, call = call)
  n_anchor <- nrow(anchors_tbl)

  if (size < n_anchor) {
    cli::cli_warn(
      c(
        "!" = "{.arg size} ({size}) is smaller than the number of anchors
               ({n_anchor}).",
        "i" = "All anchors are returned and no combinations are added."
      ),
      call = call
    )
    return(new_param_grid(anchors_tbl))
  }

  n_new <- size - n_anchor
  if (n_new == 0L) {
    return(new_param_grid(anchors_tbl))
  }

  if (!any(active)) {
    cli::cli_warn(
      c(
        "!" = "No combinations can be added.",
        "i" = "Every parameter has a single feasible value, so the anchors
               already cover the whole space."
      ),
      call = call
    )
    return(new_param_grid(anchors_tbl))
  }

  anchor_coords <- coords_from_tbl(info, active, anchors_tbl)

  cand <- make_anchor_candidates(
    info,
    active,
    anchors_tbl,
    anchor_coords,
    n_new = n_new,
    call = call
  )

  n_available <- nrow(cand$tbl)
  if (n_available < n_new) {
    cli::cli_warn(
      c(
        "!" = "Only {n_available} new combination{?s} can be added,
               not {n_new}.",
        "i" = cand$note
      ),
      call = call
    )
    n_new <- n_available
  }
  if (n_new == 0L) {
    return(new_param_grid(anchors_tbl))
  }

  chosen <- greedy_maximin(cand$coords, anchor_coords, n_new)
  if (length(chosen) == 0L) {
    return(new_param_grid(anchors_tbl))
  }

  out <- vctrs::vec_rbind(anchors_tbl, cand$tbl[chosen, , drop = FALSE])
  new_param_grid(out)
}

# Anchor validation ------------------------------------------------------------

check_anchor_grid <- function(anchors, params, info, original, call) {
  check_data_frame(anchors, arg = "anchors", call = call)
  ids <- names(params)

  if (is.null(names(anchors))) {
    cli::cli_abort(
      c(
        "{.arg anchors} must have named columns.",
        "i" = "The columns must be named with the parameter identifiers
               {.val {ids}}."
      ),
      call = call
    )
  }

  missing <- setdiff(ids, names(anchors))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg anchors} does not have a column for every parameter.",
        "x" = "Missing column{?s}: {.var {missing}}.",
        "i" = "Parameter identifiers: {.val {ids}}."
      ),
      call = call
    )
  }
  extra <- setdiff(names(anchors), ids)
  if (length(extra) > 0L) {
    cli::cli_abort(
      c(
        "{.arg anchors} has columns that are not parameters.",
        "x" = "Unexpected column{?s}: {.var {extra}}.",
        "i" = "Parameter identifiers: {.val {ids}}."
      ),
      call = call
    )
  }

  anchors <- anchors[, ids, drop = FALSE]

  na_row <- Reduce(
    `|`,
    lapply(anchors, is.na),
    init = rep(FALSE, nrow(anchors))
  )
  if (any(na_row)) {
    cli::cli_abort(
      c(
        "{.arg anchors} must not contain missing values.",
        "x" = "Missing values in row{?s} {fmt_rows(which(na_row))}."
      ),
      call = call
    )
  }

  for (nm in ids) {
    check_anchor_type(anchors[[nm]], info[[nm]], original, call = call)
  }

  for (nm in ids) {
    bad <- !col_validity(anchors[[nm]], info[[nm]])
    if (any(bad)) {
      abort_bad_anchor(anchors[[nm]], bad, nm, original, call = call)
    }
  }

  dup <- vctrs::vec_duplicate_detect(anchors)
  if (any(dup)) {
    cli::cli_abort(
      c(
        "{.arg anchors} must not contain duplicate combinations.",
        "x" = "Duplicate row{?s}: {fmt_rows(which(dup))}."
      ),
      call = call
    )
  }

  anchors
}

check_anchor_type <- function(vals, info, original, call) {
  nm <- info$nm
  ok <- if (info$kind == "qual") {
    is.character(vals) || is.logical(vals) || is.factor(vals)
  } else {
    is.numeric(vals) && is.null(dim(vals)) && !is.object(vals)
  }
  if (ok) {
    return(invisible(NULL))
  }

  what <- if (info$kind == "qual") {
    "a character or logical vector"
  } else {
    "a numeric vector"
  }
  cli::cli_abort(
    c(
      "The {.arg anchors} column {.val {nm}} must be {what}.",
      "x" = "The column has class {.cls {class(vals)}}."
    ),
    call = call
  )
}

abort_bad_anchor <- function(vals, bad, nm, original, call) {
  offenders <- unique(vals[bad])
  if (length(offenders) > 5L) {
    offenders <- c(offenders[1:5], "...")
  }
  cli::cli_abort(
    c(
      "{.arg anchors} has values that are not valid for {.val {nm}}.",
      "x" = "{length(offenders)} invalid value{?s}: {.val {offenders}}.",
      "i" = if (original) {
        "Anchors must be in the original units."
      } else {
        "Anchors must be on the transformed scale."
      }
    ),
    call = call
  )
}

fmt_rows <- function(rows) {
  paste(rows, collapse = ", ")
}

# Is each value a feasible tuning value for the parameter behind `info`?
col_validity <- function(vals, info) {
  if (info$discrete) {
    return(!is.na(value_index(vals, info$values)))
  }

  param <- info$param
  internal <- if (info$original) value_transform(param, vals) else vals
  valid <- value_validate(param, internal)
  valid[is.na(valid)] <- FALSE

  # Integer parameters must be whole on their natural scale. When the values
  # are on the transformed scale and a transformation is present, the scale is
  # genuinely continuous (as with the other `grid_*()` functions).
  if (param$type == "integer" && (info$original || is.null(param$trans))) {
    valid <- valid & is_whole(vals)
  }

  valid
}

is_whole <- function(x) {
  is.finite(x) & abs(x - round(x)) < sqrt(.Machine$double.eps)
}

# Re-check the decoded anchor values and give the columns the same types as the
# design that `grid_space_filling()` would normally return.
cast_anchor_grid <- function(anchors, info, call) {
  proto <- decode_units(
    info,
    stats::setNames(rep(0.5, length(info)), names(info))
  )
  out <- map2(anchors, names(anchors), function(col, nm) {
    if (is.factor(col) && is.character(proto[[nm]])) {
      col <- as.character(col)
    }
    vctrs::vec_cast(col, proto[[nm]])
  })
  names(out) <- names(anchors)
  tibble::as_tibble(out)
}

# Parameter information --------------------------------------------------------

param_augment_info <- function(param, original, nm, call = caller_env()) {
  if (inherits(param, "qual_param")) {
    values <- unique(param$values)
    if (length(values) == 0L) {
      cli::cli_abort("The parameter {.val {nm}} has no values.", call = call)
    }
    n <- length(values)
    coords <- if (n == 1L) 0 else (seq_len(n) - 1) / (n - 1)
    return(discrete_info(
      values,
      coords,
      active = n > 1L,
      nm = nm,
      kind = "qual"
    ))
  }

  rng <- unlist(range_get(param, original = original))

  if (!is.null(param$values)) {
    values <- if (original) value_inverse(param, param$values) else param$values
    values <- values[!is.na(values)]
    coords <- quant_coords(param, values, original)
    ord <- order(coords)
    return(discrete_info(
      values[ord],
      coords[ord],
      active = length(values) > 1L,
      nm = nm,
      kind = "quant"
    ))
  }

  if (all(is.finite(rng)) && rng[[1]] == rng[[2]]) {
    fixed <- encode_unit(param, 0.5, "backward", original = original)
    return(discrete_info(fixed, 0, active = FALSE, nm = nm, kind = "quant"))
  }

  if (any(!is.finite(rng))) {
    cli::cli_abort(
      c(
        "The parameter {.val {nm}} has an unbounded range and cannot be
         augmented with {.arg anchors}.",
        "i" = "Set a finite range for this parameter before making the grid."
      ),
      call = call
    )
  }

  integer_scale <- original || is.null(param$trans)
  if (
    param$type == "integer" &&
      integer_scale &&
      integer_span(param, rng) <= max_integer_discrete
  ) {
    values <- integer_values(param, rng)
    coords <- quant_coords(param, values, original)
    return(discrete_info(
      values,
      coords,
      active = length(values) > 1L,
      nm = nm,
      kind = "quant"
    ))
  }

  list(
    active = TRUE,
    discrete = FALSE,
    values = NULL,
    coords = NULL,
    param = param,
    original = original,
    nm = nm,
    kind = "quant"
  )
}

discrete_info <- function(values, coords, active, nm, kind) {
  if (length(values) == 0L) {
    cli::cli_abort(
      "The parameter {.val {nm}} has no feasible values to anchor.",
      call = caller_env()
    )
  }
  list(
    active = active,
    discrete = TRUE,
    values = values,
    coords = coords,
    param = NULL,
    original = NULL,
    nm = nm,
    kind = kind
  )
}

integer_span <- function(param, rng) {
  bounds <- integer_bounds(param, rng)
  if (bounds[[1]] > bounds[[2]]) {
    return(0)
  }
  bounds[[2]] - bounds[[1]] + 1
}

integer_bounds <- function(param, rng) {
  lo <- if (param$inclusive[["lower"]]) {
    ceiling(rng[[1]])
  } else {
    floor(rng[[1]]) + 1
  }
  hi <- if (param$inclusive[["upper"]]) {
    floor(rng[[2]])
  } else {
    ceiling(rng[[2]]) - 1
  }
  c(lo, hi)
}

integer_values <- function(param, rng) {
  bounds <- integer_bounds(param, rng)
  if (bounds[[1]] > bounds[[2]]) {
    return(integer(0))
  }
  seq.int(bounds[[1]], bounds[[2]])
}

quant_coords <- function(param, values, original) {
  internal <- if (original) value_transform(param, values) else values
  encode_unit(param, internal, "forward")
}

# Candidate generation ---------------------------------------------------------

make_anchor_candidates <- function(
  info,
  active,
  anchors,
  anchor_coords,
  n_new,
  call
) {
  active_info <- info[active]
  all_discrete <- all(map_lgl(active_info, function(.x) .x$discrete))

  if (all_discrete) {
    total <- prod(map_dbl(active_info, function(.x) length(.x$values)))
    if (total <= max_anchor_enumeration) {
      cand <- enumerate_candidates(info)
      return(filter_candidates(
        cand,
        info,
        active,
        anchors,
        anchor_coords,
        note = "The discrete parameter space is exhausted.",
        call = call
      ))
    }
  }

  n_cand <- min(
    max(1000, 40 * n_new, 10 * (nrow(anchors) + n_new)),
    max_anchor_candidates
  )
  cand <- halton_candidates(info, active, n_cand)
  note <- paste0(
    "The candidate pool (",
    n_cand,
    " combinations) may not cover the entire parameter space."
  )
  filter_candidates(
    cand,
    info,
    active,
    anchors,
    anchor_coords,
    note = note,
    call = call
  )
}

# Every discrete combination. `expand.grid()` recycles the single value that
# inactive (fixed) parameters contribute.
enumerate_candidates <- function(info) {
  values <- map(info, function(.x) .x$values)
  grid <- expand.grid(
    values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  names(grid) <- names(info)
  tibble::as_tibble(grid)
}

halton_candidates <- function(info, active, n_cand) {
  active_info <- info[active]
  units <- halton_sequence(n_cand, length(active_info))

  out <- map(info, function(.x) if (.x$active) NA else .x$values[1])
  names(out) <- names(info)

  for (j in seq_along(active_info)) {
    nm <- names(active_info)[[j]]
    out[[nm]] <- decode_dim(units[, j], active_info[[j]])
  }
  tibble::as_tibble(out)
}

decode_dim <- function(unit, info) {
  if (info$discrete) {
    return(info$values[snap_index(unit, info$coords)])
  }
  decode_one(unit, info$param, info$original)
}

decode_one <- function(unit, param, original) {
  value <- encode_unit(param, unit, "backward", original = original)
  # Return genuinely tunable integers on the natural scale.
  if (param$type == "integer" && (original || is.null(param$trans))) {
    rng <- unlist(range_get(param, original = original))
    bounds <- integer_bounds(param, rng)
    if (bounds[[1]] <= bounds[[2]]) {
      value <- pmin(pmax(round(value), bounds[[1]]), bounds[[2]])
    }
    value <- as.integer(value)
  }
  value
}

# Remove candidates that decode to missing values, that duplicate an anchor, or
# that duplicate each other. The returned coordinates are recomputed from the
# decoded values, which re-checks distances after integer rounding.
#
# Candidates are feasible by construction (discrete values are copied from the
# parameter's feasible set and continuous values come from `encode_unit()` on
# [0, 1]), so they are not re-validated here; anchors were validated up front.
filter_candidates <- function(
  cand,
  info,
  active,
  anchors,
  anchor_coords,
  note,
  call
) {
  empty <- list(
    tbl = cand[0, , drop = FALSE],
    coords = matrix(numeric(0), 0, length(active)),
    note = note
  )
  if (nrow(cand) == 0L) {
    return(empty)
  }

  ok <- rep(TRUE, nrow(cand))
  for (nm in names(info)) {
    ok <- ok & !is.na(cand[[nm]])
  }
  cand <- cand[ok, , drop = FALSE]
  if (nrow(cand) == 0L) {
    return(empty)
  }

  cand <- cand[vctrs::vec_unique_loc(cand), , drop = FALSE]
  coords <- coords_from_tbl(info, active, cand)

  # Anchors can be a few ulps away from a feasible value, so compare on the
  # canonical unit coordinates instead of the user-supplied values.
  dup <- !is.na(
    vctrs::vec_match(as.data.frame(coords), as.data.frame(anchor_coords))
  )
  cand <- cand[!dup, , drop = FALSE]
  coords <- coords[!dup, , drop = FALSE]
  if (nrow(cand) == 0L) {
    return(empty)
  }

  list(tbl = cand, coords = coords, note = note)
}

# Coordinates in the shared [0, 1] space, recomputed from actual values.
coords_from_tbl <- function(info, active, tbl) {
  active_info <- info[active]
  n <- nrow(tbl)
  if (length(active_info) == 0L) {
    return(matrix(numeric(0), n, 0))
  }

  coords <- matrix(0, n, length(active_info))
  for (j in seq_along(active_info)) {
    nm <- names(active_info)[[j]]
    .x <- active_info[[j]]
    coords[, j] <- if (.x$discrete) {
      .x$coords[value_index(tbl[[nm]], .x$values)]
    } else {
      quant_coords(.x$param, tbl[[nm]], .x$original)
    }
  }
  colnames(coords) <- names(active_info)
  coords
}

# Index of each value within the parameter's feasible values, or `NA`.
#
# Generated candidates are copied from `values`, so an exact `match()` finds
# them immediately. User-supplied anchors can be off by a few ulps (for example,
# after a log transform and back), so unmatched values fall back to a tolerant
# nearest-value lookup.
value_index <- function(x, values) {
  idx <- match(x, values)
  if (all(!is.na(idx)) || !is.numeric(values)) {
    return(idx)
  }

  tol <- sqrt(.Machine$double.eps)
  missing <- which(is.na(idx))
  for (k in missing) {
    d <- abs(x[[k]] - values)
    d[is.na(d)] <- Inf
    if (length(d) > 0L && min(d) <= tol * max(1, abs(x[[k]]))) {
      idx[[k]] <- which.min(d)
    }
  }
  idx
}

# Decoding helpers -------------------------------------------------------------

decode_units <- function(info, units) {
  out <- map2(info, names(info), function(.x, nm) {
    if (.x$discrete) {
      return(.x$values[snap_index(units[[nm]], .x$coords)])
    }
    decode_one(units[[nm]], .x$param, .x$original)
  })
  names(out) <- names(info)
  out
}

# Nearest feasible value in unit coordinates.
snap_index <- function(unit, coords) {
  n <- length(coords)
  if (n == 1L) {
    return(rep(1L, length(unit)))
  }
  idx <- findInterval(unit, coords)
  idx <- pmax(1L, pmin(n - 1L, idx))
  right <- idx + 1L
  pick_right <- abs(coords[right] - unit) < abs(unit - coords[idx])
  idx[pick_right] <- right[pick_right]
  idx
}

# Halton sequence --------------------------------------------------------------

halton_sequence <- function(n, dim) {
  bases <- first_primes(dim)
  out <- matrix(0, n, dim)
  for (j in seq_len(dim)) {
    out[, j] <- radical_inverse(seq_len(n), bases[[j]])
  }
  out
}

radical_inverse <- function(index, base) {
  out <- numeric(length(index))
  factor <- 1 / base
  while (any(index > 0)) {
    out <- out + factor * (index %% base)
    index <- index %/% base
    factor <- factor / base
  }
  out
}

first_primes <- function(n) {
  primes <- integer(0)
  candidate <- 2L
  while (length(primes) < n) {
    if (all(candidate %% primes != 0)) {
      primes <- c(primes, candidate)
    }
    candidate <- candidate + 1L
  }
  primes
}

# Greedy maximin ---------------------------------------------------------------

# Indices (into the rows of `cand`) of the selected new points.
greedy_maximin <- function(cand, anchors, n_new) {
  n <- nrow(cand)
  if (n == 0L || n_new == 0L || ncol(cand) == 0L) {
    return(integer(0))
  }

  min_d2 <- if (nrow(anchors) == 0L) {
    rep(Inf, n)
  } else {
    min_sq_dist(cand, anchors)
  }

  selected <- integer(0)
  for (k in seq_len(min(n_new, n))) {
    best <- which.max(min_d2)
    if (length(best) == 0L || !isTRUE(min_d2[best] > 0)) {
      break
    }
    selected <- c(selected, best)
    min_d2 <- pmin(min_d2, min_sq_dist(cand, cand[best, , drop = FALSE]))
    min_d2[best] <- -Inf
  }
  selected
}

min_sq_dist <- function(x, ref) {
  n <- nrow(x)
  out <- rep(Inf, n)
  for (j in seq_len(nrow(ref))) {
    out <- pmin(out, sq_dist(x, ref[j, ]))
  }
  out
}

# Squared Euclidean distance from every row of `x` to the single point `point`.
# Column-wise arithmetic avoids the overhead of `rowSums()` on a matrix.
sq_dist <- function(x, point) {
  d2 <- (x[, 1L] - point[[1L]])^2
  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) {
      d2 <- d2 + (x[, j] - point[[j]])^2
    }
  }
  d2
}
