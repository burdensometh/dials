# activate() rejects unknown target identifiers

    Code
      activate(test_pset(), not_a_param = prune == "cost")
    Condition
      Error in `activate()`:
      ! Activation condition references parameter not in the set: "not_a_param".
      i Parameter identifiers available: "prune", "penalty", and "mixture".

# activate() rejects unknown references

    Code
      activate(test_pset(), penalty = not_a_param == "cost")
    Condition
      Error in `activate()`:
      ! The activation condition for "penalty" references unknown parameter: "not_a_param".
      i Parameter identifiers available: "prune", "penalty", and "mixture".

# activate() rejects self references

    Code
      activate(test_pset(), penalty = penalty > 0.1)
    Condition
      Error in `activate()`:
      ! The activation condition for "penalty" cannot reference "penalty" itself.

# activate() rejects constant conditions

    Code
      activate(test_pset(), penalty = TRUE)
    Condition
      Error in `activate()`:
      ! The activation condition for "penalty" must reference at least one other parameter.

# activate() rejects unnamed and duplicated condition names

    Code
      activate(test_pset(), prune == "cost")
    Condition
      Error in `activate()`:
      ! Every activation condition must be named with a parameter `id`.

---

    Code
      activate(test_pset(), penalty = prune == "cost", penalty = prune == "none")
    Condition
      Error in `activate()`:
      x `activation condition names` must have unique values.
      i Duplicates: "penalty"

# activate() rejects dependency cycles

    Code
      activate(p, mixture = penalty > 0.5)
    Condition
      Error in `activate()`:
      ! Activation conditions contain a dependency cycle; parameters must be processable in dependency order.
      x Parameters involved in the cycle: "penalty" and "mixture".

# activate() rejects conflicting (never satisfiable) conditions

    Code
      activate(test_pset(), penalty = prune == "not_a_value")
    Condition
      Error in `activate()`:
      ! The activation condition for "penalty" can never be satisfied given the possible values of the parameters it references.
      i Check the condition and the referenced parameter values.

# activate() rejects non-logical conditions

    Code
      activate(test_pset(), penalty = prune)
    Condition
      Error in `activate()`:
      ! The activation condition for "penalty" must evaluate to a logical vector, not a character vector.

# strict probing does not drop qualitative levels

    Code
      activate(p, c = a == "cost" & b == "no_such_activation")
    Condition
      Error in `activate()`:
      ! The activation condition for "c" can never be satisfied given the possible values of the parameters it references.
      i Check the condition and the referenced parameter values.

# deactivate() rejects unknown identifiers

    Code
      deactivate(activate(test_pset(), penalty = prune == "cost"), not_a_param)
    Condition
      Error in `handlers[[2L]]()`:
      ! `ids` must reference parameter identifiers in the set: "prune", "penalty", and "mixture".

