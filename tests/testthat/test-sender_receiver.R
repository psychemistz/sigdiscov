# Tests for sender/receiver definition functions

# ==== Visium: Expression-based splitting ====

test_that("split_by_expression returns correct structure", {
    expr <- 1:100
    result <- split_by_expression(expr, percentile = 75)

    expect_type(result, "list")
    expect_true(all(c("sender_idx", "receiver_idx", "threshold",
                      "percentile", "n_senders", "n_receivers") %in% names(result)))
})

test_that("split_by_expression respects percentile", {
    set.seed(42)
    expr <- rnorm(100)
    result <- split_by_expression(expr, percentile = 75)

    # Top 25% should be senders
    expect_true(result$n_senders <= 30)  # About 25% of 100
    expect_true(result$n_receivers >= 70)
})

test_that("split_by_expression uses non-zero values by default", {
    expr <- c(rep(0, 50), 1:50)  # Half zeros

    result_nonzero <- split_by_expression(expr, percentile = 75, use_nonzero = TRUE)
    result_all <- split_by_expression(expr, percentile = 75, use_nonzero = FALSE)

    # With use_nonzero, threshold is based on non-zero values (1-50)
    # 75th percentile of 1:50 is ~37.75
    expect_gt(result_nonzero$threshold, result_all$threshold)
})

test_that("split_by_expression validates input", {
    expect_error(split_by_expression("not numeric"))
    expect_error(split_by_expression(1:10, percentile = 150))
    expect_error(split_by_expression(rep(0, 10)))  # All zeros with use_nonzero = TRUE
})

# ==== Single-cell: Cell type + expression filtering ====

test_that("split_by_celltype returns correct indices", {
    cell_types <- c(rep("A", 30), rep("B", 30), rep("C", 40))

    result <- split_by_celltype(cell_types, "A", "B", min_cells = 10)

    expect_equal(length(result$sender_type_idx), 30)
    expect_equal(length(result$receiver_idx), 30)
    expect_equal(result$sender_type, "A")
    expect_equal(result$receiver_type, "B")
})

test_that("split_by_celltype enforces min_cells", {
    cell_types <- c(rep("A", 5), rep("B", 50))  # Only 5 A cells

    expect_error(split_by_celltype(cell_types, "A", "B", min_cells = 10))
})

test_that("filter_high_expression keeps top 75% by default", {
    set.seed(42)
    expr <- rnorm(100)

    result <- filter_high_expression(expr, quantile_prob = 0.25)

    # Should keep about 75% (those above 25th percentile)
    expect_true(result$n_high >= 70 && result$n_high <= 80)
    expect_equal(result$threshold, as.numeric(quantile(expr, 0.25)))
})

test_that("define_sender_receiver_sc combines cell type and expression", {
    cell_types <- c(rep("A", 100), rep("B", 100))
    factor_expr <- c(rnorm(100, mean = 0), rnorm(100, mean = 5))  # B has higher expression

    result <- define_sender_receiver_sc(
        cell_types, "A", "B",
        factor_expr = factor_expr,
        quantile_prob = 0.25,
        min_cells = 10
    )

    # Senders are A cells with high expression
    expect_true(result$n_senders < result$n_sender_type)  # Filtered
    expect_equal(result$n_receivers, 100)  # All B cells
})

test_that("define_sender_receiver_sc works without expression filter", {
    cell_types <- c(rep("A", 50), rep("B", 50))

    result <- define_sender_receiver_sc(
        cell_types, "A", "B",
        factor_expr = NULL,  # No expression filter
        min_cells = 10
    )

    # All A cells become senders
    expect_equal(result$n_senders, 50)
    expect_true(is.na(result$threshold))
})

test_that("get_celltype_pairs generates all valid pairs", {
    cell_types <- c(rep("A", 20), rep("B", 20), rep("C", 5))

    pairs <- get_celltype_pairs(cell_types, min_cells = 10, exclude_self = TRUE)

    # Only A and B have >= 10 cells
    expect_equal(nrow(pairs), 2)  # A->B and B->A
    expect_true(all(pairs$sender %in% c("A", "B")))
    expect_true(!"C" %in% pairs$sender)
})

test_that("get_celltype_pairs can include self-pairs", {
    cell_types <- c(rep("A", 20), rep("B", 20))

    pairs_no_self <- get_celltype_pairs(cell_types, min_cells = 10, exclude_self = TRUE)
    pairs_with_self <- get_celltype_pairs(cell_types, min_cells = 10, exclude_self = FALSE)

    expect_equal(nrow(pairs_no_self), 2)  # A->B, B->A
    expect_equal(nrow(pairs_with_self), 4)  # A->A, A->B, B->A, B->B
})
