#' Parse10M Analysis Workflow for sigdiscov
#'
#' This module provides a complete workflow for analyzing Parse10M dataset
#' with sender-receiver spatial interactions. It includes spatial layout
#' generation, differential expression filtering, and I_ND computation.
#'
#' Replaces the Python parse10m_unified.py with native R implementation.
#'
#' @name parse10m_analysis
#' @keywords internal
NULL

# =============================================================================
# CONFIGURATION
# =============================================================================

#' Create Parse10M analysis configuration
#'
#' @param h5ad_path Path to Parse10M h5ad file
#' @param output_dir Output directory
#' @param source_cell_type Source/sender cell type
#' @param receiver_cell_type Receiver cell type
#' @param inert_cell_type Inert (non-responsive) cell type
#' @param cytokine_name Cytokine name for treatment
#' @param ligand_gene Ligand gene name
#' @param tissue_radius Tissue radius (um)
#' @param source_region_radius Source region radius (um)
#' @param sigma_decay Expression decay sigma
#' @param p_cytokine_at_source Probability of cytokine-treated cell at source
#' @param p_cytokine_background Probability of cytokine-treated cell in background
#' @param n_source_cells Number of source cells
#' @param n_control_sources Number of control source cells
#' @param n_receiver_cells Number of receiver cells
#' @param n_inert_cells Number of inert cells
#' @param log2fc_threshold Log2 fold change threshold for DE
#' @param pval_threshold P-value threshold for DE
#' @param radii Vector of radii for I_ND analysis
#' @param use_annular Use annular weights
#' @param annular_width Annular width
#' @param response_genes Known response genes
#' @param random_seed Random seed
#' @return Configuration list
#' @export
parse10m_config <- function(
    h5ad_path = "",
    output_dir = "./parse10m_analysis",
    source_cell_type = "CD4 Memory",
    receiver_cell_type = "CD14 Mono",
    inert_cell_type = "B Naive",
    cytokine_name = "IFN-gamma",
    ligand_gene = "IFNG",
    tissue_radius = 5000.0,
    source_region_radius = 500.0,
    sigma_decay = 1500.0,
    p_cytokine_at_source = 0.90,
    p_cytokine_background = 0.10,
    n_source_cells = 500,
    n_control_sources = 500,
    n_receiver_cells = 15000,
    n_inert_cells = 5000,
    log2fc_threshold = 0.5,
    pval_threshold = 0.05,
    radii = seq(50, 5000, by = 50),
    use_annular = TRUE,
    annular_width = 50.0,
    response_genes = c("MX1", "IFITM2", "IFITM3", "IFI6", "IFI27", "IFIT1",
                       "IFIT2", "IFIT3", "IRF1", "GBP1", "GBP2", "STAT1",
                       "B2M", "HLA-A", "HLA-B", "HLA-C"),
    random_seed = 42) {

  list(
    h5ad_path = h5ad_path,
    output_dir = output_dir,
    source_cell_type = source_cell_type,
    receiver_cell_type = receiver_cell_type,
    inert_cell_type = inert_cell_type,
    cytokine_name = cytokine_name,
    ligand_gene = ligand_gene,
    tissue_radius = tissue_radius,
    source_region_radius = source_region_radius,
    sigma_decay = sigma_decay,
    p_cytokine_at_source = p_cytokine_at_source,
    p_cytokine_background = p_cytokine_background,
    n_source_cells = n_source_cells,
    n_control_sources = n_control_sources,
    n_receiver_cells = n_receiver_cells,
    n_inert_cells = n_inert_cells,
    log2fc_threshold = log2fc_threshold,
    pval_threshold = pval_threshold,
    radii = radii,
    use_annular = use_annular,
    annular_width = annular_width,
    response_genes = response_genes,
    random_seed = random_seed
  )
}

# =============================================================================
# SPATIAL LAYOUT GENERATION
# =============================================================================

#' Generate positions in source region (central circular area)
#'
#' @param n Number of positions
#' @param source_region_radius Radius of source region
#' @return Matrix of positions (n x 2)
#' @export
generate_source_positions <- function(n, source_region_radius) {
  positions <- matrix(0, nrow = n, ncol = 2)
  i <- 1

  while (i <= n) {
    x <- runif(1, -source_region_radius, source_region_radius)
    y <- runif(1, -source_region_radius, source_region_radius)

    if (x^2 + y^2 <= source_region_radius^2) {
      positions[i, ] <- c(x, y)
      i <- i + 1
    }
  }

  colnames(positions) <- c("x", "y")
  positions
}

#' Generate positions for control sources
#'
#' @param n Number of positions
#' @param tissue_radius Full tissue radius
#' @param source_region_radius Source region radius (to exclude)
#' @param mode Placement mode: "random", "periphery", or "center"
#' @return Matrix of positions (n x 2)
#' @export
generate_control_positions <- function(n, tissue_radius, source_region_radius,
                                        mode = "random") {
  positions <- matrix(0, nrow = n, ncol = 2)

  if (mode == "random") {
    # Random positions anywhere outside source region
    i <- 1
    while (i <= n) {
      angle <- runif(1) * 2 * pi
      r <- sqrt(runif(1)) * tissue_radius

      if (r > source_region_radius) {
        positions[i, ] <- c(r * cos(angle), r * sin(angle))
        i <- i + 1
      }
    }

  } else if (mode == "periphery") {
    # Positions at the edge of tissue
    min_r <- tissue_radius * 0.8
    i <- 1
    while (i <= n) {
      angle <- runif(1) * 2 * pi
      r <- runif(1, min_r, tissue_radius)
      positions[i, ] <- c(r * cos(angle), r * sin(angle))
      i <- i + 1
    }

  } else if (mode == "center") {
    # Same as source positions (for comparison)
    positions <- generate_source_positions(n, source_region_radius)
  }

  colnames(positions) <- c("x", "y")
  positions
}

#' Generate receiver positions with distance-dependent mixing
#'
#' @param n Number of positions
#' @param tissue_radius Tissue radius
#' @param source_region_radius Source region radius
#' @return Matrix of positions (n x 2)
#' @export
generate_receiver_positions <- function(n, tissue_radius, source_region_radius) {
  positions <- matrix(0, nrow = n, ncol = 2)

  # Exclude source region
  i <- 1
  while (i <= n) {
    angle <- runif(1) * 2 * pi
    r <- sqrt(runif(1)) * tissue_radius

    if (r > source_region_radius) {
      positions[i, ] <- c(r * cos(angle), r * sin(angle))
      i <- i + 1
    }
  }

  colnames(positions) <- c("x", "y")
  positions
}

#' Generate inert cell positions
#'
#' @param n Number of positions
#' @param tissue_radius Tissue radius
#' @return Matrix of positions (n x 2)
#' @export
generate_inert_positions <- function(n, tissue_radius) {
  angles <- runif(n) * 2 * pi
  radii <- sqrt(runif(n)) * tissue_radius

  positions <- cbind(
    x = radii * cos(angles),
    y = radii * sin(angles)
  )

  positions
}

#' Calculate response probability based on distance from sources
#'
#' @param receiver_positions Receiver positions (n x 2)
#' @param source_positions Source positions (n_sources x 2)
#' @param sigma_decay Decay parameter
#' @param p_at_source Probability at source
#' @param p_background Background probability
#' @return Vector of probabilities
#' @export
calculate_response_probability <- function(receiver_positions, source_positions,
                                            sigma_decay, p_at_source, p_background) {
  n_receivers <- nrow(receiver_positions)

  # Calculate minimum distance to any source
  min_distances <- rep(Inf, n_receivers)

  for (i in seq_len(nrow(source_positions))) {
    d <- sqrt(rowSums((receiver_positions - matrix(source_positions[i, ],
                                                    nrow = n_receivers,
                                                    ncol = 2,
                                                    byrow = TRUE))^2))
    min_distances <- pmin(min_distances, d)
  }

  # Gaussian decay
  decay <- exp(-0.5 * (min_distances / sigma_decay)^2)

  # Mix probabilities
  probs <- p_background + (p_at_source - p_background) * decay

  probs
}

# =============================================================================
# DIFFERENTIAL EXPRESSION FILTERING
# =============================================================================
#' Perform differential expression analysis between two conditions
#'
#' @param expr_matrix Expression matrix (genes x cells)
#' @param condition_labels Vector of condition labels for each cell
#' @param condition1 Name of first condition (e.g., "cytokine")
#' @param condition2 Name of second condition (e.g., "PBS")
#' @param min_pct Minimum percentage of cells expressing gene
#' @return Data frame with DE results
#' @export
perform_de_analysis <- function(expr_matrix, condition_labels,
                                 condition1, condition2, min_pct = 0.1) {

  genes <- rownames(expr_matrix)
  n_genes <- length(genes)

  cells1 <- which(condition_labels == condition1)
  cells2 <- which(condition_labels == condition2)

  if (length(cells1) == 0 || length(cells2) == 0) {
    stop("No cells found for one or both conditions")
  }

  results <- data.frame(
    gene = genes,
    mean1 = numeric(n_genes),
    mean2 = numeric(n_genes),
    log2fc = numeric(n_genes),
    pct1 = numeric(n_genes),
    pct2 = numeric(n_genes),
    pval = numeric(n_genes),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_genes)) {
    expr1 <- expr_matrix[i, cells1]
    expr2 <- expr_matrix[i, cells2]

    m1 <- mean(expr1)
    m2 <- mean(expr2)

    results$mean1[i] <- m1
    results$mean2[i] <- m2
    results$log2fc[i] <- log2((m1 + 1e-10) / (m2 + 1e-10))
    results$pct1[i] <- mean(expr1 > 0)
    results$pct2[i] <- mean(expr2 > 0)

    # Wilcoxon test
    if (results$pct1[i] >= min_pct || results$pct2[i] >= min_pct) {
      tryCatch({
        wt <- wilcox.test(expr1, expr2, exact = FALSE)
        results$pval[i] <- wt$p.value
      }, error = function(e) {
        results$pval[i] <- 1.0
      })
    } else {
      results$pval[i] <- 1.0
    }
  }

  # FDR correction
  results$padj <- p.adjust(results$pval, method = "BH")

  results
}

#' Filter genes based on DE results
#'
#' @param de_results DE results from perform_de_analysis()
#' @param log2fc_threshold Log2 fold change threshold
#' @param pval_threshold Adjusted p-value threshold
#' @param direction Direction: "up", "down", or "both"
#' @return Character vector of significant genes
#' @export
filter_de_genes <- function(de_results, log2fc_threshold = 0.5,
                             pval_threshold = 0.05, direction = "up") {

  sig_mask <- de_results$padj < pval_threshold

  if (direction == "up") {
    sig_mask <- sig_mask & (de_results$log2fc > log2fc_threshold)
  } else if (direction == "down") {
    sig_mask <- sig_mask & (de_results$log2fc < -log2fc_threshold)
  } else {
    sig_mask <- sig_mask & (abs(de_results$log2fc) > log2fc_threshold)
  }

  de_results$gene[sig_mask]
}

# =============================================================================
# I_ND COMPUTATION (using sigdiscov functions)
# =============================================================================

#' Compute genome-wide I_ND for Parse10M simulation
#'
#' @param expr_matrix Expression matrix (genes x cells)
#' @param coords Cell coordinates (n_cells x 2)
#' @param sender_idx Sender cell indices (1-based)
#' @param receiver_idx Receiver cell indices (1-based)
#' @param radii Vector of radii
#' @param annular_width Annular width
#' @param verbose Print progress
#' @return Data frame with I_ND results for each gene and radius
#' @export
compute_genomewide_ind <- function(expr_matrix, coords, sender_idx, receiver_idx,
                                    radii, annular_width = 50, verbose = TRUE) {

  genes <- rownames(expr_matrix)
  n_genes <- length(genes)
  n_radii <- length(radii)

  if (verbose) {
    message(sprintf("Computing I_ND for %d genes across %d radii...", n_genes, n_radii))
  }

  # CRITICAL: Global normalization
  expr_norm <- t(scale(t(expr_matrix)))  # Gene-wise z-score

  # Extract sender/receiver positions
  sender_pos <- coords[sender_idx, , drop = FALSE]
  receiver_pos <- coords[receiver_idx, , drop = FALSE]

  # Compute pairwise distances once
  all_pos <- rbind(sender_pos, receiver_pos)
  n_senders <- length(sender_idx)
  n_receivers <- length(receiver_idx)

  dists <- as.matrix(dist(all_pos))[1:n_senders, (n_senders + 1):(n_senders + n_receivers)]

  results <- vector("list", n_genes * n_radii)
  idx <- 1

  for (g in seq_len(n_genes)) {
    gene <- genes[g]

    # Get normalized expression
    z_sender <- expr_norm[gene, sender_idx]
    z_receiver <- expr_norm[gene, receiver_idx]

    for (r_idx in seq_along(radii)) {
      r <- radii[r_idx]

      # Build annular weight matrix
      inner <- max(0, r - annular_width)
      outer <- r
      W <- ((dists > inner) & (dists <= outer)) * 1.0

      total_connections <- sum(W)

      if (total_connections == 0) {
        I_ND <- 0.0
      } else {
        # Row-normalize
        row_sums <- rowSums(W)
        row_sums[row_sums == 0] <- 1.0
        W_tilde <- W / row_sums

        # Compute I_ND
        mean_neighbor_z <- as.vector(W_tilde %*% z_receiver)
        numerator <- sum(z_sender * mean_neighbor_z)
        norm_z_s <- sqrt(sum(z_sender^2))
        norm_W_z_r <- sqrt(sum(mean_neighbor_z^2))

        if (norm_z_s > 1e-10 && norm_W_z_r > 1e-10) {
          I_ND <- numerator / (norm_z_s * norm_W_z_r)
        } else {
          I_ND <- 0.0
        }
      }

      results[[idx]] <- data.frame(
        gene = gene,
        radius = r,
        I_ND = I_ND,
        n_connections = as.integer(total_connections),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }

    if (verbose && g %% 100 == 0) {
      message(sprintf("  Processed %d/%d genes", g, n_genes))
    }
  }

  do.call(rbind, results)
}

#' Compute delta I_ND between active and control conditions
#'
#' @param ind_active I_ND results for active condition
#' @param ind_control I_ND results for control condition
#' @return Data frame with delta I_ND
#' @export
compute_delta_ind <- function(ind_active, ind_control) {
  # Merge by gene and radius
  merged <- merge(ind_active, ind_control,
                  by = c("gene", "radius"),
                  suffixes = c("_active", "_control"))

  merged$delta_I_ND <- merged$I_ND_active - merged$I_ND_control

  merged
}

# =============================================================================
# EVALUATION METRICS
# =============================================================================

#' Compute ROC and PR metrics
#'
#' @param scores Numeric vector of scores (higher = more likely positive)
#' @param labels Binary vector of true labels (1 = positive)
#' @return List with ROC and PR curve data and metrics
#' @export
compute_roc_pr_metrics <- function(scores, labels) {
  if (length(scores) != length(labels)) {
    stop("scores and labels must have same length")
  }

  # Remove NA
  valid <- !is.na(scores) & !is.na(labels)
  scores <- scores[valid]
  labels <- labels[valid]

  if (length(unique(labels)) < 2) {
    warning("Only one class present in labels")
    return(list(
      auroc = NA,
      auprc = NA,
      roc_curve = NULL,
      pr_curve = NULL
    ))
  }

  # Sort by score descending
  ord <- order(scores, decreasing = TRUE)
  scores_sorted <- scores[ord]
  labels_sorted <- labels[ord]

  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)
  n_total <- length(labels)

  # Compute ROC curve
  tpr <- cumsum(labels_sorted == 1) / n_pos
  fpr <- cumsum(labels_sorted == 0) / n_neg

  # Add (0,0) point
  tpr <- c(0, tpr)
  fpr <- c(0, fpr)

  # AUROC (trapezoidal rule)
  auroc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)

  # Compute PR curve
  tp <- cumsum(labels_sorted == 1)
  fp <- cumsum(labels_sorted == 0)
  precision <- tp / (tp + fp)
  recall <- tp / n_pos

  # Handle edge cases
  precision[is.nan(precision)] <- 1

  # Add (0, 1) point for recall
  precision <- c(1, precision)
  recall <- c(0, recall)

  # AUPRC (trapezoidal rule)
  auprc <- sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)

  list(
    auroc = auroc,
    auprc = auprc,
    roc_curve = data.frame(fpr = fpr, tpr = tpr),
    pr_curve = data.frame(recall = recall, precision = precision)
  )
}

#' Evaluate I_ND results against known response genes
#'
#' @param ind_results I_ND results data frame
#' @param response_genes Vector of known response genes
#' @param score_column Column to use as score (default: "delta_I_ND" or "I_ND")
#' @return Evaluation metrics at each radius
#' @export
evaluate_ind_results <- function(ind_results, response_genes,
                                  score_column = NULL) {

  # Determine score column
  if (is.null(score_column)) {
    if ("delta_I_ND" %in% colnames(ind_results)) {
      score_column <- "delta_I_ND"
    } else {
      score_column <- "I_ND"
    }
  }

  radii <- unique(ind_results$radius)
  results <- vector("list", length(radii))

  for (i in seq_along(radii)) {
    r <- radii[i]
    subset <- ind_results[ind_results$radius == r, ]

    # Create labels
    labels <- as.integer(subset$gene %in% response_genes)
    scores <- subset[[score_column]]

    metrics <- compute_roc_pr_metrics(scores, labels)

    results[[i]] <- data.frame(
      radius = r,
      auroc = metrics$auroc,
      auprc = metrics$auprc,
      n_positive = sum(labels),
      n_total = length(labels),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

#' Run Parse10M simulation workflow (with synthetic expression)
#'
#' This function creates a synthetic spatial simulation based on Parse10M-style
#' parameters without requiring the actual h5ad file.
#'
#' @param config Configuration from parse10m_config()
#' @param expr_matrix Optional expression matrix (genes x cells). If NULL, generates synthetic.
#' @param verbose Print progress messages
#' @return List with simulation results
#' @export
run_parse10m_simulation <- function(config, expr_matrix = NULL, verbose = TRUE) {

  set.seed(config$random_seed)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message(paste(rep("=", 70), collapse = ""))
    message("PARSE10M SPATIAL SIMULATION")
    message(paste(rep("=", 70), collapse = ""))
  }

  # Generate spatial positions
  if (verbose) message("\n1. Generating spatial layout...")

  source_pos <- generate_source_positions(
    config$n_source_cells, config$source_region_radius
  )
  control_pos <- generate_control_positions(
    config$n_control_sources, config$tissue_radius,
    config$source_region_radius, mode = "random"
  )
  receiver_pos <- generate_receiver_positions(
    config$n_receiver_cells, config$tissue_radius,
    config$source_region_radius
  )
  inert_pos <- generate_inert_positions(
    config$n_inert_cells, config$tissue_radius
  )

  if (verbose) {
    message(sprintf("  Source cells: %d", config$n_source_cells))
    message(sprintf("  Control sources: %d", config$n_control_sources))
    message(sprintf("  Receiver cells: %d", config$n_receiver_cells))
    message(sprintf("  Inert cells: %d", config$n_inert_cells))
  }

  # Combine all positions
  all_positions <- rbind(source_pos, control_pos, receiver_pos, inert_pos)
  n_total <- nrow(all_positions)

  # Define indices (1-based)
  source_idx <- 1:config$n_source_cells
  control_idx <- (config$n_source_cells + 1):(config$n_source_cells + config$n_control_sources)
  receiver_idx <- (config$n_source_cells + config$n_control_sources + 1):
                  (config$n_source_cells + config$n_control_sources + config$n_receiver_cells)
  inert_idx <- (n_total - config$n_inert_cells + 1):n_total

  # Calculate response probabilities
  response_probs <- calculate_response_probability(
    receiver_pos, source_pos,
    config$sigma_decay, config$p_cytokine_at_source, config$p_cytokine_background
  )

  # Generate or use provided expression matrix
  if (is.null(expr_matrix)) {
    if (verbose) message("\n2. Generating synthetic expression...")

    # Create synthetic expression for response genes + background
    n_response <- length(config$response_genes)
    n_background <- 500  # Additional background genes
    all_genes <- c(config$response_genes, paste0("BG_", seq_len(n_background)))

    expr_matrix <- matrix(
      rlnorm(length(all_genes) * n_total, 0, 1),
      nrow = length(all_genes),
      ncol = n_total,
      dimnames = list(all_genes, paste0("cell_", seq_len(n_total)))
    )

    # Add spatial signal for response genes in receivers
    for (g in config$response_genes) {
      if (g %in% rownames(expr_matrix)) {
        responding <- runif(config$n_receiver_cells) < response_probs
        expr_matrix[g, receiver_idx[responding]] <-
          expr_matrix[g, receiver_idx[responding]] * 5  # Fold change
      }
    }

    if (verbose) {
      message(sprintf("  Generated %d genes x %d cells", nrow(expr_matrix), ncol(expr_matrix)))
    }
  }

  # Compute I_ND for active (source) condition
  if (verbose) message("\n3. Computing I_ND for active condition...")

  ind_active <- compute_genomewide_ind(
    expr_matrix, all_positions,
    source_idx, receiver_idx,
    config$radii, config$annular_width, verbose
  )
  ind_active$condition <- "active"

  # Compute I_ND for control condition
  if (verbose) message("\n4. Computing I_ND for control condition...")

  ind_control <- compute_genomewide_ind(
    expr_matrix, all_positions,
    control_idx, receiver_idx,
    config$radii, config$annular_width, verbose
  )
  ind_control$condition <- "control"

  # Compute delta I_ND
  if (verbose) message("\n5. Computing delta I_ND...")

  delta_ind <- compute_delta_ind(ind_active, ind_control)

  # Evaluate against known response genes
  if (verbose) message("\n6. Evaluating against known response genes...")

  eval_results <- evaluate_ind_results(delta_ind, config$response_genes)

  # Find best radius
  best_idx <- which.max(eval_results$auroc)
  best_radius <- eval_results$radius[best_idx]
  best_auroc <- eval_results$auroc[best_idx]

  if (verbose) {
    message(sprintf("  Best AUROC: %.3f at radius %.0f", best_auroc, best_radius))
  }

  # Save results
  if (verbose) message("\n7. Saving results...")

  write.csv(ind_active, file.path(config$output_dir, "ind_active.csv"), row.names = FALSE)
  write.csv(ind_control, file.path(config$output_dir, "ind_control.csv"), row.names = FALSE)
  write.csv(delta_ind, file.path(config$output_dir, "delta_ind.csv"), row.names = FALSE)
  write.csv(eval_results, file.path(config$output_dir, "evaluation_metrics.csv"), row.names = FALSE)

  # Save positions
  positions_df <- data.frame(
    cell_id = paste0("cell_", seq_len(n_total)),
    x = all_positions[, 1],
    y = all_positions[, 2],
    cell_type = c(
      rep("source", config$n_source_cells),
      rep("control_source", config$n_control_sources),
      rep("receiver", config$n_receiver_cells),
      rep("inert", config$n_inert_cells)
    )
  )
  write.csv(positions_df, file.path(config$output_dir, "cell_positions.csv"), row.names = FALSE)

  # Save config
  saveRDS(config, file.path(config$output_dir, "config.rds"))

  if (verbose) {
    message(sprintf("\nSimulation complete. Results in %s", config$output_dir))
  }

  list(
    positions = all_positions,
    positions_df = positions_df,
    source_idx = source_idx,
    control_idx = control_idx,
    receiver_idx = receiver_idx,
    inert_idx = inert_idx,
    response_probs = response_probs,
    expr_matrix = expr_matrix,
    ind_active = ind_active,
    ind_control = ind_control,
    delta_ind = delta_ind,
    eval_results = eval_results,
    best_radius = best_radius,
    best_auroc = best_auroc,
    config = config
  )
}

# =============================================================================
# VISUALIZATION
# =============================================================================

#' Plot spatial layout
#'
#' @param results Results from run_parse10m_simulation()
#' @param output_path Output file path (PNG)
#' @export
plot_parse10m_layout <- function(results, output_path = NULL) {
  pos <- results$positions
  pos_df <- results$positions_df

  if (!is.null(output_path)) {
    png(output_path, width = 10, height = 10, units = "in", res = 300)
  }

  # Color mapping
  colors <- c(
    "source" = "red",
    "control_source" = "orange",
    "receiver" = "blue",
    "inert" = "lightgray"
  )

  plot(NULL, xlim = range(pos[, 1]) * 1.1, ylim = range(pos[, 2]) * 1.1,
       xlab = "X (um)", ylab = "Y (um)", main = "Parse10M Spatial Layout",
       asp = 1)

  # Plot in order: inert, receiver, control, source
  for (ct in c("inert", "receiver", "control_source", "source")) {
    idx <- which(pos_df$cell_type == ct)
    if (length(idx) > 0) {
      cex <- ifelse(ct %in% c("source", "control_source"), 1, 0.2)
      pch <- ifelse(ct %in% c("source", "control_source"), 17, 16)
      points(pos[idx, 1], pos[idx, 2],
             col = colors[ct], pch = pch, cex = cex)
    }
  }

  legend("topright",
         legend = names(colors),
         col = colors,
         pch = c(17, 17, 16, 16),
         pt.cex = c(1, 1, 0.5, 0.5))

  if (!is.null(output_path)) {
    dev.off()
    message(sprintf("Saved plot to %s", output_path))
  }
}

#' Plot AUROC/AUPRC across radii
#'
#' @param eval_results Evaluation results from run_parse10m_simulation()
#' @param output_path Output file path (PNG)
#' @export
plot_evaluation_curves <- function(eval_results, output_path = NULL) {
  if (!is.null(output_path)) {
    png(output_path, width = 10, height = 6, units = "in", res = 300)
  }

  par(mfrow = c(1, 2))

  # AUROC
  plot(eval_results$radius, eval_results$auroc,
       type = "l", lwd = 2, col = "blue",
       xlab = "Radius (um)", ylab = "AUROC",
       main = "AUROC vs Radius", ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "gray")
  grid()

  best_idx <- which.max(eval_results$auroc)
  points(eval_results$radius[best_idx], eval_results$auroc[best_idx],
         pch = 19, col = "red", cex = 1.5)
  text(eval_results$radius[best_idx], eval_results$auroc[best_idx],
       sprintf("%.3f", eval_results$auroc[best_idx]),
       pos = 3, col = "red")

  # AUPRC
  plot(eval_results$radius, eval_results$auprc,
       type = "l", lwd = 2, col = "darkgreen",
       xlab = "Radius (um)", ylab = "AUPRC",
       main = "AUPRC vs Radius", ylim = c(0, 1))
  grid()

  best_idx <- which.max(eval_results$auprc)
  points(eval_results$radius[best_idx], eval_results$auprc[best_idx],
         pch = 19, col = "red", cex = 1.5)
  text(eval_results$radius[best_idx], eval_results$auprc[best_idx],
       sprintf("%.3f", eval_results$auprc[best_idx]),
       pos = 3, col = "red")

  par(mfrow = c(1, 1))

  if (!is.null(output_path)) {
    dev.off()
    message(sprintf("Saved plot to %s", output_path))
  }
}

#' Plot delta I_ND profile
#'
#' @param delta_ind Delta I_ND results
#' @param response_genes Known response genes
#' @param output_path Output file path (PNG)
#' @export
plot_delta_ind_profile <- function(delta_ind, response_genes, output_path = NULL) {
  if (!is.null(output_path)) {
    png(output_path, width = 12, height = 6, units = "in", res = 300)
  }

  radii <- unique(delta_ind$radius)

  # Calculate mean delta I_ND for response vs background genes
  response_means <- sapply(radii, function(r) {
    subset <- delta_ind[delta_ind$radius == r & delta_ind$gene %in% response_genes, ]
    mean(subset$delta_I_ND, na.rm = TRUE)
  })

  background_means <- sapply(radii, function(r) {
    subset <- delta_ind[delta_ind$radius == r & !(delta_ind$gene %in% response_genes), ]
    mean(subset$delta_I_ND, na.rm = TRUE)
  })

  plot(radii, response_means, type = "l", lwd = 2, col = "red",
       xlab = "Radius (um)", ylab = "Mean Delta I_ND",
       main = "Delta I_ND: Response vs Background Genes",
       ylim = range(c(response_means, background_means), na.rm = TRUE))
  lines(radii, background_means, lwd = 2, col = "gray50")
  abline(h = 0, lty = 2)
  grid()

  legend("topright",
         legend = c("Response genes", "Background genes"),
         col = c("red", "gray50"),
         lwd = 2)

  if (!is.null(output_path)) {
    dev.off()
    message(sprintf("Saved plot to %s", output_path))
  }
}
