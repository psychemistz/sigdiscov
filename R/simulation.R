#' Unified Spatial Simulation Framework for sigdiscov
#'
#' This module provides a complete simulation framework for generating
#' spatial transcriptomics data with sender-receiver cell interactions.
#' It replaces the Python unified_simulation.py with native R implementation.
#'
#' @name simulation
#' @keywords internal
NULL

# =============================================================================
# CONFIGURATION CLASSES (using lists as R doesn't have native dataclasses)
# =============================================================================

#' Create domain configuration
#'
#' @param n_cells Number of cells in simulation
#' @param max_radius Maximum radius of circular domain
#' @param center Center coordinates (x, y)
#' @param random_seed Random seed for reproducibility
#' @return Domain configuration list
#' @export
domain_config <- function(n_cells = 100000,
                          max_radius = 5000.0,
                          center = c(0.0, 0.0),
                          random_seed = 42) {
  list(
    n_cells = n_cells,
    max_radius = max_radius,
    center = center,
    random_seed = random_seed
  )
}

#' Create cell type configuration
#'
#' @param n_active_senders Number of active sender cells
#' @param n_silent_senders Number of silent sender cells
#' @param receiver_fractions Vector of receiver fractions to simulate
#' @return Cell type configuration list
#' @export
celltype_config <- function(n_active_senders = 20,
                            n_silent_senders = 0,
                            receiver_fractions = c(0.1, 0.2, 0.3, 0.5)) {
  list(
    n_active_senders = n_active_senders,
    n_silent_senders = n_silent_senders,
    receiver_fractions = receiver_fractions
  )
}

#' Create sender position configuration
#'
#' @param mode Position mode: "center", "fixed_5", or "random"
#' @param n_positions Number of positions for random mode
#' @param offset_distance Offset distance for fixed_5 mode
#' @param min_separation Minimum separation for random mode
#' @return Position configuration list
#' @export
position_config <- function(mode = "center",
                            n_positions = 1,
                            offset_distance = 3000.0,
                            min_separation = 500.0) {
  list(
    mode = mode,
    n_positions = n_positions,
    offset_distance = offset_distance,
    min_separation = min_separation
  )
}

#' Create diffusion configuration
#'
#' @param D Diffusion coefficient (um^2/s)
#' @param k_max Maximum uptake rate
#' @param Kd Dissociation constant
#' @param secretion_rate Expression to concentration scaling
#' @return Diffusion configuration list
#' @export
diffusion_config <- function(D = 100.0,
                             k_max = 10.0,
                             Kd = 1.0,
                             secretion_rate = 1.0) {
  list(
    D = D,
    k_max = k_max,
    Kd = Kd,
    secretion_rate = secretion_rate
  )
}

#' Create expression configuration
#'
#' @param F_basal Basal factor expression
#' @param F_high High factor expression for active senders
#' @param R_basal Basal response expression
#' @param fold_change Response fold change
#' @param sigma_f Lognormal sigma for expressing senders
#' @param sigma_f_basal Lognormal sigma for non-expressing senders
#' @param sigma_r Lognormal sigma for responding receivers
#' @param sigma_r_basal Lognormal sigma for non-responding receivers
#' @return Expression configuration list
#' @export
expression_config <- function(F_basal = 0.1,
                              F_high = 10.0,
                              R_basal = 0.1,
                              fold_change = 10.0,
                              sigma_f = 0.1,
                              sigma_f_basal = 0.1,
                              sigma_r = 0.1,
                              sigma_r_basal = 0.1) {
  list(
    F_basal = F_basal,
    F_high = F_high,
    R_basal = R_basal,
    fold_change = fold_change,
    sigma_f = sigma_f,
    sigma_f_basal = sigma_f_basal,
    sigma_r = sigma_r,
    sigma_r_basal = sigma_r_basal
  )
}

#' Create stochastic configuration
#'
#' @param p_sender_express Probability sender is "on"
#' @param p_receiver_respond_max Max probability receiver responds
#' @param hill_coefficient Hill coefficient for response probability
#' @param zero_inflate_factor Fraction of non-sender cells with zero expression
#' @param zero_inflate_response Fraction of non-receiver cells with zero expression
#' @return Stochastic configuration list
#' @export
stochastic_config <- function(p_sender_express = 0.9,
                              p_receiver_respond_max = 1.0,
                              hill_coefficient = 1.0,
                              zero_inflate_factor = 0.0,
                              zero_inflate_response = 0.0) {
  list(
    p_sender_express = p_sender_express,
    p_receiver_respond_max = p_receiver_respond_max,
    hill_coefficient = hill_coefficient,
    zero_inflate_factor = zero_inflate_factor,
    zero_inflate_response = zero_inflate_response
  )
}

#' Create normalization configuration
#'
#' @param method Normalization method: "zscore", "vst_log1p", "vst_pearson", "vst_shifted_log"
#' @param vst_theta Theta for Pearson residuals
#' @param vst_clip_value Clip value for Pearson residuals
#' @param vst_pseudocount Pseudocount for shifted log
#' @return Normalization configuration list
#' @export
normalization_config <- function(method = "zscore",
                                 vst_theta = 100.0,
                                 vst_clip_value = 10.0,
                                 vst_pseudocount = 1.0) {
  list(
    method = method,
    vst_theta = vst_theta,
    vst_clip_value = vst_clip_value,
    vst_pseudocount = vst_pseudocount
  )
}

#' Create analysis configuration
#'
#' @param radii Vector of radii to evaluate
#' @param bandwidth Ring bandwidth
#' @param weight_type Weight type: "ring", "gaussian", or "annular"
#' @param annular_width Annular width
#' @return Analysis configuration list
#' @export
analysis_config <- function(radii = seq(50, 5000, by = 50),
                            bandwidth = 100.0,
                            weight_type = "ring",
                            annular_width = 50.0) {
  list(
    radii = radii,
    bandwidth = bandwidth,
    weight_type = weight_type,
    annular_width = annular_width
  )
}

#' Create complete simulation configuration
#'
#' @param domain Domain configuration (from domain_config())
#' @param cell_types Cell type configuration (from celltype_config())
#' @param position Position configuration (from position_config())
#' @param diffusion Diffusion configuration (from diffusion_config())
#' @param expression Expression configuration (from expression_config())
#' @param stochastic Stochastic configuration (from stochastic_config())
#' @param normalization Normalization configuration (from normalization_config())
#' @param analysis Analysis configuration (from analysis_config())
#' @param output_dir Output directory
#' @return Complete simulation configuration
#' @export
simulation_config <- function(domain = domain_config(),
                              cell_types = celltype_config(),
                              position = position_config(),
                              diffusion = diffusion_config(),
                              expression = expression_config(),
                              stochastic = stochastic_config(),
                              normalization = normalization_config(),
                              analysis = analysis_config(),
                              output_dir = "./output") {
  list(
    domain = domain,
    cell_types = cell_types,
    position = position,
    diffusion = diffusion,
    expression = expression,
    stochastic = stochastic,
    normalization = normalization,
    analysis = analysis,
    output_dir = output_dir
  )
}

# =============================================================================
# PRESET CONFIGURATIONS
# =============================================================================

#' Get preset simulation configuration
#'
#' Available presets:
#' - "basic": Basic simulation with senders at center
#' - "stochastic_sender": Stochastic sender expression
#' - "stochastic_hill": Stochastic Hill-function receiver response
#' - "fixed_5": Senders at 5 fixed positions
#' - "random_multi": Senders at multiple random positions
#' - "full_stochastic": Full stochastic model
#' - "vst": VST-normalized simulation
#' - "annular": Annular weight matrices
#'
#' @param name Preset name
#' @return Simulation configuration
#' @export
get_simulation_preset <- function(name) {
  config <- simulation_config()

  if (name == "basic") {
    config$position$mode <- "center"
    config$stochastic$p_sender_express <- 1.0
    config$stochastic$p_receiver_respond_max <- 1.0

  } else if (name == "stochastic_sender") {
    config$stochastic$p_sender_express <- 0.7
    config$stochastic$p_receiver_respond_max <- 1.0

  } else if (name == "stochastic_hill") {
    config$stochastic$p_sender_express <- 0.9
    config$stochastic$p_receiver_respond_max <- 0.9
    config$stochastic$hill_coefficient <- 1.0

  } else if (name == "fixed_5") {
    config$position$mode <- "fixed_5"
    config$position$offset_distance <- 3000.0
    config$cell_types$n_active_senders <- 100

  } else if (name == "random_multi") {
    config$position$mode <- "random"
    config$position$n_positions <- 5
    config$position$min_separation <- 500.0

  } else if (name == "full_stochastic") {
    config$position$mode <- "random"
    config$position$n_positions <- 1
    config$cell_types$n_active_senders <- 200
    config$cell_types$n_silent_senders <- 100
    config$stochastic$p_sender_express <- 0.9
    config$stochastic$p_receiver_respond_max <- 1.0
    config$stochastic$hill_coefficient <- 1.0
    config$analysis$bandwidth <- 20.0

  } else if (name == "vst") {
    config$normalization$method <- "vst_log1p"
    config$expression$F_basal <- 0.5
    config$expression$F_high <- 50.0
    config$expression$R_basal <- 0.5
    config$expression$fold_change <- 5.0
    config$expression$sigma_f <- 0.8
    config$expression$sigma_r <- 0.8
    config$stochastic$zero_inflate_factor <- 0.7
    config$stochastic$zero_inflate_response <- 0.5

  } else if (name == "annular") {
    config$analysis$weight_type <- "annular"
    config$analysis$annular_width <- 50.0

  } else {
    stop(paste("Unknown preset:", name,
               "\nAvailable: basic, stochastic_sender, stochastic_hill, fixed_5,",
               "random_multi, full_stochastic, vst, annular"))
  }

  config
}

# =============================================================================
# SPATIAL DOMAIN GENERATION
# =============================================================================

#' Generate cell positions in circular domain
#'
#' @param config Domain configuration
#' @return Matrix of cell positions (n_cells x 2)
#' @export
generate_cell_positions <- function(config) {
  set.seed(config$random_seed)

  n <- config$n_cells
  center <- config$center
  max_r <- config$max_radius

  # Uniform distribution in circle using polar coordinates
  angles <- runif(n) * 2 * pi
  radii <- sqrt(runif(n)) * max_r

  positions <- cbind(
    x = center[1] + radii * cos(angles),
    y = center[2] + radii * sin(angles)
  )

  positions
}

#' Generate sender position dictionary
#'
#' @param config Position configuration
#' @param domain_config Domain configuration
#' @return Named list of position coordinates
#' @export
generate_sender_positions <- function(config, domain_config) {
  center <- domain_config$center

  if (config$mode == "center") {
    return(list(C = center))

  } else if (config$mode == "fixed_5") {
    offset <- config$offset_distance
    return(list(
      C = center,
      W = c(center[1] - offset, center[2]),
      E = c(center[1] + offset, center[2]),
      N = c(center[1], center[2] + offset),
      S = c(center[1], center[2] - offset)
    ))

  } else if (config$mode == "random") {
    return(.generate_random_positions(config, domain_config))

  } else {
    stop(paste("Unknown position mode:", config$mode))
  }
}

#' Generate random positions with minimum separation
#' @keywords internal
.generate_random_positions <- function(config, domain_config) {
  center <- domain_config$center
  effective_radius <- domain_config$max_radius * 0.8

  positions <- list()
  coords_list <- list()

  attempts <- 0
  max_attempts <- 1000

  while (length(positions) < config$n_positions && attempts < max_attempts) {
    angle <- runif(1) * 2 * pi
    r <- sqrt(runif(1)) * effective_radius

    new_pos <- c(
      center[1] + r * cos(angle),
      center[2] + r * sin(angle)
    )

    if (length(coords_list) == 0) {
      valid <- TRUE
    } else {
      distances <- sapply(coords_list, function(existing) {
        sqrt(sum((new_pos - existing)^2))
      })
      valid <- min(distances) >= config$min_separation
    }

    if (valid) {
      pos_label <- paste0("P", length(positions) + 1)
      positions[[pos_label]] <- new_pos
      coords_list[[length(coords_list) + 1]] <- new_pos
    }

    attempts <- attempts + 1
  }

  if (length(positions) < config$n_positions) {
    warning(sprintf("Only placed %d/%d positions", length(positions), config$n_positions))
  }

  positions
}

#' Distribute senders across positions
#'
#' @param n_senders Number of senders to distribute
#' @param position_dict Named list of positions
#' @return List with assignments and sender_positions
#' @export
distribute_senders <- function(n_senders, position_dict) {
  position_names <- names(position_dict)
  n_pos <- length(position_names)

  if (n_senders < n_pos) {
    stop(sprintf("n_senders (%d) must be >= number of positions (%d)", n_senders, n_pos))
  }

  # Each position gets at least one sender
  assignments <- setNames(rep(1, n_pos), position_names)
  sender_positions <- lapply(position_names, function(name) {
    list(name = name, coord = position_dict[[name]])
  })

  # Distribute remaining senders randomly
  remaining <- n_senders - n_pos
  for (i in seq_len(remaining)) {
    chosen <- sample(position_names, 1)
    assignments[chosen] <- assignments[chosen] + 1
    sender_positions[[length(sender_positions) + 1]] <- list(
      name = chosen,
      coord = position_dict[[chosen]]
    )
  }

  list(
    assignments = assignments,
    sender_positions = sender_positions
  )
}

# =============================================================================
# DIFFUSION SOLVER
# =============================================================================

#' Calculate characteristic diffusion length
#'
#' @param diffusion_config Diffusion configuration
#' @param n_receivers Receiver density
#' @param p_r Probability of responding
#' @return Lambda (characteristic diffusion length)
#' @export
calculate_lambda <- function(diffusion_config, n_receivers, p_r = 1.0) {
  n_eff <- n_receivers * p_r
  if (n_eff <= 0) return(Inf)
  sqrt(diffusion_config$D * diffusion_config$Kd / (n_eff * diffusion_config$k_max))
}

#' Solve concentration field from sender cells
#'
#' @param sender_positions Sender cell positions (n_senders x 2)
#' @param sender_expression Sender expression values
#' @param cell_positions All cell positions (n_cells x 2)
#' @param diffusion_config Diffusion configuration
#' @param n_receivers_density Receiver density for lambda calculation
#' @param position_dict Optional position dictionary for grouping
#' @param p_r Probability of responding
#' @return List with concentrations and lambda
#' @export
solve_diffusion <- function(sender_positions, sender_expression,
                            cell_positions, diffusion_config,
                            n_receivers_density, position_dict = NULL,
                            p_r = 1.0) {

  lambda_val <- calculate_lambda(diffusion_config, n_receivers_density, p_r)
  n_cells <- nrow(cell_positions)
  concentrations <- rep(0, n_cells)

  # Identify active senders (expression > basal threshold)
  active_mask <- sender_expression > 0.2
  active_expr <- sender_expression[active_mask]
  active_pos <- sender_positions[active_mask, , drop = FALSE]

  if (sum(active_mask) == 0) {
    return(list(concentrations = concentrations, lambda = lambda_val))
  }

  # Group by unique positions if using position_dict
  if (!is.null(position_dict) && length(position_dict) > 1) {
    # Group active senders by their position
    unique_positions <- list()

    for (i in seq_len(nrow(active_pos))) {
      pos_key <- paste(round(active_pos[i, ], 6), collapse = "_")
      if (is.null(unique_positions[[pos_key]])) {
        unique_positions[[pos_key]] <- list(pos = active_pos[i, ], total_expr = 0)
      }
      unique_positions[[pos_key]]$total_expr <-
        unique_positions[[pos_key]]$total_expr + active_expr[i]
    }

    for (pos_key in names(unique_positions)) {
      source_pos <- unique_positions[[pos_key]]$pos
      total_factor <- unique_positions[[pos_key]]$total_expr * diffusion_config$secretion_rate

      for (i in seq_len(n_cells)) {
        cell_pos <- cell_positions[i, ]
        r <- sqrt(sum((cell_pos - source_pos)^2))

        if (r < 1e-3) {
          concentrations[i] <- concentrations[i] + total_factor * 100
        } else {
          concentrations[i] <- concentrations[i] + total_factor * exp(-r / lambda_val) / sqrt(r)
        }
      }
    }
  } else {
    # Single source at center
    center_pos <- active_pos[1, ]
    total_factor <- sum(active_expr) * diffusion_config$secretion_rate

    for (i in seq_len(n_cells)) {
      r <- sqrt(sum((cell_positions[i, ] - center_pos)^2))

      if (r < 1e-3) {
        concentrations[i] <- total_factor * 100
      } else {
        concentrations[i] <- total_factor * exp(-r / lambda_val) / sqrt(r)
      }
    }
  }

  list(concentrations = concentrations, lambda = lambda_val)
}

# =============================================================================
# EXPRESSION GENERATION
# =============================================================================

#' Generate stochastic factor expression
#'
#' Model:
#'   F_i = S_i * F_high * LogNormal(0, sigma_f^2) + (1-S_i) * F_basal * LogNormal(0, sigma_f_basal^2)
#'   S_i ~ Bernoulli(p_sender_express)
#'
#' @param n_total Total number of cells
#' @param n_active Number of active senders
#' @param active_indices Indices of active sender cells
#' @param expr_config Expression configuration
#' @param stoch_config Stochastic configuration
#' @return List with factor_expr and expressing_mask
#' @export
generate_factor_expression <- function(n_total, n_active, active_indices,
                                        expr_config, stoch_config) {

  # Initialize with basal expression
  factor_expr <- expr_config$F_basal * rlnorm(n_total, 0, expr_config$sigma_f_basal)

  # Zero-inflation for non-senders
  if (stoch_config$zero_inflate_factor > 0) {
    zero_mask <- runif(n_total) < stoch_config$zero_inflate_factor
    zero_mask[active_indices] <- FALSE
    factor_expr[zero_mask] <- 0.0
  }

  # Stochastic sender expression
  expressing_mask <- runif(n_active) < stoch_config$p_sender_express
  n_expressing <- sum(expressing_mask)

  if (n_expressing > 0) {
    expressing_indices <- active_indices[expressing_mask]
    factor_expr[expressing_indices] <-
      expr_config$F_high * rlnorm(n_expressing, 0, expr_config$sigma_f)
  }

  list(
    factor_expr = factor_expr,
    expressing_mask = expressing_mask
  )
}

#' Generate stochastic response expression
#'
#' Hybrid Model:
#'   1. B_j ~ Bernoulli(p_max * C^n / (Kd^n + C^n))
#'   2. Act_j = C / (Kd + C)
#'   3. If responding: R_j = R_basal * (1 + FC * Act_j) * LogNormal(0, sigma_r^2)
#'      If not: R_j = R_basal * LogNormal(0, sigma_r_basal^2)
#'
#' @param n_total Total number of cells
#' @param receiver_indices Indices of receiver cells
#' @param concentrations Concentration values at all cells
#' @param expr_config Expression configuration
#' @param stoch_config Stochastic configuration
#' @return List with response_expr, responding_mask, response_probs
#' @export
generate_response_expression <- function(n_total, receiver_indices, concentrations,
                                          expr_config, stoch_config) {

  n_receivers <- length(receiver_indices)

  # Initialize with basal expression
  response_expr <- expr_config$R_basal * rlnorm(n_total, 0, expr_config$sigma_r_basal)

  # Zero-inflation for non-receivers
  if (stoch_config$zero_inflate_response > 0) {
    zero_mask <- runif(n_total) < stoch_config$zero_inflate_response
    zero_mask[receiver_indices] <- FALSE
    response_expr[zero_mask] <- 0.0
  }

  # Get concentrations at receiver positions
  C <- concentrations[receiver_indices]

  # Calculate response probability (Hill function)
  n <- stoch_config$hill_coefficient
  Kd <- 1.0  # Normalized Kd for response probability
  C_n <- C^n
  Kd_n <- Kd^n
  response_probs <- stoch_config$p_receiver_respond_max * C_n / (Kd_n + C_n)

  # Stochastic binary decision
  responding_mask <- runif(n_receivers) < response_probs

  # Non-responding receivers
  non_responding_indices <- receiver_indices[!responding_mask]
  n_non_responding <- length(non_responding_indices)
  if (n_non_responding > 0) {
    response_expr[non_responding_indices] <-
      expr_config$R_basal * rlnorm(n_non_responding, 0, expr_config$sigma_r_basal)
  }

  # Responding receivers
  responding_indices <- receiver_indices[responding_mask]
  n_responding <- length(responding_indices)
  if (n_responding > 0) {
    C_responding <- C[responding_mask]
    activation <- C_responding / (Kd + C_responding)
    mean_expr <- expr_config$R_basal * (1 + expr_config$fold_change * activation)
    response_expr[responding_indices] <-
      mean_expr * rlnorm(n_responding, 0, expr_config$sigma_r)
  }

  list(
    response_expr = response_expr,
    responding_mask = responding_mask,
    response_probs = response_probs
  )
}

# =============================================================================
# NORMALIZATION
# =============================================================================

#' Normalize expression data
#'
#' @param raw_expr Raw expression vector
#' @param norm_config Normalization configuration
#' @return Normalized expression vector
#' @export
normalize_expression <- function(raw_expr, norm_config) {
  method <- norm_config$method

  if (method == "zscore") {
    mu <- mean(raw_expr)
    sigma <- sd(raw_expr)
    return((raw_expr - mu) / (sigma + 1e-10))

  } else if (method == "vst_log1p") {
    log_expr <- log1p(raw_expr)
    return(log_expr - mean(log_expr))

  } else if (method == "vst_pearson") {
    mu <- mean(raw_expr) + 1e-10
    variance <- mu + (mu^2) / norm_config$vst_theta
    residuals <- (raw_expr - mu) / sqrt(variance)
    return(pmax(pmin(residuals, norm_config$vst_clip_value), -norm_config$vst_clip_value))

  } else if (method == "vst_shifted_log") {
    log_expr <- log(raw_expr + norm_config$vst_pseudocount)
    baseline <- log(norm_config$vst_pseudocount)
    shift <- baseline + 0.5
    return(log_expr - shift)

  } else {
    stop(paste("Unknown normalization method:", method))
  }
}

# =============================================================================
# I_ND COMPUTATION
# =============================================================================

#' Compute I_ND (Normalized Directional Moran's I) at a single radius
#'
#' Formula: I_ND = (z_f . W.z_r) / (||z_f|| . ||W.z_r||)
#'
#' @param sender_indices Indices of sender cells
#' @param receiver_indices Indices of receiver cells
#' @param all_positions All cell positions (n_cells x 2)
#' @param factor_expr Factor expression (all cells)
#' @param response_expr Response expression (all cells)
#' @param radius Radius for weight matrix
#' @param analysis_config Analysis configuration
#' @return List with I_ND and n_connections
#' @export
compute_ind_single <- function(sender_indices, receiver_indices, all_positions,
                                factor_expr, response_expr, radius, analysis_config) {

  # Global normalization (CRITICAL)
  mu_f <- mean(factor_expr)
  sigma_f <- sd(factor_expr)
  mu_r <- mean(response_expr)
  sigma_r <- sd(response_expr)

  z_s <- (factor_expr[sender_indices] - mu_f) / (sigma_f + 1e-10)
  z_r <- (response_expr[receiver_indices] - mu_r) / (sigma_r + 1e-10)

  sender_pos <- all_positions[sender_indices, , drop = FALSE]
  receiver_pos <- all_positions[receiver_indices, , drop = FALSE]
  n_senders <- length(sender_indices)
  n_receivers <- length(receiver_indices)

  # Compute cross-distance matrix directly (memory efficient)
  dists <- matrix(0, nrow = n_senders, ncol = n_receivers)
  for (i in seq_len(n_senders)) {
    dists[i, ] <- sqrt(rowSums((receiver_pos - matrix(sender_pos[i, ],
                                                       nrow = n_receivers,
                                                       ncol = 2,
                                                       byrow = TRUE))^2))
  }

  # Build weight matrix
  W <- .build_weight_matrix(dists, radius, analysis_config)

  # Row-normalize
  row_sums <- rowSums(W)
  total_connections <- sum(W)

  if (total_connections == 0) {
    return(list(I_ND = 0.0, n_connections = 0))
  }

  row_sums[row_sums == 0] <- 1.0
  W_tilde <- W / row_sums

  # Compute I_ND
  mean_neighbor_z <- as.vector(W_tilde %*% z_r)
  numerator <- sum(z_s * mean_neighbor_z)
  norm_z_s <- sqrt(sum(z_s^2))
  norm_W_z_r <- sqrt(sum(mean_neighbor_z^2))

  if (norm_z_s > 1e-10 && norm_W_z_r > 1e-10) {
    I_ND <- numerator / (norm_z_s * norm_W_z_r)
  } else {
    I_ND <- 0.0
  }

  list(I_ND = I_ND, n_connections = as.integer(total_connections))
}

#' Build weight matrix based on weight type
#' @keywords internal
.build_weight_matrix <- function(dists, radius, analysis_config) {
  weight_type <- analysis_config$weight_type

  if (weight_type == "ring") {
    half_bw <- analysis_config$bandwidth / 2.0
    lower <- radius - half_bw
    upper <- radius + half_bw
    return(((dists > lower) & (dists <= upper)) * 1.0)

  } else if (weight_type == "annular") {
    inner <- max(0, radius - analysis_config$annular_width)
    outer <- radius
    return(((dists > inner) & (dists <= outer)) * 1.0)

  } else if (weight_type == "gaussian") {
    sigma <- analysis_config$bandwidth
    return(exp(-0.5 * (dists / sigma)^2))

  } else {
    stop(paste("Unknown weight type:", weight_type))
  }
}

#' Compute I_ND across multiple radii
#'
#' @param sender_indices Indices of sender cells
#' @param receiver_indices Indices of receiver cells
#' @param all_positions All cell positions (n_cells x 2)
#' @param factor_expr Factor expression (all cells)
#' @param response_expr Response expression (all cells)
#' @param analysis_config Analysis configuration
#' @param verbose Print progress
#' @return Data frame with radius, I_ND, n_connections
#' @export
compute_ind_radial <- function(sender_indices, receiver_indices, all_positions,
                                factor_expr, response_expr, analysis_config,
                                verbose = TRUE) {

  radii <- analysis_config$radii
  n_radii <- length(radii)

  if (verbose) {
    message(sprintf("Computing I_ND across %d radii...", n_radii))
  }

  results <- data.frame(
    radius = numeric(n_radii),
    I_ND = numeric(n_radii),
    n_connections = integer(n_radii)
  )

  for (i in seq_along(radii)) {
    r <- radii[i]

    result <- compute_ind_single(
      sender_indices, receiver_indices, all_positions,
      factor_expr, response_expr, r, analysis_config
    )

    results$radius[i] <- r
    results$I_ND[i] <- result$I_ND
    results$n_connections[i] <- result$n_connections

    if (verbose && i %% 10 == 0) {
      message(sprintf("  Processed radius %d/%d (r=%.0f, I_ND=%.4f)",
                      i, n_radii, r, result$I_ND))
    }
  }

  results
}

# =============================================================================
# MAIN SIMULATION FUNCTION
# =============================================================================

#' Run unified simulation
#'
#' @param config Simulation configuration (from simulation_config() or get_simulation_preset())
#' @param verbose Print progress messages
#' @return List with results for each receiver fraction
#' @export
run_simulation <- function(config, verbose = TRUE) {

  if (verbose) {
    message(paste(rep("=", 70), collapse = ""))
    message("UNIFIED SIMULATION")
    message(paste(rep("=", 70), collapse = ""))
  }

  set.seed(config$domain$random_seed)

  # Create output directory
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Generate positions
  all_positions <- generate_cell_positions(config$domain)
  position_dict <- generate_sender_positions(config$position, config$domain)

  if (verbose) {
    message(sprintf("Generated %d cell positions", nrow(all_positions)))
    message(sprintf("Sender positions: %s", paste(names(position_dict), collapse = ", ")))
  }

  # Assign senders
  n_senders <- config$cell_types$n_active_senders + config$cell_types$n_silent_senders
  n_active <- config$cell_types$n_active_senders

  all_indices <- seq_len(nrow(all_positions))
  sender_indices <- sample(all_indices, n_senders, replace = FALSE)
  active_indices <- sender_indices[1:n_active]
  silent_indices <- if (n_senders > n_active) sender_indices[(n_active + 1):n_senders] else integer(0)

  # Distribute active senders
  dist_result <- distribute_senders(n_active, position_dict)
  assignments <- dist_result$assignments
  sender_pos_list <- dist_result$sender_positions

  # Set sender positions
  for (i in seq_along(sender_pos_list)) {
    all_positions[active_indices[i], ] <- sender_pos_list[[i]]$coord
  }

  if (verbose) {
    message(sprintf("Active sender distribution: %s",
                    paste(names(assignments), assignments, sep = ":", collapse = ", ")))
  }

  # Save config
  config_list <- config
  saveRDS(config_list, file.path(config$output_dir, "config.rds"))

  results <- list()

  for (frac in config$cell_types$receiver_fractions) {
    if (verbose) {
      message(sprintf("\nProcessing %.0f%% Receivers...", frac * 100))
    }

    # Assign receivers
    non_sender_indices <- setdiff(all_indices, sender_indices)
    n_receivers <- as.integer(nrow(all_positions) * frac)
    receiver_indices <- sample(non_sender_indices, n_receivers, replace = FALSE)

    # Generate expression
    factor_result <- generate_factor_expression(
      nrow(all_positions), n_active, active_indices,
      config$expression, config$stochastic
    )
    factor_expr <- factor_result$factor_expr
    expressing_mask <- factor_result$expressing_mask

    n_expressing <- sum(expressing_mask)
    if (verbose) {
      message(sprintf("  Expressing senders: %d/%d", n_expressing, n_active))
    }

    # Solve diffusion
    domain_area <- pi * config$domain$max_radius^2
    n_density <- n_receivers / domain_area

    diffusion_result <- solve_diffusion(
      all_positions[sender_indices, , drop = FALSE],
      factor_expr[sender_indices],
      all_positions,
      config$diffusion,
      n_density,
      position_dict,
      config$stochastic$p_receiver_respond_max
    )
    concentrations <- diffusion_result$concentrations
    lambda_val <- diffusion_result$lambda

    # Generate response
    response_result <- generate_response_expression(
      nrow(all_positions), receiver_indices, concentrations,
      config$expression, config$stochastic
    )
    response_expr <- response_result$response_expr
    responding_mask <- response_result$responding_mask

    n_responding <- sum(responding_mask)
    if (verbose) {
      message(sprintf("  Responding receivers: %d/%d", n_responding, n_receivers))
    }

    # Compute I_ND
    ind_results <- compute_ind_radial(
      sender_indices, receiver_indices, all_positions,
      factor_expr, response_expr, config$analysis, verbose
    )

    results[[as.character(frac)]] <- list(
      lambda = lambda_val,
      ind_curve = ind_results,
      n_expressing = n_expressing,
      n_responding = n_responding,
      assignments = assignments,
      sender_indices = sender_indices,
      receiver_indices = receiver_indices,
      factor_expr = factor_expr,
      response_expr = response_expr,
      positions = all_positions,
      concentrations = concentrations
    )

    if (verbose) {
      message(sprintf("  Lambda: %.0f um", lambda_val))
    }
  }

  # Save results
  .save_simulation_results(results, config)

  if (verbose) {
    message(sprintf("\nSimulation complete. Results in %s", config$output_dir))
  }

  results
}

#' Save simulation results to files
#' @keywords internal
.save_simulation_results <- function(results, config) {
  output_dir <- config$output_dir

  # I_ND curves
  all_rows <- list()
  for (frac_str in names(results)) {
    data <- results[[frac_str]]
    ind_df <- data$ind_curve
    ind_df$receiver_fraction <- as.numeric(frac_str)
    ind_df$lambda <- data$lambda
    all_rows[[frac_str]] <- ind_df
  }

  if (length(all_rows) > 0) {
    combined_df <- do.call(rbind, all_rows)
    write.csv(combined_df, file.path(output_dir, "ind_results.csv"), row.names = FALSE)
  }

  # Summary
  summary_df <- data.frame(
    receiver_fraction = as.numeric(names(results)),
    lambda = sapply(results, function(x) x$lambda),
    n_expressing = sapply(results, function(x) x$n_expressing),
    n_responding = sapply(results, function(x) x$n_responding)
  )
  write.csv(summary_df, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)

  # Full results as RDS
  saveRDS(results, file.path(output_dir, "simulation_results.rds"))
}

# =============================================================================
# VISUALIZATION
# =============================================================================

#' Plot I_ND curves
#'
#' @param results Simulation results from run_simulation()
#' @param output_path Output file path (PNG)
#' @param smooth Apply Gaussian smoothing
#' @export
plot_ind_curves <- function(results, output_path = NULL, smooth = TRUE) {
  fractions <- as.numeric(names(results))
  colors <- grDevices::colorRampPalette(c("blue", "green", "orange", "red"))(length(fractions))

  if (!is.null(output_path)) {
    png(output_path, width = 10, height = 8, units = "in", res = 300)
  }

  plot(NULL, xlim = c(0, 5000), ylim = c(-1.1, 1.1),
       xlab = "Distance from Senders (um)", ylab = "I_ND",
       main = "I_ND vs Distance", cex.lab = 1.2, cex.main = 1.4)
  grid()
  abline(h = 0, col = "gray50", lty = 2)

  for (i in seq_along(fractions)) {
    frac <- fractions[i]
    frac_str <- as.character(frac)
    data <- results[[frac_str]]

    radii <- data$ind_curve$radius
    vals <- data$ind_curve$I_ND
    lambda_v <- data$lambda

    if (smooth && length(vals) > 3) {
      vals <- stats::filter(vals, rep(1/3, 3), sides = 2)
      vals[is.na(vals)] <- data$ind_curve$I_ND[is.na(vals)]
    }

    lines(radii, vals, col = colors[i], lwd = 2)
    abline(v = lambda_v, col = colors[i], lty = 3, lwd = 1.5)
  }

  legend("topright",
         legend = sprintf("%.0f%% Receivers (lambda=%.0f)",
                          fractions * 100, sapply(results, function(x) x$lambda)),
         col = colors, lwd = 2, cex = 0.9)

  if (!is.null(output_path)) {
    dev.off()
    message(sprintf("Saved plot to %s", output_path))
  }
}

#' Plot cell distribution
#'
#' @param results Simulation results (single receiver fraction)
#' @param output_path Output file path (PNG)
#' @export
plot_cell_distribution <- function(results, output_path = NULL) {
  positions <- results$positions
  sender_idx <- results$sender_indices
  receiver_idx <- results$receiver_indices

  all_idx <- seq_len(nrow(positions))
  other_idx <- setdiff(all_idx, c(sender_idx, receiver_idx))

  if (!is.null(output_path)) {
    png(output_path, width = 10, height = 10, units = "in", res = 300)
  }

  plot(positions[other_idx, 1], positions[other_idx, 2],
       pch = 16, cex = 0.1, col = "lightgray",
       xlab = "X (um)", ylab = "Y (um)", main = "Cell Distribution",
       asp = 1)

  points(positions[receiver_idx, 1], positions[receiver_idx, 2],
         pch = 16, cex = 0.2, col = rgb(0, 0, 1, 0.5))

  points(positions[sender_idx, 1], positions[sender_idx, 2],
         pch = 17, cex = 1, col = "red")

  legend("topright",
         legend = c("Other", "Receivers", "Senders"),
         col = c("lightgray", "blue", "red"),
         pch = c(16, 16, 17),
         cex = 0.9)

  if (!is.null(output_path)) {
    dev.off()
    message(sprintf("Saved plot to %s", output_path))
  }
}
