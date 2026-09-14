#' Maximum-entropy ecological network model
#'
#' Constructs a directed network inference backend that generates binary
#' adjacency matrices under hard structural constraints and returns their
#' marginal species-pair probabilities. Rows are consumers, columns are
#' resources, and `matrix[i, j] = 1` means consumer `i` uses resource `j`.
#'
#' For `constraint = "connectance"` or `"links"`, networks are sampled
#' uniformly from all allowed adjacency matrices with exactly the requested
#' number of links. This is the maximum-entropy distribution under a hard link
#' count. For `constraint = "degree_sequence"`, a feasible matrix is found by
#' maximum flow and randomized with degree-preserving checkerboard swaps. This
#' is an approximate uniform MCMC ensemble conditional on the supplied margins.
#'
#' The implementation follows the constrained, minimally biased ensemble
#' rationale of Banville, Gravel & Poisot (2023), but does not yet implement
#' their SVD-entropy objective or simulated-annealing optimizer.
#'
#' @param row_group Assemblage group represented by matrix rows (consumers).
#' @param column_group Assemblage group represented by matrix columns
#'   (resources). May equal `row_group` for a directed unipartite food web.
#' @param constraint One of `"connectance"`, `"links"`, or
#'   `"degree_sequence"`.
#' @param connectance Desired connectance in `[0, 1]`. Converted to an exact
#'   local link count with `round(connectance * allowed_dyads)`.
#' @param links Exact number of interactions in every local ensemble member.
#' @param out_degree Row sums (numbers of resources per consumer).
#' @param in_degree Column sums (numbers of consumers per resource).
#'   Degree constraints may be named numeric vectors or functions accepting the
#'   locally present species names and returning a numeric vector.
#' @param ensemble_size Number of adjacency matrices used for marginal
#'   probabilities.
#' @param self_links Whether a species may interact with itself when its name
#'   appears in both groups.
#' @param seed Integer seed. A cell identifier is added to this seed so spatial
#'   cells receive distinct reproducible ensembles. The previous global random
#'   state is restored after inference. Use `NULL` to use normal R RNG behavior.
#' @param burn_in Number of degree-preserving swap proposals before recording a
#'   degree-constrained ensemble.
#' @param thin Number of swap proposals between recorded matrices.
#' @param keep_ensemble Store binary ensemble members in network metadata.
#' @return An `inference_model` compatible with [infer_network()] and
#'   [infer_networks()].
#' @references Banville F, Gravel D, Poisot T (2023). What constrains food
#'   webs? A maximum entropy framework for predicting their structure with
#'   minimal biases. PLOS Computational Biology 19: e1011458.
#'   \doi{10.1371/journal.pcbi.1011458}
#' @examples
#' assemblage <- new_assemblage(list(food_web = c("a", "b", "c")))
#' model <- maxent_model(
#'   row_group = "food_web",
#'   constraint = "connectance",
#'   connectance = 0.25,
#'   ensemble_size = 50,
#'   self_links = FALSE,
#'   seed = 42
#' )
#' network <- infer_network(assemblage, model)
#' network$matrix
#' @export
maxent_model <- function(row_group, column_group = row_group,
                         constraint = c("connectance", "links", "degree_sequence"),
                         connectance = NULL, links = NULL,
                         out_degree = NULL, in_degree = NULL,
                         ensemble_size = 100L, self_links = FALSE,
                         seed = 1L, burn_in = 1000L, thin = 100L,
                         keep_ensemble = FALSE) {
  constraint <- match.arg(constraint)
  .maxent_scalar_name(row_group, "row_group")
  .maxent_scalar_name(column_group, "column_group")
  .maxent_count(ensemble_size, "ensemble_size", minimum = 1L)
  .maxent_count(burn_in, "burn_in", minimum = 0L)
  .maxent_count(thin, "thin", minimum = 0L)
  if (!is.logical(self_links) || length(self_links) != 1L || is.na(self_links)) {
    stop("`self_links` must be TRUE or FALSE.")
  }
  if (!is.logical(keep_ensemble) || length(keep_ensemble) != 1L || is.na(keep_ensemble)) {
    stop("`keep_ensemble` must be TRUE or FALSE.")
  }
  if (!is.null(seed)) {
    .maxent_count(seed, "seed", minimum = 0L)
    if (seed > .Machine$integer.max) stop("`seed` must not exceed `.Machine$integer.max`.")
  }

  if (constraint == "connectance") {
    if (!is.numeric(connectance) || length(connectance) != 1L || !is.finite(connectance) ||
        connectance < 0 || connectance > 1) {
      stop("`connectance` must be one finite value in [0, 1].")
    }
    if (!is.null(links) || !is.null(out_degree) || !is.null(in_degree)) {
      stop("Connectance constraints cannot be combined with links or degree sequences.")
    }
  } else if (constraint == "links") {
    .maxent_count(links, "links", minimum = 0L)
    if (!is.null(connectance) || !is.null(out_degree) || !is.null(in_degree)) {
      stop("Link-count constraints cannot be combined with connectance or degree sequences.")
    }
  } else {
    if (is.null(out_degree) || is.null(in_degree)) {
      stop("Degree-sequence constraints require both `out_degree` and `in_degree`.")
    }
    if (!is.null(connectance) || !is.null(links)) {
      stop("Degree-sequence constraints cannot be combined with connectance or links.")
    }
    .maxent_degree_input(out_degree, "out_degree")
    .maxent_degree_input(in_degree, "in_degree")
  }

  predict_maxent <- function(model, assemblage) {
    p <- model$parameters
    required <- unique(c(p$row_group, p$column_group))
    if (!all(required %in% names(assemblage$species))) {
      stop("Assemblage lacks model group(s): ", paste(setdiff(required, names(assemblage$species)), collapse = ", "), ".")
    }
    rows <- unique(assemblage$species[[p$row_group]])
    columns <- unique(assemblage$species[[p$column_group]])
    if (!length(rows) || !length(columns)) stop("Maximum-entropy inference requires non-empty local groups.")

    allowed <- matrix(TRUE, length(rows), length(columns), dimnames = list(rows, columns))
    if (!p$self_links) allowed[outer(rows, columns, `==`)] <- FALSE
    allowed_dyads <- sum(allowed)

    local_seed <- .maxent_cell_seed(p$seed, assemblage$cell_id)
    sampled <- .maxent_with_seed(local_seed, function() {
      if (p$constraint %in% c("connectance", "links")) {
        link_count <- if (p$constraint == "connectance") {
          as.integer(round(p$connectance * allowed_dyads))
        } else {
          as.integer(p$links)
        }
        if (link_count > allowed_dyads) {
          stop("Requested ", link_count, " links, but only ", allowed_dyads,
               " dyads are allowed in this assemblage.")
        }
        .maxent_sample_links(allowed, link_count, p$ensemble_size, p$keep_ensemble)
      } else {
        out_degree <- .maxent_resolve_degree(p$out_degree, rows, "out_degree")
        in_degree <- .maxent_resolve_degree(p$in_degree, columns, "in_degree")
        .maxent_sample_degrees(allowed, out_degree, in_degree,
                               p$ensemble_size, p$burn_in, p$thin, p$keep_ensemble)
      }
    })
    if (p$constraint == "connectance") {
      sampled$constraints$connectance <- p$connectance
    }

    metadata <- list(
      cell_id = assemblage$cell_id,
      model = "maximum entropy",
      constraint = p$constraint,
      ensemble_size = p$ensemble_size,
      seed = local_seed,
      self_links = p$self_links,
      sampler = sampled$sampler,
      allowed_dyads = allowed_dyads,
      constraints = sampled$constraints
    )
    if (p$keep_ensemble) metadata$ensemble <- sampled$ensemble

    new_ecological_network(
      sampled$probability,
      type = "probability",
      metadata = metadata,
      row_group = p$row_group,
      column_group = p$column_group,
      directed = TRUE,
      self_links = p$self_links
    )
  }

  new_inference_model(
    predict_maxent,
    name = "maximum entropy",
    parameters = list(
      row_group = row_group,
      column_group = column_group,
      constraint = constraint,
      connectance = connectance,
      links = links,
      out_degree = out_degree,
      in_degree = in_degree,
      ensemble_size = as.integer(ensemble_size),
      self_links = self_links,
      seed = if (is.null(seed)) NULL else as.integer(seed),
      burn_in = as.integer(burn_in),
      thin = as.integer(thin),
      keep_ensemble = keep_ensemble
    ),
    subclass = "maxent_model"
  )
}

.maxent_scalar_name <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", label, "` must be one non-empty group name.")
  }
  invisible(TRUE)
}

.maxent_count <- function(x, label, minimum = 0L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum || x != floor(x)) {
    stop("`", label, "` must be one integer >= ", minimum, ".")
  }
  invisible(TRUE)
}

.maxent_degree_input <- function(x, label) {
  if (!is.function(x) && !is.numeric(x)) stop("`", label, "` must be numeric or a function.")
  invisible(TRUE)
}

.maxent_resolve_degree <- function(x, species, label) {
  values <- if (is.function(x)) x(species) else x
  if (!is.numeric(values) || any(!is.finite(values))) stop("`", label, "` must return finite numeric values.")
  if (!is.null(names(values))) {
    missing <- setdiff(species, names(values))
    if (length(missing)) stop("`", label, "` is missing species: ", paste(missing, collapse = ", "), ".")
    values <- values[species]
  } else if (length(values) != length(species)) {
    stop("Unnamed `", label, "` must have one value per locally present species.")
  }
  if (any(values < 0 | values != floor(values))) stop("`", label, "` must contain non-negative integers.")
  setNames(as.integer(values), species)
}

.maxent_cell_seed <- function(seed, cell_id) {
  if (is.null(seed)) return(NULL)
  offset <- if (length(cell_id) && !is.na(cell_id)) as.double(cell_id) else 0
  value <- (as.double(seed) + offset) %% .Machine$integer.max
  if (value == 0) value <- 1
  as.integer(value)
}

.maxent_with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (existed) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fun()
}

.maxent_sample_links <- function(allowed, links, ensemble_size, keep) {
  candidates <- which(allowed)
  accumulator <- matrix(0, nrow(allowed), ncol(allowed), dimnames = dimnames(allowed))
  ensemble <- if (keep) vector("list", ensemble_size) else NULL
  for (r in seq_len(ensemble_size)) {
    adjacency <- matrix(0, nrow(allowed), ncol(allowed), dimnames = dimnames(allowed))
    if (links > 0L) {
      selected <- candidates[sample.int(length(candidates), links, replace = FALSE)]
      adjacency[selected] <- 1
    }
    accumulator <- accumulator + adjacency
    if (keep) ensemble[[r]] <- adjacency
  }
  list(probability = accumulator / ensemble_size, ensemble = ensemble,
       sampler = "uniform fixed-link ensemble",
       constraints = list(links = links))
}

.maxent_sample_degrees <- function(allowed, out_degree, in_degree,
                                   ensemble_size, burn_in, thin, keep) {
  if (sum(out_degree) != sum(in_degree)) {
    stop("Degree sequences are inconsistent: sums of `out_degree` and `in_degree` differ.")
  }
  row_capacity <- rowSums(allowed)
  column_capacity <- colSums(allowed)
  if (any(out_degree > row_capacity)) {
    stop("Infeasible `out_degree`: at least one consumer requests more allowed resources than exist.")
  }
  if (any(in_degree > column_capacity)) {
    stop("Infeasible `in_degree`: at least one resource requests more allowed consumers than exist.")
  }

  adjacency <- .maxent_feasible_matrix(allowed, out_degree, in_degree)
  if (burn_in > 0L) adjacency <- .maxent_swap_chain(adjacency, allowed, burn_in)
  accumulator <- matrix(0, nrow(allowed), ncol(allowed), dimnames = dimnames(allowed))
  ensemble <- if (keep) vector("list", ensemble_size) else NULL
  for (r in seq_len(ensemble_size)) {
    if (thin > 0L) adjacency <- .maxent_swap_chain(adjacency, allowed, thin)
    accumulator <- accumulator + adjacency
    if (keep) ensemble[[r]] <- adjacency
  }
  list(probability = accumulator / ensemble_size, ensemble = ensemble,
       sampler = "degree-preserving swap MCMC",
       constraints = list(out_degree = out_degree, in_degree = in_degree))
}

.maxent_feasible_matrix <- function(allowed, out_degree, in_degree) {
  nr <- nrow(allowed); nc <- ncol(allowed)
  source <- 1L
  row_nodes <- seq_len(nr) + 1L
  column_nodes <- seq_len(nc) + 1L + nr
  sink <- nr + nc + 2L
  capacity <- matrix(0, sink, sink)
  capacity[source, row_nodes] <- out_degree
  capacity[row_nodes, column_nodes] <- allowed * 1
  capacity[column_nodes, sink] <- in_degree
  residual <- capacity
  flow <- 0L
  target <- sum(out_degree)

  while (flow < target) {
    parent <- rep(NA_integer_, sink)
    parent[source] <- 0L
    queue <- source
    head <- 1L
    while (head <= length(queue) && is.na(parent[sink])) {
      u <- queue[head]; head <- head + 1L
      next_nodes <- which(residual[u, ] > 0 & is.na(parent))
      if (length(next_nodes)) {
        parent[next_nodes] <- u
        queue <- c(queue, next_nodes)
      }
    }
    if (is.na(parent[sink])) {
      stop("Degree sequences are infeasible for the locally allowed dyads.")
    }
    increment <- Inf
    v <- sink
    while (v != source) {
      u <- parent[v]
      increment <- min(increment, residual[u, v])
      v <- u
    }
    v <- sink
    while (v != source) {
      u <- parent[v]
      residual[u, v] <- residual[u, v] - increment
      residual[v, u] <- residual[v, u] + increment
      v <- u
    }
    flow <- flow + increment
  }

  adjacency <- capacity[row_nodes, column_nodes, drop = FALSE] -
    residual[row_nodes, column_nodes, drop = FALSE]
  storage.mode(adjacency) <- "numeric"
  dimnames(adjacency) <- dimnames(allowed)
  adjacency
}

.maxent_swap_chain <- function(adjacency, allowed, steps) {
  if (nrow(adjacency) < 2L || ncol(adjacency) < 2L) return(adjacency)
  for (step in seq_len(steps)) {
    rows <- sample.int(nrow(adjacency), 2L, replace = FALSE)
    left <- which(adjacency[rows[1], ] == 1 & adjacency[rows[2], ] == 0 & allowed[rows[2], ])
    right <- which(adjacency[rows[2], ] == 1 & adjacency[rows[1], ] == 0 & allowed[rows[1], ])
    if (!length(left) || !length(right)) next
    c1 <- left[sample.int(length(left), 1L)]
    c2 <- right[sample.int(length(right), 1L)]
    if (c1 == c2) next
    adjacency[rows[1], c1] <- 0
    adjacency[rows[2], c2] <- 0
    adjacency[rows[1], c2] <- 1
    adjacency[rows[2], c1] <- 1
  }
  adjacency
}
