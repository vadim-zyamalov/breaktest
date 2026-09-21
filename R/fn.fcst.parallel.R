# Scoped doSNOW execution, shared by comparison and Monte Carlo simulation.
# The sequential path needs no parallel backend and uses the same task code.

.fcst.n.cores <- function(n.cores) {
  if (!is.numeric(n.cores) || length(n.cores) != 1L || !is.finite(n.cores) ||
      n.cores < 1 || n.cores != floor(n.cores) || n.cores > .Machine$integer.max) {
    stop("ERROR! fcst.parallel: n.cores must be a positive integer")
  }
  as.integer(n.cores)
}


.fcst.seed <- function(seed) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed != floor(seed) || abs(seed) > .Machine$integer.max)) {
    stop("ERROR! fcst.parallel: seed must be NULL or a finite integer")
  }
  seed
}


.fcst.rng.state <- function() {
  list(kind = RNGkind(),
       seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
         get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)
}


.fcst.rng.restore <- function(state) {
  if (!identical(RNGkind(), state$kind)) do.call(RNGkind, as.list(state$kind))
  if (is.null(state$seed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  }
}


# Assign streams to tasks, not to workers: load balancing and core count
# must not change a stochastic forecast callback's random numbers.
.fcst.rng.streams <- function(n, seed) {
  if (is.null(seed)) return(NULL)
  state <- .fcst.rng.state()
  on.exit(.fcst.rng.restore(state), add = TRUE)
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
  set.seed(.fcst.seed(seed))
  streams <- vector("list", n)
  stream <- get(".Random.seed", .GlobalEnv)
  for (i in seq_len(n)) {
    streams[[i]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  streams
}


.fcst.parallel.task <- function(job, FUN, context) {
  if (!is.null(job$seed)) assign(".Random.seed", job$seed, envir = .GlobalEnv)
  FUN(job$task, context)
}


# Loading namespaces / allocating a socket port must not consume the stream
# used to generate Monte Carlo innovations in the calling process.
.fcst.parallel.start <- function(n.cores) {
  state <- .fcst.rng.state()
  on.exit(.fcst.rng.restore(state), add = TRUE)
  for (package in c("foreach", "doSNOW", "iterators")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("ERROR! fcst.parallel: install ", package, " to use n.cores > 1")
    }
  }
  parallel::makePSOCKcluster(n.cores)
}


# foreach has a public setter but no getter for the complete registration.
# Keep the only access to its private registration store in these two helpers.
# Restoring only registerDoSEQ() would destroy a caller's existing backend.
.fcst.backend.state <- function() {
  env <- get(".foreachGlobals", asNamespace("foreach"))
  fields <- intersect(c("fun", "data", "info"), ls(env, all.names = TRUE))
  mget(fields, env, inherits = FALSE)
}


.fcst.backend.restore <- function(state) {
  env <- get(".foreachGlobals", asNamespace("foreach"))
  fields <- intersect(c("fun", "data", "info"), ls(env, all.names = TRUE))
  if (length(fields)) rm(list = fields, envir = env)
  list2env(state, envir = env)
  invisible(NULL)
}


# Installed-package functions resolve through their namespace. The standalone
# loader records exactly which top-level bindings came from its source files,
# so workers need neither a separately installed breaktest nor the user's
# entire global workspace.
.fcst.worker.setup <- function(cluster, export, packages) {
  env <- environment(.fcst.worker.setup)
  if (isNamespace(env)) {
    packages <- union(getNamespaceName(env), packages)
  } else {
    if (!exists(".fcst.source.names", env, inherits = FALSE)) {
      stop("ERROR! fcst.parallel: use bt.load() or install the module in breaktest")
    }
    bindings <- get(".fcst.source.names", env, inherits = FALSE)
    if (length(intersect(names(export), bindings))) {
      stop("ERROR! fcst.parallel: parallel.export must not overwrite module bindings")
    }
    parallel::clusterExport(cluster, bindings, envir = env)
    packages <- union("Rfast", packages)
  }
  if (length(export)) {
    parallel::clusterExport(cluster, names(export),
                             envir = list2env(export, parent = emptyenv()))
  }
  packages
}


.fcst.parallel.options <- function(export, packages) {
  if (!is.list(export) || (length(export) &&
      (is.null(names(export)) || anyNA(names(export)) ||
       any(names(export) == "") || anyDuplicated(names(export))))) {
    stop("ERROR! fcst.parallel: parallel.export must be a distinctly named list")
  }
  if (!is.character(packages) || anyNA(packages) || any(packages == "")) {
    stop("ERROR! fcst.parallel: parallel.packages must contain package names")
  }
}


# make.task is evaluated lazily on the master. In the simulator it supplies
# bounded blocks of random innovations in the legacy RNG order.
# receive, when supplied, assembles ordered results without retaining a second
# full copy of the simulation arrays. Errors in either tasks or the collector
# propagate after cluster cleanup.
.fcst.parallel.map <- function(
  n, FUN, context = NULL, n.cores = 1L, make.task = identity,
  receive = NULL, seed = NULL, verbose = FALSE, label = "fcst.parallel",
  export = list(), packages = character()
) {
  n.cores <- min(.fcst.n.cores(n.cores), n)
  .fcst.parallel.options(export, packages)
  streams <- .fcst.rng.streams(n, .fcst.seed(seed))
  if (!is.null(streams)) {
    state <- .fcst.rng.state()
    on.exit(.fcst.rng.restore(state), add = TRUE)
  }
  progress <- function(done) {
    if (verbose && (done %% max(1L, round(n / 10)) == 0L || done == n)) {
      message(label, ": ", done, "/", n)
    }
  }
  job <- function(i) list(task = make.task(i), seed = if (!is.null(streams)) streams[[i]])
  if (n.cores == 1L) {
    result <- if (is.null(receive)) vector("list", n) else NULL
    for (i in seq_len(n)) {
      value <- .fcst.parallel.task(job(i), FUN, context)
      if (is.null(receive)) result[i] <- list(value) else receive(value)
      progress(i)
    }
    return(result)
  }

  cluster <- .fcst.parallel.start(n.cores)
  backend <- NULL
  on.exit({
    try(parallel::stopCluster(cluster), silent = TRUE)
    if (!is.null(backend)) .fcst.backend.restore(backend)
  }, add = TRUE)
  backend <- .fcst.backend.state()
  packages <- .fcst.worker.setup(cluster, export, packages)
  doSNOW::registerDoSNOW(cluster)
  i <- 0L
  tasks <- iterators::iter(function() {
    i <<- i + 1L
    if (i > n) stop("StopIteration")
    job(i)
  })
  collection.error <- NULL
  combine <- if (is.null(receive)) function(acc, value) c(acc, list(value)) else {
    function(acc, value) {
      if (!inherits(value, "error") && is.null(collection.error)) {
        tryCatch(receive(value), error = function(e) collection.error <<- e)
      }
      NULL
    }
  }
  task <- NULL  # foreach iteration variable, for static code checking.
  result <- foreach::"%dopar%"(
    foreach::foreach(task = tasks, .combine = combine, .init = list(),
                     .inorder = TRUE, .multicombine = FALSE, .errorhandling = "stop",
                     .packages = packages, .options.snow = list(progress = progress)),
    .fcst.parallel.task(task, FUN, context)
  )
  if (!is.null(collection.error)) stop(collection.error)
  result
}
