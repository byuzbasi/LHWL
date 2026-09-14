
###############################################################################
# LHWL package review-ready audit
#
# This script records the local checks used before submitting a revision,
# creating a review supplement, or preparing a package release. It does not
# install or update dependencies and refuses to overwrite an audit directory.
###############################################################################

lhwl_utc_time <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

lhwl_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(unname(digest::digest(path, algo = "sha256", file = TRUE)))
  }

  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    value <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
    if (is.null(attr(value, "status")) || identical(attr(value, "status"), 0L)) {
      return(strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]])
    }
  }

  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    value <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
    if (is.null(attr(value, "status")) || identical(attr(value, "status"), 0L)) {
      return(strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]])
    }
  }

  stop("SHA-256 requires the digest package, sha256sum, or shasum.", call. = FALSE)
}

lhwl_write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    na = ""
  )
}

lhwl_source_manifest <- function(pkg_dir) {
  files <- list.files(
    pkg_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  info <- file.info(files)
  files <- files[!is.na(info$isdir) & !info$isdir]
  relative <- substring(files, nchar(pkg_dir) + 2L)
  keep <- grepl(
    "^(DESCRIPTION|NAMESPACE|NEWS[.]md|REVIEW_READY_CHECKLIST[.]md|[.]Rbuildignore|R/|src/|tests/|man/|vignettes/|tools/|data/|inst/)",
    relative
  )
  generated <- grepl("[.](o|so|dll|dylib)$", relative, ignore.case = TRUE)
  files <- files[keep & !generated]
  relative <- relative[keep & !generated]
  ord <- order(relative)
  files <- files[ord]
  relative <- relative[ord]
  info <- file.info(files)

  data.frame(
    path = relative,
    bytes = unname(info$size),
    sha256 = vapply(files, lhwl_sha256, character(1)),
    stringsAsFactors = FALSE
  )
}

lhwl_r_cmd_config <- function(key) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "config", key),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  data.frame(
    key = key,
    status = as.integer(status),
    value = paste(output, collapse = " "),
    stringsAsFactors = FALSE
  )
}

lhwl_library_dirs <- function(toolchain) {
  flag_rows <- toolchain$key %in% c(
    "LAPACK_LIBS", "BLAS_LIBS", "FLIBS", "SHLIB_OPENMP_CXXFLAGS"
  ) & toolchain$status == 0L
  values <- toolchain$value[flag_rows]
  keys <- toolchain$key[flag_rows]
  out <- list()

  for (i in seq_along(values)) {
    tokens <- tryCatch(
      scan(text = values[[i]], what = character(), quiet = TRUE),
      error = function(e) character()
    )
    j <- 1L
    while (j <= length(tokens)) {
      path <- NULL
      if (identical(tokens[[j]], "-L") && j < length(tokens)) {
        path <- tokens[[j + 1L]]
        j <- j + 1L
      } else if (startsWith(tokens[[j]], "-L") && nchar(tokens[[j]]) > 2L) {
        path <- substring(tokens[[j]], 3L)
      }
      if (!is.null(path)) {
        path <- path.expand(path)
        out[[length(out) + 1L]] <- data.frame(
          source = keys[[i]],
          path = path,
          exists = dir.exists(path),
          stringsAsFactors = FALSE
        )
      }
      j <- j + 1L
    }
  }

  if (!length(out)) {
    return(data.frame(
      source = character(), path = character(), exists = logical()
    ))
  }
  unique(do.call(rbind, out))
}

lhwl_capture_step <- function(label, expression, log_dir) {
  started <- Sys.time()
  transcript <- character()
  error_message <- NULL
  value <- withCallingHandlers(
    tryCatch(
      {
        printed <- capture.output(result <- force(expression), type = "output")
        transcript <- c(transcript, printed)
        result
      },
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    ),
    message = function(m) {
      transcript <<- c(transcript, paste0("MESSAGE: ", conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      transcript <<- c(transcript, paste0("WARNING: ", conditionMessage(w)))
    }
  )
  ended <- Sys.time()
  if (!is.null(error_message)) {
    transcript <- c(transcript, paste0("ERROR: ", error_message))
  }
  writeLines(transcript, file.path(log_dir, paste0(label, ".log")), useBytes = TRUE)

  list(
    value = value,
    record = data.frame(
      step = label,
      status = if (is.null(error_message)) "PASSED" else "FAILED",
      started_utc = format(started, tz = "UTC", usetz = TRUE),
      ended_utc = format(ended, tz = "UTC", usetz = TRUE),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      detail = if (is.null(error_message)) "" else error_message,
      stringsAsFactors = FALSE
    )
  )
}

run_lhwr_review_audit <- function(pkg_dir = "LHWL_package",
                                  run_check = FALSE,
                                  run_examples = FALSE,
                                  document = TRUE,
                                  output_dir = NULL,
                                  check_args = character()) {
  if (!dir.exists(pkg_dir)) {
    stop("Package directory not found: ", pkg_dir, call. = FALSE)
  }
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("The audit requires the already-installed devtools package.", call. = FALSE)
  }

  pkg_dir <- normalizePath(pkg_dir, winslash = "/", mustWork = TRUE)
  description <- read.dcf(file.path(pkg_dir, "DESCRIPTION"))
  package_version <- unname(description[1L, "Version"])
  if (is.null(output_dir)) {
    stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
    output_dir <- file.path(
      dirname(pkg_dir),
      "validation",
      paste0("LHWL_", package_version, "_review_audit_", stamp)
    )
  }
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("Refusing to overwrite audit output: ", output_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  log_dir <- file.path(output_dir, "logs")
  dir.create(log_dir, showWarnings = FALSE)

  started <- Sys.time()
  audit_status <- "RUNNING"
  failure <- ""
  steps <- data.frame(
    step = character(), status = character(), started_utc = character(),
    ended_utc = character(), elapsed_seconds = numeric(), detail = character(),
    stringsAsFactors = FALSE
  )

  finalize <- function() {
    ended <- Sys.time()
    if (identical(audit_status, "RUNNING")) {
      audit_status <<- "PASSED"
    }
    lhwl_write_tsv(steps, file.path(output_dir, "audit_steps.tsv"))
    lhwl_write_tsv(
      data.frame(
        package = unname(description[1L, "Package"]),
        version = package_version,
        package_directory = pkg_dir,
        output_directory = output_dir,
        status = audit_status,
        failure = failure,
        started_utc = format(started, tz = "UTC", usetz = TRUE),
        ended_utc = format(ended, tz = "UTC", usetz = TRUE),
        elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
        stringsAsFactors = FALSE
      ),
      file.path(output_dir, "audit_summary.tsv")
    )
    lhwl_write_tsv(
      lhwl_source_manifest(pkg_dir),
      file.path(output_dir, "source_manifest_after.tsv")
    )
    artifacts <- list.files(
      output_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE
    )
    artifacts <- artifacts[!file.info(artifacts)$isdir]
    artifacts <- artifacts[basename(artifacts) != "SHA256SUMS"]
    artifact_relative <- substring(artifacts, nchar(output_dir) + 2L)
    ord <- order(artifact_relative)
    lines <- paste0(
      vapply(artifacts[ord], lhwl_sha256, character(1)),
      "  ",
      artifact_relative[ord]
    )
    writeLines(lines, file.path(output_dir, "SHA256SUMS"), useBytes = TRUE)
  }
  on.exit(finalize(), add = TRUE)

  lhwl_write_tsv(
    lhwl_source_manifest(pkg_dir),
    file.path(output_dir, "source_manifest_before.tsv")
  )

  makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = "")
  if (!nzchar(makevars_user)) {
    makevars_user <- path.expand("~/.R/Makevars")
  } else {
    makevars_user <- path.expand(makevars_user)
  }
  makevars_exists <- file.exists(makevars_user)
  makevars_sha256 <- if (makevars_exists) lhwl_sha256(makevars_user) else NA_character_
  makevars_lines <- if (makevars_exists) {
    readLines(makevars_user, warn = FALSE)
  } else {
    character()
  }
  relevant_makevars <- grep(
    "^[[:space:]]*(FC|F77|FLIBS|SHLIB_OPENMP_CFLAGS|SHLIB_OPENMP_CXXFLAGS|CPPFLAGS|LDFLAGS)[[:space:]]*[+:?]?=",
    makevars_lines
  )
  lhwl_write_tsv(
    data.frame(
      path = makevars_user,
      exists = makevars_exists,
      sha256 = makevars_sha256,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "makevars_user.tsv")
  )
  lhwl_write_tsv(
    data.frame(
      line_number = relevant_makevars,
      assignment = makevars_lines[relevant_makevars],
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "makevars_user_assignments.tsv")
  )

  config_keys <- c(
    "CC", "CXX", "CXX17", "FC", "F77", "CPPFLAGS", "CFLAGS",
    "CXXFLAGS", "CXX17FLAGS", "SHLIB_OPENMP_CFLAGS",
    "SHLIB_OPENMP_CXXFLAGS", "LAPACK_LIBS", "BLAS_LIBS", "FLIBS"
  )
  toolchain <- do.call(rbind, lapply(config_keys, lhwl_r_cmd_config))
  lhwl_write_tsv(toolchain, file.path(output_dir, "r_cmd_config.tsv"))
  library_dirs <- lhwl_library_dirs(toolchain)
  lhwl_write_tsv(library_dirs, file.path(output_dir, "link_library_dirs.tsv"))

  compiler_keys <- c("CC", "CXX17", "FC")
  compiler_config <- toolchain[match(compiler_keys, toolchain$key), , drop = FALSE]
  compiler_command <- vapply(
    compiler_config$value,
    function(value) {
      tokens <- tryCatch(
        scan(text = value, what = character(), quiet = TRUE),
        error = function(e) character()
      )
      if (length(tokens)) tokens[[1L]] else ""
    },
    character(1)
  )
  compiler_path <- unname(Sys.which(compiler_command))
  compiler_path[
    !nzchar(compiler_path) & file.exists(compiler_command)
  ] <- compiler_command[
    !nzchar(compiler_path) & file.exists(compiler_command)
  ]
  compiler_inventory <- data.frame(
    key = compiler_keys,
    command = compiler_command,
    resolved_path = compiler_path,
    exists = compiler_config$status == 0L & nzchar(compiler_path),
    stringsAsFactors = FALSE
  )
  lhwl_write_tsv(
    compiler_inventory,
    file.path(output_dir, "compiler_executables.tsv")
  )

  environment <- data.frame(
    variable = c(
      "R_VERSION", "R_PLATFORM", "R_HOME", "R_MAKEVARS_USER",
      "R_MAKEVARS_USER_SHA256", ".libPaths", "OMP_NUM_THREADS",
      "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS"
    ),
    value = c(
      R.version.string, R.version$platform, R.home(), makevars_user,
      makevars_sha256, paste(.libPaths(), collapse = .Platform$path.sep),
      Sys.getenv("OMP_NUM_THREADS"), Sys.getenv("OPENBLAS_NUM_THREADS"),
      Sys.getenv("MKL_NUM_THREADS"), Sys.getenv("BLIS_NUM_THREADS"),
      Sys.getenv("VECLIB_MAXIMUM_THREADS")
    ),
    exists = c(
      TRUE, TRUE, dir.exists(R.home()), makevars_exists, makevars_exists,
      all(dir.exists(.libPaths())), rep(NA, 5L)
    ),
    stringsAsFactors = FALSE
  )
  lhwl_write_tsv(environment, file.path(output_dir, "runtime_environment.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

  packages <- c(
    "Rcpp", "RcppArmadillo", "MASS", "testthat", "devtools", "roxygen2",
    "pkgbuild", "rcmdcheck", "digest", "knitr", "rmarkdown"
  )
  package_versions <- data.frame(
    package = packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(
      packages,
      function(x) if (requireNamespace(x, quietly = TRUE)) {
        as.character(utils::packageVersion(x))
      } else {
        NA_character_
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  lhwl_write_tsv(package_versions, file.path(output_dir, "package_versions.tsv"))

  required_config <- c("CC", "CXX17", "FC", "LAPACK_LIBS", "BLAS_LIBS", "FLIBS")
  bad_config <- toolchain$key %in% required_config & toolchain$status != 0L
  missing_library_dirs <- library_dirs[!library_dirs$exists, , drop = FALSE]
  preflight_failures <- character()
  if (any(bad_config)) {
    preflight_failures <- c(
      preflight_failures,
      paste0(
        "R CMD config failed for: ",
        paste(toolchain$key[bad_config], collapse = ", ")
      )
    )
  }
  if (any(!compiler_inventory$exists)) {
    preflight_failures <- c(
      preflight_failures,
      paste0(
        "Compiler executable not found for: ",
        paste(compiler_inventory$key[!compiler_inventory$exists], collapse = ", ")
      )
    )
  }
  if (nrow(missing_library_dirs)) {
    preflight_failures <- c(
      preflight_failures,
      paste0(
        "Link-library directories do not exist: ",
        paste(unique(missing_library_dirs$path), collapse = ", ")
      )
    )
  }
  preflight_time <- lhwl_utc_time()
  preflight_detail <- if (length(preflight_failures)) {
    paste(preflight_failures, collapse = "; ")
  } else {
    "required compiler and link-library paths exist"
  }
  steps <- rbind(steps, data.frame(
    step = "toolchain_preflight",
    status = if (length(preflight_failures)) "FAILED" else "PASSED",
    started_utc = preflight_time,
    ended_utc = preflight_time,
    elapsed_seconds = 0,
    detail = preflight_detail,
    stringsAsFactors = FALSE
  ))
  if (length(preflight_failures)) {
    audit_status <- "FAILED"
    failure <- preflight_detail
    stop("Toolchain preflight failed: ", failure, call. = FALSE)
  }

  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(pkg_dir)

  run_step <- function(label, expression) {
    result <- lhwl_capture_step(label, expression, log_dir)
    steps <<- rbind(steps, result$record)
    if (!identical(result$record$status, "PASSED")) {
      audit_status <<- "FAILED"
      failure <<- result$record$detail
      stop("Audit step failed: ", label, ": ", failure, call. = FALSE)
    }
    result$value
  }

  if (isTRUE(document)) {
    message("1. Documenting package...")
    run_step("document", devtools::document())
  } else {
    steps <- rbind(steps, data.frame(
      step = "document", status = "SKIPPED", started_utc = lhwl_utc_time(),
      ended_utc = lhwl_utc_time(), elapsed_seconds = 0,
      detail = "document = FALSE", stringsAsFactors = FALSE
    ))
  }

  message("2. Running tests...")
  run_step("test", devtools::test(stop_on_failure = TRUE, reporter = "summary"))

  if (isTRUE(run_examples)) {
    message("3. Running examples...")
    run_step("examples", devtools::run_examples())
  } else {
    steps <- rbind(steps, data.frame(
      step = "examples", status = "SKIPPED", started_utc = lhwl_utc_time(),
      ended_utc = lhwl_utc_time(), elapsed_seconds = 0,
      detail = "run_examples = FALSE", stringsAsFactors = FALSE
    ))
  }

  if (isTRUE(run_check)) {
    message("4. Running devtools::check()...")
    check <- run_step(
      "check",
      devtools::check(document = FALSE, args = check_args, error_on = "never")
    )
    findings <- c(
      errors = length(check$errors),
      warnings = length(check$warnings),
      notes = length(check$notes)
    )
    steps$detail[steps$step == "check"] <- paste(
      names(findings), findings, collapse = "; "
    )
    if (findings[["errors"]] > 0L) {
      audit_status <- "FAILED"
      failure <- steps$detail[steps$step == "check"]
      stop("R CMD check reported errors: ", failure, call. = FALSE)
    }
    if (sum(findings) > 0L) {
      audit_status <- "PASSED_WITH_FINDINGS"
    }
  } else {
    message("4. Skipping devtools::check(); set run_check = TRUE before release/revision.")
    steps <- rbind(steps, data.frame(
      step = "check", status = "SKIPPED", started_utc = lhwl_utc_time(),
      ended_utc = lhwl_utc_time(), elapsed_seconds = 0,
      detail = "run_check = FALSE", stringsAsFactors = FALSE
    ))
  }

  if (identical(audit_status, "RUNNING")) {
    audit_status <- "PASSED"
  }

  invisible(list(
    status = audit_status,
    output_dir = output_dir,
    steps = steps,
    missing_link_library_dirs = library_dirs[!library_dirs$exists, , drop = FALSE]
  ))
}

message("Loaded run_lhwr_review_audit(). Example:")
message(paste0(
  "run_lhwr_review_audit(pkg_dir = 'LHWL_package', run_check = FALSE, ",
  "document = FALSE)"
))
