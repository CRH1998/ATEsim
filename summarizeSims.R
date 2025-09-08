summarize_sim_ATE <- function(sim_ests,
                              truths,                 # either named vector by target (e.g., c("treat0"=..,"treat1"=..,"treat:1-0"=..))
                              # or list(model = named vector by target)
                              ci_z = 1.96,
                              name_prefix = "ATEests_",
                              method_tokens = c("coef","G","DR"),
                              collapse_fun = function(x) mean(x, na.rm = TRUE)) {
  ## helpers
  tok_pat <- paste(method_tokens, collapse = "|")
  norm_target <- function(x) {
    x <- sub("^treat:?", "treat", x)
    gsub("\\s+", "", x)
  }
  
  parse_names <- function(nm) {
    # model = anything between prefix and _(coef|G|DR)_
    model  <- sub(paste0("^", name_prefix, "(.+?)_(?:", tok_pat, ")_.*$"), "\\1", nm, perl = TRUE)
    bad_m  <- (model == nm)
    method <- sub(paste0("^", name_prefix, ".+?_(", tok_pat, ")_.*$"), "\\1", nm, perl = TRUE)
    method <- toupper(method)
    ok_method <- grepl(paste0("_(?:", tok_pat, ")_"), nm, perl = TRUE)
    method[!ok_method] <- NA
    
    stat   <- ifelse(grepl("_est_", nm), "est",
                     ifelse(grepl("_se_", nm), "se", NA))
    
    tr_raw <- sub(".*_treat(.*)$", "\\1", nm, perl = TRUE)
    tr_bad <- !grepl("_treat", nm, perl = TRUE)
    tr_raw[tr_bad] <- NA
    target <- ifelse(is.na(tr_raw), NA, paste0("treat", sub("^:", "", tr_raw)))
    
    data.frame(name = nm, model = model, method = method, stat = stat, target = target,
               stringsAsFactors = FALSE)
  }
  
  collapse_wrapper <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    collapse_fun(x)
  }
  
  agg_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  agg_sd   <- function(x) { x <- x[!is.na(x)]; if (length(x) <= 1) return(NA_real_); stats::sd(x) }
  first_non_na <- function(x) { i <- match(TRUE, !is.na(x)); if (is.na(i)) NA else x[i] }
  
  truth_lookup <- function(models, targets) {
    out <- rep(NA_real_, length(models))
    if (is.list(truths)) {
      for (i in seq_along(models)) {
        m <- models[i]; t <- norm_target(targets[i])
        if (!is.na(m) && !is.na(t) && !is.null(truths[[m]]) && !is.null(truths[[m]][t])) {
          out[i] <- suppressWarnings(as.numeric(truths[[m]][t]))
        }
      }
    } else {
      # a single vector by target, applied to all models
      tv <- truths; names(tv) <- norm_target(names(tv))
      for (i in seq_along(targets)) {
        t <- norm_target(targets[i])
        if (!is.na(t) && !is.null(tv[t])) out[i] <- suppressWarnings(as.numeric(tv[t]))
      }
    }
    out
  }
  
  ## 1) Long format
  if (inherits(sim_ests, "data.frame")) {
    cols <- names(sim_ests)[startsWith(names(sim_ests), name_prefix)]
    if (!length(cols)) stop("No columns starting with '", name_prefix, "'.")
    n <- nrow(sim_ests)
    name_vec <- rep(cols, each = n)
    sim_vec  <- rep(seq_len(n), times = length(cols))
    val_vec  <- as.vector(as.matrix(sim_ests[, cols, drop = FALSE]))
    meta <- parse_names(cols)
    meta_rep <- meta[match(name_vec, meta$name), ]
    long <- data.frame(model  = meta_rep$model,
                       method = meta_rep$method,
                       stat   = meta_rep$stat,
                       target = meta_rep$target,
                       .sim   = sim_vec,
                       val    = suppressWarnings(as.numeric(val_vec)),
                       stringsAsFactors = FALSE)
  } else {
    nms <- names(sim_ests); if (is.null(nms)) stop("Named list expected.")
    keep <- startsWith(nms, name_prefix)
    if (!any(keep)) stop("No list entries starting with '", name_prefix, "'.")
    nms  <- nms[keep]
    lens <- lengths(sim_ests[keep])
    name_vec <- rep(nms, times = lens)
    sim_vec  <- unlist(lapply(lens, function(k) seq_len(k)), use.names = FALSE)
    val_vec  <- suppressWarnings(as.numeric(unlist(sim_ests[keep], use.names = FALSE)))
    meta <- parse_names(nms)
    meta_rep <- meta[match(name_vec, meta$name), ]
    long <- data.frame(model  = meta_rep$model,
                       method = meta_rep$method,
                       stat   = meta_rep$stat,
                       target = meta_rep$target,
                       .sim   = sim_vec,
                       val    = val_vec,
                       stringsAsFactors = FALSE)
  }
  
  ok <- !is.na(long$model) & !is.na(long$method) & !is.na(long$stat) & !is.na(long$target)
  long <- long[ok, , drop = FALSE]
  long$target_n <- norm_target(long$target)
  
  ## 2) Collapse duplicates
  agg <- aggregate(val ~ model + method + target_n + .sim + stat,
                   data = long, FUN = collapse_wrapper)
  
  ## 3) Pivot to est/se
  wide <- reshape(agg,
                  idvar = c("model","method","target_n",".sim"),
                  timevar = "stat",
                  direction = "wide")
  names(wide) <- sub("^val\\.", "", names(wide))
  wide <- wide[!is.na(wide$est) & !is.na(wide$se), , drop = FALSE]
  
  ## 4) Truths, CIs, coverage
  wide$true <- truth_lookup(wide$model, wide$target_n)
  wide$LCI  <- wide$est - ci_z * wide$se
  wide$UCI  <- wide$est + ci_z * wide$se
  wide$covered <- ifelse(is.na(wide$true), NA,
                         (wide$LCI < wide$true) & (wide$UCI > wide$true))
  
  ## 5) Summaries
  key <- c("model","method","target_n")
  
  MeanEst <- aggregate(est ~ model + method + target_n, wide, agg_mean); names(MeanEst)[4] <- "MeanEst"
  SD      <- aggregate(est ~ model + method + target_n, wide, agg_sd);   names(SD)[4]      <- "SD"
  MeanSE  <- aggregate(se  ~ model + method + target_n, wide, agg_mean); names(MeanSE)[4]  <- "MeanSE"
  Coverage<- aggregate(covered ~ model + method + target_n, wide, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  names(Coverage)[4] <- "Coverage"
  Truth   <- aggregate(true ~ model + method + target_n, wide, first_non_na); names(Truth)[4] <- "Truth"
  n_sims  <- aggregate(est ~ model + method + target_n, wide, function(x) sum(!is.na(x))); names(n_sims)[4] <- "n_sims"
  
  summ <- Reduce(function(a,b) merge(a,b, by = key, all = TRUE),
                 list(merge(MeanEst, SD, by = key, all = TRUE),
                      merge(MeanSE, Coverage, by = key, all = TRUE),
                      merge(Truth, n_sims, by = key, all = TRUE)))
  names(summ)[names(summ) == "target_n"] <- "target"
  summ <- summ[order(summ$model, summ$target, summ$method), ]
  rownames(summ) <- NULL
  
  list(summary = summ, draws = wide)
}


# ---- end function ------------------------------------------------------------








summarize_sim_coefs <- function(sim_ests,
                                truths,                 # named vector by term, or list(model = named vector by term)
                                ci_z = 1.96,
                                name_prefix = "ATEests_",
                                collapse_fun = function(x) mean(x, na.rm = TRUE)) {
  ## helpers
  parse_names <- function(nm) {
    # model = anything between prefix and _coef_
    model <- sub(paste0("^", name_prefix, "(.+?)_coef_.*$"), "\\1", nm, perl = TRUE)
    bad_m <- (model == nm)
    
    stat  <- ifelse(grepl("_coef_est_", nm), "est",
                    ifelse(grepl("_coef_se_", nm), "se", NA))
    
    term  <- sub(".*_coef_(?:est|se)_(.*)$", "\\1", nm, perl = TRUE)
    bad_t <- !grepl("_coef_(?:est|se)_", nm, perl = TRUE)
    term[bad_t] <- NA
    
    data.frame(name = nm, model = ifelse(bad_m, NA, model), stat = stat, term = term,
               stringsAsFactors = FALSE)
  }
  collapse_wrapper <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    collapse_fun(x)
  }
  agg_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  agg_sd   <- function(x) { x <- x[!is.na(x)]; if (length(x) <= 1) return(NA_real_); stats::sd(x) }
  first_non_na <- function(x) { i <- match(TRUE, !is.na(x)); if (is.na(i)) NA else x[i] }
  
  truth_lookup <- function(models, terms) {
    out <- rep(NA_real_, length(models))
    if (is.list(truths)) {
      for (i in seq_along(models)) {
        m <- models[i]; t <- terms[i]
        if (!is.na(m) && !is.na(t) && !is.null(truths[[m]]) && !is.null(truths[[m]][t])) {
          out[i] <- suppressWarnings(as.numeric(truths[[m]][t]))
        }
      }
    } else {
      for (i in seq_along(terms)) {
        t <- terms[i]
        if (!is.na(t) && !is.null(truths[t])) out[i] <- suppressWarnings(as.numeric(truths[t]))
      }
    }
    out
  }
  
  ## 1) Long
  if (inherits(sim_ests, "data.frame")) {
    cols <- names(sim_ests)[startsWith(names(sim_ests), name_prefix)]
    if (!length(cols)) stop("No columns starting with '", name_prefix, "'.")
    n <- nrow(sim_ests)
    name_vec <- rep(cols, each = n)
    sim_vec  <- rep(seq_len(n), times = length(cols))
    val_vec  <- as.vector(as.matrix(sim_ests[, cols, drop = FALSE]))
    meta <- parse_names(cols)
    meta_rep <- meta[match(name_vec, meta$name), ]
    long <- data.frame(model = meta_rep$model,
                       stat  = meta_rep$stat,
                       term  = meta_rep$term,
                       .sim  = sim_vec,
                       val   = suppressWarnings(as.numeric(val_vec)),
                       stringsAsFactors = FALSE)
  } else {
    nms <- names(sim_ests); if (is.null(nms)) stop("Named list expected.")
    keep <- startsWith(nms, name_prefix)
    if (!any(keep)) stop("No list entries starting with '", name_prefix, "'.")
    nms  <- nms[keep]
    lens <- lengths(sim_ests[keep])
    name_vec <- rep(nms, times = lens)
    sim_vec  <- unlist(lapply(lens, function(k) seq_len(k)), use.names = FALSE)
    val_vec  <- suppressWarnings(as.numeric(unlist(sim_ests[keep], use.names = FALSE)))
    meta <- parse_names(nms)
    meta_rep <- meta[match(name_vec, meta$name), ]
    long <- data.frame(model = meta_rep$model,
                       stat  = meta_rep$stat,
                       term  = meta_rep$term,
                       .sim  = sim_vec,
                       val   = val_vec,
                       stringsAsFactors = FALSE)
  }
  
  ok <- !is.na(long$model) & !is.na(long$stat) & !is.na(long$term)
  long <- long[ok, , drop = FALSE]
  
  ## 2) Collapse duplicates
  agg <- aggregate(val ~ model + term + .sim + stat,
                   data = long, FUN = collapse_wrapper)
  
  ## 3) Pivot to est/se
  wide <- reshape(agg,
                  idvar = c("model","term",".sim"),
                  timevar = "stat",
                  direction = "wide")
  names(wide) <- sub("^val\\.", "", names(wide))
  wide <- wide[!is.na(wide$est) & !is.na(wide$se), , drop = FALSE]
  
  ## 4) Truths, CIs, coverage
  wide$true <- truth_lookup(wide$model, wide$term)
  wide$LCI  <- wide$est - ci_z * wide$se
  wide$UCI  <- wide$est + ci_z * wide$se
  wide$covered <- ifelse(is.na(wide$true), NA,
                         (wide$LCI < wide$true) & (wide$UCI > wide$true))
  
  ## 5) Summaries
  key <- c("model","term")
  MeanEst <- aggregate(est ~ model + term, wide, agg_mean); names(MeanEst)[3] <- "MeanEst"
  SD      <- aggregate(est ~ model + term, wide, agg_sd);   names(SD)[3]      <- "SD"
  MeanSE  <- aggregate(se  ~ model + term, wide, agg_mean); names(MeanSE)[3]  <- "MeanSE"
  Coverage<- aggregate(covered ~ model + term, wide, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  names(Coverage)[3] <- "Coverage"
  Truth   <- aggregate(true ~ model + term, wide, first_non_na); names(Truth)[3] <- "Truth"
  n_sims  <- aggregate(est ~ model + term, wide, function(x) sum(!is.na(x))); names(n_sims)[3] <- "n_sims"
  
  summ <- Reduce(function(a,b) merge(a,b, by = key, all = TRUE),
                 list(merge(MeanEst, SD, by = key, all = TRUE),
                      merge(MeanSE, Coverage, by = key, all = TRUE),
                      merge(Truth, n_sims, by = key, all = TRUE)))
  summ <- summ[order(summ$model, summ$term), ]
  rownames(summ) <- NULL
  
  list(summary = summ, draws = wide)
}

# ---- end function ------------------------------------------------------------








# ---- Example usage -----------------------------------------------------------
# Provide the true values once. Name flexibility is allowed:
# truths <- c(`treat0` = true_treat0, `treat1` = true_treat1, `treat:1-0` = true_treat01)
# res <- summarize_sim_ATE(sim_ests, truths = truths)
# 
# 
# 
# truths <- c(Intercept = -1, Z21 = 0.3, Z31 = -0.3)
# 
# res <- summarize_sim_coefs(sim_ests, truths = truths)
# res$summary

# View numeric summary table
# res$summary
# print(res$plot_dist)
# print(res$plot_summary)
