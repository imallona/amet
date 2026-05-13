## Authors: Atreya Choudhury, Izaskun Mallona
## Upstream: https://github.com/imallona/yamet/blob/main/workflow/rules/src/diff_testing.R

#' Differential entropy testing with flexible formula and locations.
#'
#' The model fits i_total against meth + I(meth^2) + loc + patient per window
#' using rowwise lm and limma::squeezeVar moderation. The 'i_total' / 'meth'
#' column names are kept internally so the formula text passed by callers
#' stays unchanged.
#'
#' @param sub_i_total matrix of per-cell i_total values (rows = regions, cols = cells)
#' @param sub_meth   matrix of per-cell methylation values (same dims)
#' @param groups      data.frame with columns: subloc, patient (length = ncol)
#' @param formula     model formula (string or formula)
#' @param loc_levels  character vector of two location codes to compare
#' @param ref_level   reference level for loc factor
#' @param contrast    coefficient name to extract
#' @param param       BiocParallel parameter object
#' @param top_n       number of top regions to return
#' @param out_file    file path to save coefs/ results (RDS)
#'
#' @return list with coefs_df, top_entropy, top_meth
diff_entropy_test <- function(sub_i_total, sub_meth, groups,
                              formula = "i_total ~ meth + I(meth^2) + loc + patient",
                              loc_levels = c("PT","NC"),
                              ref_level = "NC",
                              contrast = "locPT",
                              param,
                              top_n = 2000,
                              out_file = "diff_entropy_coefs.rds") {

  if (is.character(formula)) formula <- as.formula(formula)

  rowwise_lm <- function(i) {
    df <- data.frame(
      i_total = sub_i_total[i, ],
      meth = sub_meth[i, ],
      loc = factor(substr(groups$subloc, 1, 2)),
      patient = groups$patient
    )
    df <- df[df$loc %in% loc_levels, ]
    df <- df[stats::complete.cases(df), ]
    if (nrow(df) == 0 || !(ref_level %in% df$loc)) return(rep(NA_real_, 5))

    df$loc <- relevel(factor(df$loc, levels = loc_levels), ref = ref_level)

    fit <- try(lm(formula, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(rep(NA_real_, 5))

    s <- summary(fit)$coefficients
    if (!(contrast %in% rownames(s))) return(rep(NA_real_, 5))

    c(
      estimate = s[contrast, "Estimate"],
      std_error = s[contrast, "Std. Error"],
      t_value = s[contrast, "t value"],
      p_value = s[contrast, "Pr(>|t|)"],
      df = df.residual(fit)
    )
  }

  coefs_list <- BiocParallel::bplapply(seq_len(nrow(sub_i_total)), rowwise_lm,
                                       BPPARAM = param)
  coefs_df <- as.data.frame(do.call(rbind, coefs_list))
  colnames(coefs_df) <- c("estimate","std_error","t_value","p_value","df")
  coefs_df[] <- lapply(coefs_df, as.numeric)

  coefs_df$moderated_t <- NA_real_
  coefs_df$moderated_p <- NA_real_
  coefs_df$adj_p <- NA_real_

  valid <- complete.cases(coefs_df[, c("estimate","std_error","t_value","p_value","df")])
  ## squeezeVar requires at least one valid row; with sparse proto data every
  ## row can be NA and we then leave moderated_* / adj_p as NA.
  if (any(valid)) {
    coefs_valid <- coefs_df[valid, ]
    squeezed <- limma::squeezeVar(var = coefs_valid$std_error^2,
                                  df = coefs_valid$df)
    moderated_t <- coefs_valid$estimate / sqrt(squeezed$var.post)
    moderated_p <- 2 * pt(-abs(moderated_t), df = squeezed$df.prior + coefs_valid$df)
    adj_p <- p.adjust(moderated_p, method = "BH")
    coefs_df[valid, c("moderated_t","moderated_p","adj_p")] <-
      cbind(moderated_t, moderated_p, adj_p)
  }

  coefs_df$region <- rownames(sub_i_total)

  saveRDS(coefs_df, file = out_file)

  sorted_idx <- order(coefs_df$adj_p, na.last = NA)
  top_idx <- head(sorted_idx, top_n)
  top_entropy <- sub_i_total[top_idx, , drop = FALSE]
  top_meth <- sub_meth[top_idx, , drop = FALSE]

  list(
    coefs_df = coefs_df,
    top_entropy = top_entropy,
    top_meth = top_meth
  )
}
