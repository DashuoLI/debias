# ------------------------------
# Statistical Utilities
# ------------------------------
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(effectsize))
suppressPackageStartupMessages(library(infotheo))

safe_ttest <- function(x, y) {
  if (length(unique(x)) > 1 && length(unique(y)) > 1) {
    return(t.test(x, y)$p.value)
  }
  NA_real_
}

safe_wilcox <- function(x, y) {
  if (length(unique(x)) > 1 && length(unique(y)) > 1) {
    return(wilcox.test(x, y)$p.value)
  }
  NA_real_
}

safe_shapiro <- function(x) {
  if (length(unique(x)) > 1) {
    return(shapiro.test(x)$p.value)
  }
  NA_real_
}

safe_kruskal <- function(value, group){
  if (min(table(group)) > 1){ 
    return(kruskal.test(value ~ group)$p.value)
  }
  NA_real_
}

safe_anova <- function(value, group){
  if (min(table(group)) > 1){
    anova_model <- aov(value ~ group)
    return(summary(anova_model)[[1]][["Pr(>F)"]][1])
  }
  NA_real_
}

#safe_fisher_pitman <- function(value, group){
#  if (min(table(group)) > 1){
#    group <- factor(group)
#    res <- oneway_test(value ~ group, distribution = approximate(nresample = 10000))
#    return(pvalue(res))
#  }
#  NA_real_
#}

safe_fisher_pitman <- function(value, group) {
  # remove NAs
  keep <- !(is.na(value) | is.na(group))
  value <- value[keep]
  group <- group[keep]
  
  # must have at least 2 groups with > 1 obs each
  if (length(unique(group)) < 2L) return(NA_real_)
  if (min(table(group)) <= 1L) return(NA_real_)
  
  # must have variation in data
  if (all(tapply(value, group, function(x) var(x) == 0, simplify = TRUE))) {
    return(NA_real_)
  }
  
  # try-catch in case coin still fails
  out <- try(
    {
      group <- factor(group)
      res <- oneway_test(value ~ group, distribution = coin::approximate(nresample = 10000))
      pvalue(res)
    },
    silent = TRUE
  )
  
  if (inherits(out, "try-error")) NA_real_ else out
}


safe_cohens_d <- function(x, y, pooled_sd = TRUE) {
  sd_x <- sd(x)
  sd_y <- sd(y)
  
  # Case 1: both groups constant and equal
  if (sd_x == 0 && sd_y == 0 && mean(x) == mean(y)) {
    return(0)  # no difference
  }
  
  # Case 2: both groups constant but different
  if (sd_x == 0 && sd_y == 0 && mean(x) != mean(y)) {
    return(Inf)  # perfectly separated
  }
  
  # Case 3: one group has zero variance
  if (sd_x == 0 || sd_y == 0) {
#    warning("One group has zero variance; Cohen's d may not be meaningful.")
    return(NA)
  }
  
  # Otherwise: compute Cohen's d normally
  res <- cohens_d(x, y, pooled_sd = pooled_sd)
  return(res$Cohens_d)  # extract numeric value
}


safe_cohens_f <- function(y, group) {
  # Ensure group is a factor
  group <- as.factor(group)
  n_groups <- nlevels(group)

  # Edge case: only one group
  if (n_groups < 2) {
    warning("Less than 2 groups; Cohen's f is undefined.")
    return(NA)
  }

  # Compute group means and SDs
  group_means <- tapply(y, group, mean)
  group_sds <- tapply(y, group, sd)
  
  # Check for all groups constant
  if (all(group_sds == 0)) {
    if (length(unique(group_means)) == 1) {
      return(0)  # no difference
    } else {
      return(Inf)  # perfectly separated
    }
  }
  
  # Check for any group with zero variance
  if (any(group_sds == 0)) {
    # Cohen's f may not be meaningful
    return(NA)
  }
  
  # Fit ANOVA
  fit <- aov(y ~ group)
  ss <- summary(fit)[[1]][, "Sum Sq"]
  ss_between <- ss[1]
  ss_total <- sum(ss)
  
  # Compute eta-squared and Cohen's f
  eta2 <- ss_between / ss_total
  f <- sqrt(eta2 / (1 - eta2))
  
  return(f)
}

