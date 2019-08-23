run_ebrct <- function(data, experimental_design) {
  
  groups <- unique(experimental_design)
  
  Y1 <- as.data.frame(data[, experimental_design == groups[1]])
  Y2 <- as.data.frame(data[, experimental_design == groups[2]])
  
  too_many_missing <- rowSums(is.na(data[, experimental_design == groups[1]])) == sum(experimental_design == groups[1]) | 
    rowSums(is.na(data[, experimental_design == groups[2]]))  == sum(experimental_design == groups[2])
  res <- RanCenHier2Fun_mod(Y1[! too_many_missing, ], Y2[! too_many_missing, ])
  
  diff_mean <- res$post.mu.x1.mean - res$post.mu.x2.mean
  pooled_var <- 1/(1/res$post.mu.x1.sd^2 + 1/res$post.mu.x2.sd^2)
  pvals <- pnorm(diff_mean / sqrt(pooled_var))
  pvals <- pmin(pvals, 1-pvals) * 2
  tibble(diff = diff_mean,
         pval = pvals,
         adj_pval = p.adjust(pval, method="BH"),
         name = rownames(Y1[! too_many_missing, ]))
  
  
}