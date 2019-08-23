

run_dep <- function(data, experimental_design, imputation_method, ...) {
  
  dep_col_annot <- data.frame(ID = colnames(data),
                              label = colnames(data),
                              condition = experimental_design, 
                              replicate = proDD:::as_replicate(experimental_design),
                              stringsAsFactors = FALSE) 
  
  dep_row_annot <- data.frame(ID=rownames(data),
                              name=rownames(data),
                              stringsAsFactors = FALSE) 
  
  se_dep <-  SummarizedExperiment(data, colData = dep_col_annot, rowData = dep_row_annot)
  
  se_dep_imp <- DEP::impute(se_dep, fun = imputation_method, ...)
  groups <- unique(experimental_design)
  
  prefix <- paste0(groups[1], "_vs_", groups[2])
  suppressWarnings({
    DEP::test_diff(se_dep_imp, type = "manual", test = prefix) %>%
      rowData() %>%
      as_tibble() %>%
      arrange(name) %>%
      transmute(name, diff = !! sym(paste0(prefix, "_diff")), 
                adj_pval = p.adjust( !! sym(paste0(prefix, "_p.val")), "BH"),
                pval = !! sym(paste0(prefix, "_p.val")), imputed, num_NAs) 
    
  })
}
