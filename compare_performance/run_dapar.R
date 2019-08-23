run_dapar <- function(data, experimental_design, imputation_method, ...) {

  pData <- AnnotatedDataFrame(data.frame(Condition = experimental_design,
                     Bio.Rep = seq_along(experimental_design),
                     Tech.Rep = seq_along(experimental_design),
                     row.names = colnames(data)))
  fData <- AnnotatedDataFrame(data.frame(bla = seq_len(nrow(data)),
                                         row.names = rownames(data)))
  mse <- MSnbase::MSnSet(data, fData = fData, pData = pData)
  
  imputed_mse <- wrapper.impute.slsa(mse)  # Bo et al., 2004. The method mentioned by the reviewer
  
  groups <- unique(experimental_design)
  prefix <- paste0(groups[1], "_vs_", groups[2])
  suppressWarnings(
    pre_dapar_res <- limmaCompleteTest(exprs(imputed_mse), pData(imputed_mse))
  )
  dapar_res <- as_tibble(as.data.frame(pre_dapar_res), rownames = "name") %>%
    dplyr::select(name, diff = !! sym(paste0(prefix, "_logFC")), pval = !! sym(paste0(prefix, "_pval")))
  
  dapar_res %>%
    mutate(adj_pval = p.adjust(pval, "BH"))
  
}