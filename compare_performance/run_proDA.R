

run_proDA <- function(data, experimental_design, verbose=FALSE){
  
  fit <- proDA(data, design = experimental_design,
               n_subsample = 1000, verbose=verbose)
  groups <- unique(experimental_design)
  result <- proDA::test_diff(fit, contrast = paste0(groups[1], " - ", groups[2]))
  
  return(list(
    result = result,
    fit = fit
  ))
}

