

run_qprot <- function(X, Y,
                      nburnin=2000,
                      niter=10000,
                      normalize=TRUE,
                      paired=FALSE) {
  stopifnot(rownames(X) == rownames(Y))
  # browser()
  data_mat <- as.data.frame(cbind(rownames(X), X, Y))
  colnames(data_mat) <- c("Protein", rep(0, ncol(X)), rep(1, ncol(Y)))
  
  tf <- tempfile()
  write.table(data_mat, file=tf, sep="\t", quote=FALSE, row.names = FALSE)
  
  if(paired){
    program <- "qprot-paired"
  }else{
    program <- "qprot-param"
  }
  
  if(!file.exists("/home/ahlmann/prog/qprot_v1.3.5")){
    stop("qprot needs to be installed. I had it in the /home/ahlmann/prog/qprot_v1.3.5 folder ", 
         "but you can put it anywhere just remember to edit the ",
         "compare_performance/run_qprot.R script.")
  }
  system2(paste0("/home/ahlmann/prog/qprot_v1.3.5/bin/", program),
          c(tf, nburnin, niter, normalize * 1.0))
  system2("/home/ahlmann/prog/qprot_v1.3.5/bin/getfdr",
          paste0(tf, "_qprot"))
  read.delim(paste0(tf, "_qprot_fdr"))
}
