

run_perseus <- function(data, tag="", return_null_if_missing = FALSE, dataset_name = "de_graaf") {
  
  hash <- digest::digest(data)  
  prefix_in <- paste0("tmp/", dataset_name, "/perseus/input")
  prefix_out <- paste0("tmp/", dataset_name, "/perseus/output")
  table_file_name <- paste0("perseus_data_table_", tag, "-hash_", hash, ".tsv")
  result_file_name <- paste0("perseus_result_table_", tag, "-hash_", hash, ".txt")
  if(!file.exists(file.path(prefix_in, table_file_name))){
    data %>%
      as.data.frame() %>%
      rownames_to_column("Name") %>%
      write_tsv(file.path(prefix_in, table_file_name))
  }
  if(dataset_name == "de_graaf"){
    pval_column_name <- "Student's T-test p-value Syn_Cond1_Syn_Cond2"
    qval_column_name <- "Student's T-test q-value Syn_Cond1_Syn_Cond2"
    diff_column_name <- "Student's T-test Difference Syn_Cond1_Syn_Cond2"
  }else if(dataset_name == "spikein_benchmark"){
    pval_column_name <- "Student's T-test p-value H_L"
    qval_column_name <- "Student's T-test q-value H_L"
    diff_column_name <- "Student's T-test Difference H_L"
  }else{
    stop(paste0("Illegal dataset_name: ", dataset_name))
  }
  if(!file.exists(file.path(prefix_out, result_file_name))){
    warning(paste0("Couldn't find perseus result file: ",
                   file.path(prefix_out, result_file_name)))
    if(return_null_if_missing){
      return(NULL)
    }
    break_loop <- FALSE
    while(! break_loop){
      resp <- readline(prompt="Abort? (y/n)")
      if(resp == "y"){
        stop("No result from Perseus available yet")
      }else{
        if(!file.exists(file.path(prefix_out, result_file_name))){
          warning(paste0(file.path(prefix_out, result_file_name), " still does not exist"))
        }else{
          break_loop <- TRUE
        }
      }
    }
    read_tsv(paste0(prefix_out, "/", result_file_name), 
             comment = "#") %>%
      mutate(pval = as.numeric(str_replace(!! sym(pval_column_name), ",", "."))) %>%
      mutate(adj_pval = as.numeric(str_replace(!! sym(qval_column_name), ",", "."))) %>%
      mutate(diff = as.numeric(str_replace(!! sym(diff_column_name), ",", "."))) %>%
      transmute(name=Name, pval, adj_pval, diff) 
  }else{
    read_tsv(paste0(prefix_out, "/", result_file_name), 
             comment = "#") %>%
      mutate(pval = as.numeric(str_replace(!! sym(pval_column_name), ",", "."))) %>%
      mutate(adj_pval = as.numeric(str_replace(!! sym(qval_column_name), ",", "."))) %>%
      mutate(diff = as.numeric(str_replace(!! sym(diff_column_name), ",", "."))) %>%
      transmute(name=Name, pval, adj_pval, diff) 
  }
}


