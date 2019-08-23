


parse_triqler_results_for_cox_benchmark <- function(){

  
  
  list <- read_lines("tmp/triqler_results/triqler_result.tsv") %>%
    str_split("\t") 
  headers <- list[[1]]  
  
  triqler_result_df_raw <- list[-1] %>%
    map_dfr(function(chr){
      c(set_names(as.list(chr[which(headers != "peptides")]), headers[headers != "peptides"]),
        list(peptides = paste0(chr[-which(headers != "peptides")], collapse=";")))
    }) %>%
    transmute(protein = protein, peptides = peptides,
              q_value = parse_number(q_value),
              posterior_error_prob = parse_number(posterior_error_prob),
              num_peptides = parse_number(num_peptides),
              log2_fold_change = parse_number(log2_fold_change),
              diff_exp_prob_0.0 = parse_number(diff_exp_prob_0.0))
  
  all_prots_species_assignment <- bind_rows(
    read_lines("../data/cox_proteome_benchmark/species_fastas/UP000005640_9606.fasta") %>%
      keep(function(x) str_starts(x, ">")) %>%
      str_split_fixed("\\|", n=3) %>%
      as_tibble() %>%
      dplyr::rename_all(~ c("DB", "Protein", "Description")) %>%
      mutate(DB  = str_sub(DB, 2)) %>%
      mutate(Origin = "HS"),
    read_lines("../data/cox_proteome_benchmark/species_fastas/UP000000625_83333.fasta") %>%
      keep(function(x) str_starts(x, ">")) %>%
      str_split_fixed("\\|", n=3) %>%
      as_tibble() %>%
      dplyr::rename_all(~ c("DB", "Protein", "Description")) %>%
      mutate(DB  = str_sub(DB, 2)) %>%
      mutate(Origin = "EC"))
  
  triqler_res <- triqler_result_df_raw %>%
    left_join(all_prots_species_assignment, by=c(protein = "Protein")) %>%
    mutate(Origin = case_when(
      str_starts(protein, "REV__") ~ "Decoy",
      str_starts(protein, "CON__") ~ "Contamination",
      TRUE ~ Origin
    )) %>%
    filter(Origin != "Decoy" & Origin != "Contamination")  
    
  triqler_res %>%
    ungroup() %>%
    filter(Origin == "HS" | Origin == "EC") %>%
    transmute(name = protein, pval = NA, adj_pval = q_value, 
              diff = -log2_fold_change,
              organism_short = ifelse(Origin == "HS", "H. sapiens", "E. coli"))
    
}