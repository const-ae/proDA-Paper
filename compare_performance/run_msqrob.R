


run_msqrob_on_cox_benchmark <- function(){
  
  library(MSqRob)
  library(MSnbase)
  
  mse <- MSqRob::import2MSnSet("../data/cox_proteome_benchmark/peptides.txt",
                               filetype = "MaxQuant")
  
  
  exprs(mse) <- ifelse(exprs(mse) == 0, NA, log2(exprs(mse)))
    
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
  
  annotated_fdata <- fData(mse) %>%
    as_tibble() %>%
    mutate(id = seq_len(n())) %>%
    mutate(Protein = str_split(Proteins, ";")) %>%
    unnest() %>%
    left_join(all_prots_species_assignment, by="Protein") %>%
    dplyr::select(-c(DB, Description)) %>%
    group_by_at(vars(-c(Protein, Origin))) %>%
    summarize(Origin = case_when(
        all(is.na(Origin)) ~ "Unknown",
        all(Origin[!is.na(Origin)] == "HS") ~ "HS",
        all(Origin[!is.na(Origin)] == "EC") ~ "EC",
        any(Origin[!is.na(Origin)] == "HS") & any(Origin[!is.na(Origin)] == "EC") ~ "Both",
        TRUE ~ "Bug"
    )) %>%
    as.data.frame()
  
  rownames(annotated_fdata) <- rownames(mse)
  fData(mse) <- as.data.frame(annotated_fdata)
  
  mse_norm <- proDA::median_normalization(mse, fData(mse)$Origin == "HS")
  
  runs <- sampleNames(mse_norm)
  experimental_design <- str_sub(runs, 1,1)
  
  exp_annotation <- data.frame(run = runs, experimental_design=experimental_design,
                               row.names = runs, stringsAsFactors = FALSE)
  pData(mse_norm) <- exp_annotation
  
  proteins_mse_norm <- MSnSet2protdata(mse_norm, accession = "Proteins",
                                       printProgress = TRUE)
  # Filter out empty proteins
  proteins_mse_norm_filtered <- proteins_mse_norm[which(sapply(seq_len(length(proteins_mse_norm)), function(x) nrow(getData(proteins_mse_norm[x]))) != 0)]
  
  prot_fit <- fit.model(proteins_mse_norm_filtered, response = "quant_value",
                        fixed = "experimental_design", random = c("run", "Sequence"),
                        add.intercept = TRUE, weights="Huber")
  
  cntrst <- makeContrast("experimental_designH - experimental_designL", c("experimental_designH", "experimental_designL"))
  prot_result <- test.protLMcontrast(prot_fit, cntrst)
  prot_result_annotated <- prot_result %>% 
    as_tibble(rownames = "name") %>%
    mutate(adj_pval = p.adjust(pval, "BH")) %>%
    mutate(diff = estimate) %>%
    mutate(id = seq_len(n())) %>%
    mutate(Protein = str_split(name, ";")) %>%
    unnest() %>%
    left_join(all_prots_species_assignment, by = "Protein") %>%
    dplyr::select(-c(DB, Description)) %>%
    group_by_at(vars(-c(Protein, Origin))) %>%
    summarize(Origin = case_when(
      all(is.na(Origin)) ~ "Unknown",
      all(Origin[!is.na(Origin)] == "HS") ~ "HS",
      all(Origin[!is.na(Origin)] == "EC") ~ "EC",
      any(Origin[!is.na(Origin)] == "HS") & any(Origin[!is.na(Origin)] == "EC") ~ "Both",
      TRUE ~ "Bug"
    )) 
  
  prot_result_annotated %>%
    ungroup() %>%
    filter(Origin == "HS" | Origin == "EC") %>%
    transmute(name, pval, adj_pval, diff,
              organism_short = ifelse(Origin == "HS", "H. sapiens", "E. coli"))
}
