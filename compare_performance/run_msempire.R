



run_msempire_on_cox_benchmark <- function(){
  library(msEmpiRe)

  # Make experimental design file
  tribble(
    ~Name, ~Condition,
    "H1",   "H",
    "H2",   "H",
    "H3",   "H",
    "L1",   "L",
    "L2",   "L",
    "L3",   "L",
  ) %>%
    write_tsv("/tmp/condition_sample_mapping.tsv")
  
  # Load Data
  msemp <- read.MaxQuant("../data/cox_proteome_benchmark/peptides.txt", "/tmp/condition_sample_mapping.tsv")
  
  # Filter out stuff based on the proteinGroups.txt file
  msemp <- filter_MaxQuant(msemp, "../data/cox_proteome_benchmark/proteinGroups.txt", filter_on = c("Contaminant", "Reverse"))
  
  # do log transformation
  msemp_cp <- msemp
  exprs(msemp_cp) <- ifelse(exprs(msemp_cp) == 0, NA, log2(exprs(msemp_cp)))
  # Assign correct species
  org_assig <- fData(msemp_cp) %>% 
    as_tibble() %>%
    dplyr::select(Sequence, prot.id) %>%
    left_join(species_assignment, by=c("prot.id"="ProteinID")) %>%
    pull(organism_short)
  # Normalize data based on species
  msemp_norm3 <- proDA::median_normalization(msemp_cp, org_assig == "H. sapiens")
  fData(msemp_norm3)$Origin <- org_assig
  
  # Transform back to un-logarithmized scale
  msemp_cp2 <- msemp_norm3
  exprs(msemp_cp2) <- ifelse(is.na(exprs(msemp_cp2)), 0, 2^(exprs(msemp_cp2)))
  # Filter out proteins that are completely missing in one condition
  cond_compl_mis <- rowSums(exprs(msemp_cp2)[,pData(msemp_cp2)$condition == "H"] == 0) == 3 |
    rowSums(exprs(msemp_cp2)[,pData(msemp_cp2)$condition == "L"] == 0) == 3
  
  # Do the actual test
  dir.create("/tmp/msempire_plots")
  msemp_res <- de.ana(msemp_cp2[! cond_compl_mis, ], out.dir = "/tmp/msempire_plots/")
  
  
  as_tibble(msemp_res) %>%
    transmute(name = prot.id, pval = p.val, adj_pval = p.adj, diff=-log2FC) %>%
    left_join(species_assignment, by=c("name"="ProteinID"))
    
}
