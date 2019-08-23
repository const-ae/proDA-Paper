#!/bin/bash

# Remove rows that come from Match-between-runs
cat data/cox_proteome_benchmark/evidence.txt | \
  grep --invert-match MULTI-MATCH > \
  compare_performance/tmp/triqler_input/evidence_wo_mbr.txt


# Convert to triqler input format  
# (takes approx. 3 min)
python -m triqler.convert.maxquant \
   --file_list_file compare_performance/tmp/triqler_input/experimentalDesign_mod.txt \
   --skip_normalization \
   --out_file compare_performance/tmp/triqler_input/triqler_input.tsv \
   compare_performance/tmp/triqler_input/evidence_wo_mbr.txt


python -m triqler \
  --fold_change_eval 0 \
  --decoy_pattern REV__  \
  --out_file compare_performance/tmp/triqler_results/triqler_result.tsv \
  compare_performance/tmp/triqler_input/triqler_input.tsv


