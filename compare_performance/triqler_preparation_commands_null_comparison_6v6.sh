#!/bin/bash
module load Python/3.6.4-foss-2017b

# Remove rows that come from Match-between-runs
cat data/degraaf_null_comparison/evidence.txt | \
    grep --invert-match MULTI-MATCH | \
    grep -e "Raw file\|OR8_121123_EDG_PG_phosJurkat_0min_A3\|OR8_121123_EDG_PG_phosJurkat_0min_B1\|OR8_121123_EDG_PG_phosJurkat_0min_B2\|OR8_121123_EDG_PG_phosJurkat_0min_B3\|OR8_121123_EDG_PG_phosJurkat_0min_C1\|OR8_121123_EDG_PG_phosJurkat_0min_C2\|OR9_20121214_EDG_PG_phosJurkat_0min_A1\|OR9_20121214_EDG_PG_phosJurkat_0min_A2\|OR9_20121214_EDG_PG_phosJurkat_0min_A3\|OR9_20121214_EDG_PG_phosJurkat_0min_B1\|OR9_20121214_EDG_PG_phosJurkat_0min_B2\|OR9_20121214_EDG_PG_phosJurkat_0min_C3" > \
    compare_performance/tmp/de_graaf/triqler_input/evidence_wo_mbr_6v6.txt

# Convert to triqler input format  
# (takes approx. 3 min)
python -m triqler.convert.maxquant \
   --file_list_file compare_performance/tmp/de_graaf/triqler_input/experimentalDesign_mod_6v6.txt \
   --skip_normalization \
   --out_file compare_performance/tmp/de_graaf/triqler_input/triqler_input_6v6.tsv \
   compare_performance/tmp/de_graaf/triqler_input/evidence_wo_mbr_6v6.txt


python -m triqler \
  --num_threads 12 \
  --fold_change_eval 0 \
  --decoy_pattern REV__  \
  --out_file compare_performance/tmp/de_graaf/triqler_results/triqler_result_6v6.tsv \
  compare_performance/tmp/de_graaf/triqler_input/triqler_input_6v6.tsv


