# 2023_Bat1Kimmunity
Customs scripts used for analyses of 2022-2023 Bat1K immunity project
No new software was developed for this study, thus we provide example commands for analyses or provide links to other sources employed.


## analyses
- Genome assembly

- Repeat masking

- Pairwise genome alignemnts
  https://github.com/hillerlab/GenomeAlignmentTools

- Genome annotation using TOGA
  https://github.com/hillerlab/TOGA
  toga.py ${chains.ref.query} ${annotation_ref.bed} ${genome_ref.2bit} ${genome_query.2bit} -i ${isoforms.-reftxt} --cb 3,5 --cjn 500 --u12 ${U12sites_ref.tsv} --ms
  
 
- Phylogenetic inference
  - Gene trees
  - Species trees with ASTRAL
  - Species trees with iq-tree


- Selection analyses
  - HYPHY

- Enrichment analyses
  - gprofiler

- Linear regression



