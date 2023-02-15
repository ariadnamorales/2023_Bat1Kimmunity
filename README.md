# 2023_Bat1Kimmunity
Customs scripts used for analyses of 2022-2023 Bat1K immunity project
No new software was developed for this study, thus we provide example commands for analyses or provide links to other sources employed.


## Analyses
#### - Genome assembly

#### - Repeat masking

#### - Pairwise genome alignemnts
  https://github.com/hillerlab/GenomeAlignmentTools

#### - Genome annotation using [TOGA](https://github.com/hillerlab/TOGA)
  ```
  toga.py ${chains.ref.query} ${annotation_ref.bed} ${genome_ref.2bit} ${genome_query.2bit} -i ${isoforms.-reftxt} --cb 3,5 --cjn 500 --u12 ${U12sites_ref.tsv} --ms
  ```
#### - Exon-by-exon codon alignments and cleaning using [extract_codon_alignment.py](https://github.com/hillerlab/TOGA/blob/master/supply/extract_codon_alignment.py) and [HmmCleaner](https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner/view/bin/HmmCleaner.pl)
   ```
   ## exon-by-exon alignment
   extract_codon_alignments_from_toga.py ${f_listSP} ${annotation_ref.bed} ${trasncriptID} -o ${trasncriptID}_raw.fa -s
   
   ## cleaning
   HmmCleaner.pl -costs ${c1} ${c2} ${c3} ${c4} ${transcriptID}_raw.fa
   ```
  
#### - Phylogenetic and Divergence Time Estimation
  - Gene trees using [raxml](https://cme.h-its.org/exelixis/web/software/raxml/)
  ```
  raxmlHPC-PTHREADS -T ${nThreads} -s ${transcriptID}.fa -m ${sustMod} -N ${reps} -p ${seedSearch} -w ${P_out} -f a -N ${bootstrapReps} --bootstop-perms=${bootstrapReps} -n ${transcriptID} -x ${seedSearch} 
  ```
  
  - Coalescent-based species trees with [ASTRAL](https://github.com/smirarab/ASTRAL)
  ```
  ## concatenate all raxml gene trees
  cat ${P_out_raxml}/RAxML_bestTree.${transcriptID}" >> ${all_raxml_inTrees}
  
  ##run ASTRAL
  ${astral_call} -i ${all_raxml_inTrees} -o ${astral_outTree}
  ```

  - Concatenation-based species trees with [iq-tree](http://www.iqtree.org/)
  ```
  ## concatenate all genes
  ## use https://github.com/ballesterus/Utensils/blob/master/geneStitcher.py
  geneStitcher.py -in ${Path_aln_fasta}/*.fasta
  
  ## run iq-tree
  iqtree -s SuperMatrix.fas -spp Partition.txt -bb 1000  -bnni  -m GTR+I+G -nt 16  -mem 50G
   ```
   
  - Divergence time estimation [treePL](https://github.com/blackrim/treePL)
   ```
   ## use config file1 to optimize parameters
   treePL config.1

   ## see output and generate config.2
   treePL config.2 treefile=tree_rooted.dated.config2.tre
   ```
   

#### - Genome-wide Unbiased Selection Screen using [HYPHY-aBSREL](https://stevenweaver.github.io/hyphy-site/methods/selection-methods/)
  For each transcript, the species tree should be trimmed to kepp only branches represented in alignemnt. We used [tree_doctor](https://github.com/UCSantaCruzComputationalGenomicsLab/phast/blob/master/src/util/tree_doctor.c).
   ```
  ## trimm input tree 
  tree_doctor -a -n -P $(grep '>' ${P_out}/${transcriptID}.fa | sed 's/>//g'| cut -f1 | tr '\n' ',') ${tree} > ${P_out}/${trasncriptID}.prunnedTree.tre

  ## run absrel
  hyphy absrel --alignment ${P_out}/${trasncriptID}.fa --tree ${P_out}/${trasncriptID}.prunnedTree.tre --output ${P_out}/${trasncriptID}.ABSREL.json | tee -a ${P_out}/${trasncriptID}.ABSREL.log
   ```

#### - Gene enrichment analyses using [gprofiler2](https://biit.cs.ut.ee/gprofiler/page/r)
   ```R
   ## Load R and library
   library(gprofiler2)
   
   ## run enrichment per mammalian group
   multi_sig_noBG_tmp_gostres = gost(listGenes_groups, organism = "hsapiens", significant = TRUE)
   
   ## save enrichment summary as tab file
   write.table(sapply(multi_sig_noBG_tmp_gostres$result, FUN = paste), file=out_summaryErich, sep="\t", quote=FALSE, row.names=FALSE)
```

#### - Linear regression



