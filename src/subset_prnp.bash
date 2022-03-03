#!/bin/bash

gzcat GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | tail -n +3 | head -1 > data/gtex/prnp_tpm.txt
gzcat GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | grep PRNP >> data/gtex/prnp_tpm.txt
src/transpose.bash data/gtex/prnp_tpm.txt > data/gtex/prnp_tpm_t.txt
