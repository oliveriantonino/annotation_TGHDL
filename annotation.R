library(data.table)
library(magrittr)
library(stringr)
library(tidyverse)

filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename
slice <- dplyr::slice
count <- dplyr::count


# Input and output --------------------------------------------------------

input_independentSNPs_filename="independent_SNPs.txt"

input_depict_gene_priorization_filename="DEPICT_geneprioritization.txt"

input_exLD_filename="res_exonic_variants_in_LD.txt"
input_independent_synonymous_SNPs_filename <- "ID_synonymous_SNPs.txt"

input_res_gtex_filename = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

input_res_starnet_filename = "res_starnet.txt"

output_folder = "/nfs/corenfs/INTMED-Speliotes-data/Projects/UK_ATLAS/IndivProj/Antonino/annotation_TGHDL/output"

# Start of the script -----------------------------------------------------

if(!dir.exists(output_folder)){dir.create(output_folder)}

res_starnet = 
  as_tibble(fread(paste0("input/",input_res_starnet_filename), na=".")) %>%
  mutate(Description = ifelse(is.na(`Gene symbol`), `Ensembl gene ID`, `Gene symbol`)) %>%
  select(-starts_with("p_"), -`Ensembl gene ID`, -`Gene symbol`) %>%
  select(rsID, Description, everything()) %>%
  mutate(across(.cols = where(is.numeric), .fns=function(v){ifelse(v<0.05, TRUE, FALSE)}))

colnames(res_starnet) = gsub("adj_", "", colnames(res_starnet))
res_starnet[is.na(res_starnet)] <- FALSE

res_gtex = 
  as_tibble(fread(paste0("input/",input_res_gtex_filename))) %>%
  mutate(
    Description = ifelse(Description=="-", Name, Description)
  )

input_independent_synonymous_SNPs <- scan(paste0("input/",input_independent_synonymous_SNPs_filename), sep = "\n", what=c("character"))

res_exLD =
  as_tibble(fread(paste0("input/",input_exLD_filename))) %>%
  select(SNPID, starts_with("ex_"), change, r2=R2) 

#independentSNPs
res_independentSNPs =
  as_tibble(
    fread(
      paste0("input/",input_independentSNPs_filename)
    )
  ) %>%
  mutate(
    location = ifelse(SNPID %in% input_independent_synonymous_SNPs_filename, "synonymous_SNV", location)
    ) %>%
  select(
    chr=CHR, pos.hg19= POS, pos.hg38=POS.b38, rsID=SNPID,
    OA, EA, EAF=freqEA, N,
    location, nearestGene
  ) %>%
  left_join(res_exLD, by=c("rsID"="SNPID")) 

# DEPICT

gene_priorization = 
  as_tibble(fread(paste0("input/",input_depict_gene_priorization_filename))) %>%
  mutate(
    `Gene symbol` = ifelse(`Gene symbol`=="-", `Ensembl gene ID`, `Gene symbol`)
  )

locus.split= 
  gene_priorization %>% 
  .$Locus %>%
  str_split(pattern = "[;]")

gpNew <- tibble(rsID = unlist(locus.split),
                rsID.orig = rep(gene_priorization$Locus, lengths(locus.split))) %>%
  filter(rsID %in% res_independentSNPs$rsID) %>%
  unique()

temp_gene_priorization = 
  gpNew %>%
  left_join(gene_priorization, by=c("rsID.orig"="Locus"), relationship = "many-to-many") %>%
  select(rsID, rsID.orig, DEPICTprioritizedGene=`Gene symbol`,
         FDR=`False discovery rate`,
         NgenesInLocus=`Nr of genes in locus`,
         DEPICTNominalPValue=`Nominal P value`) %>%
  filter(FDR %in% c('<=0.01','<0.05'))

temp_gene_priorization_list = 
  temp_gene_priorization %>%
  arrange(rsID, DEPICTNominalPValue) %>%
  group_by(rsID) %>%
  mutate(
    DEPICTprioritizedGenes = paste0(DEPICTprioritizedGene, collapse = "|"),
    DEPICTFDRs = paste0(FDR, collapse = "|"),
    DEPICTNominalPValues = paste0(DEPICTNominalPValue, collapse = "|")
  ) %>%
  slice_head(n=1) %>%
  select(-rsID.orig, -DEPICTprioritizedGene, -FDR, -NgenesInLocus, - DEPICTNominalPValue) %>%
  ungroup()

res_independentSNPs_depict =
  res_independentSNPs %>%
  left_join(temp_gene_priorization_list, by=c("rsID"))

# STARNET

temp_genes = unlist(str_split(c(res_independentSNPs_depict$nearestGene, res_independentSNPs_depict$ex_gene, res_independentSNPs_depict$DEPICTprioritizedGenes), "[,;|]"))
temp_genes = temp_genes[!is.na(temp_genes)]

filt_res_starnet =
  res_starnet %>%
  filter(Description %in% temp_genes, rsID %in% res_independentSNPs_depict$rsID) %>%
  arrange(Description, rsID) %>%
  set_colnames(c("rsID", "gene", "QAW", "QB", "QIMA", "QL", "QMS", "QAS", "QAV"))

filt_res_starnet %<>%
  select(rsID, gene, QAS, QMS, QL, QAV, everything())

# GTEX

res_gtex_median =
  bind_cols(res_gtex[,c(1:2)], median=apply(res_gtex[,-c(1:2)],1,median))  

res_gtex_expr = 
  bind_cols(res_gtex[,c(1:2)], res_gtex[,-c(1:2)] > 2*res_gtex_median$median)

gene_list = list(
  AS = res_gtex_expr$Description [ res_gtex_expr$`Adipose - Subcutaneous`],
  P = res_gtex_expr$Description [ res_gtex_expr$Pancreas],
  AG = res_gtex_expr$Description [ res_gtex_expr$`Adrenal Gland`],
  MS = res_gtex_expr$Description [ res_gtex_expr$`Muscle - Skeletal`],
  U = res_gtex_expr$Description [ res_gtex_expr$Uterus],
  L = res_gtex_expr$Description [ res_gtex_expr$Liver],
  AV = res_gtex_expr$Description [ res_gtex_expr$`Adrenal Gland`]
)

all_genes = unique(c(res_independentSNPs_depict$nearestGene, res_independentSNPs_depict$ex_gene, res_independentSNPs_depict$DEPICTprioritizedGenes))
all_genes = sort(unique(unlist(str_split(all_genes, "[,;|]"))))

gene_expr_tissue = 
  tibble(
    gene = all_genes,
    EAS = logical(length(all_genes)),
    EP = logical(length(all_genes)),
    EAG = logical(length(all_genes)),
    EMS = logical(length(all_genes)),
    EU = logical(length(all_genes)),
    EL = logical(length(all_genes)),
    EAV = logical(length(all_genes))
  ) %>%
  mutate(
    EAS = ifelse(gene %in% gene_list$AS, TRUE, FALSE),
    EP = ifelse(gene %in% gene_list$P, TRUE, FALSE),
    EAG = ifelse(gene %in% gene_list$AG, TRUE, FALSE),
    EMS = ifelse(gene %in% gene_list$AS, TRUE, FALSE),
    EU = ifelse(gene %in% gene_list$U, TRUE, FALSE),
    EL = ifelse(gene %in% gene_list$L, TRUE, FALSE),
    EAV = ifelse(gene %in% gene_list$AV, TRUE, FALSE),
  )

temp_letters = c("N","X", "X*","D")

res_independentSNPs_depict$gene_annotation = character(nrow(res_independentSNPs_depict))
res_independentSNPs_depict$Locus = character(nrow(res_independentSNPs_depict))

for(i in seq.int(nrow(res_independentSNPs_depict))){
  temp_genes = unlist(str_split(c(res_independentSNPs_depict$nearestGene[i], res_independentSNPs_depict$ex_gene[i], res_independentSNPs_depict$DEPICTprioritizedGenes[i]), "[,;|]"))
  temp_genes = temp_genes[!is.na(temp_genes)]
  temp_tab = table(temp_genes)
  temp_tab = temp_tab[order(temp_tab, decreasing = T)]
  temp_genes = names(temp_tab)
  max_temp_tab = temp_tab[which(temp_tab==max(temp_tab))]
  
  if(!is.na(res_independentSNPs_depict$ex_gene[i])){
    
    res_independentSNPs_depict$Locus[i] = res_independentSNPs_depict$ex_gene[i]
    
  }else{
    
    if(max(max_temp_tab)==2 & length(max_temp_tab)==1){
      
      res_independentSNPs_depict$Locus[i] = names(max_temp_tab)
      
    }else{
      
      if(max(max_temp_tab)==2 & length(max_temp_tab)>1){
        
        res_independentSNPs_depict$Locus[i] = paste(names(max_temp_tab), collapse=",")
        
      }else{
        
        temp_gene_tissue_expr =
          gene_expr_tissue %>%
          filter(gene %in% names(max_temp_tab)) %>%
          select(-gene) %>%
          rowSums(na.rm=T)
        
        names(temp_gene_tissue_expr) = names(max_temp_tab)
        
        temp_starnet =
          filt_res_starnet %>%
          filter(gene %in% names(max_temp_tab),
                 rsID %in% res_independentSNPs_depict$rsID[i])
        
        temp_starnet_expr =
          temp_starnet %>%
          select(-rsID, -gene) %>%
          rowSums(na.rm = T)
        
        names(temp_starnet_expr) = temp_starnet$gene
        
        genes_not_in_temp_starnet_expr = names(max_temp_tab)[!(names(max_temp_tab) %in% names(temp_starnet_expr))]
        
        if(length(genes_not_in_temp_starnet_expr)>0){
          fake_sum_of_genes_not_in_temp_starnet_expr = integer(length = length(genes_not_in_temp_starnet_expr))
          names(fake_sum_of_genes_not_in_temp_starnet_expr) = genes_not_in_temp_starnet_expr
          temp_starnet_expr = c(temp_starnet_expr, fake_sum_of_genes_not_in_temp_starnet_expr)
        }
        
        temp_starnet_expr = temp_starnet_expr[match(names(temp_gene_tissue_expr), names(temp_starnet_expr))]
        
        temp_gene_tissue_expr = temp_gene_tissue_expr + temp_starnet_expr
        
        max_temp_gene_tissue_expr = temp_gene_tissue_expr[which(temp_gene_tissue_expr==max(temp_gene_tissue_expr))]
        
        res_independentSNPs_depict$Locus[i] = paste(names(max_temp_gene_tissue_expr), collapse=",")
        
      }
      
    }
    
  }
  
  temp_label = character(length(temp_genes))
  for(j in seq.int(temp_label)){
    
    label_fp =
      temp_letters[
        c(ifelse(sum(temp_genes[j] %in% unlist(str_split(c(res_independentSNPs_depict$nearestGene[i]), "[,;|]")))>0, TRUE, FALSE),
          ifelse((sum(temp_genes[j] %in% unlist(str_split(c(res_independentSNPs_depict$ex_gene[i]), "[,;|]")))>0) & res_independentSNPs_depict$location[i]=="exonic", TRUE, FALSE),
          ifelse((sum(temp_genes[j] %in% unlist(str_split(c(res_independentSNPs_depict$ex_gene[i]), "[,;|]")))>0) & res_independentSNPs_depict$location[i]!="exonic", TRUE, FALSE),
          ifelse(sum(temp_genes[j] %in% unlist(str_split(c(res_independentSNPs_depict$DEPICTprioritizedGenes[i]), "[,;|]")))>0, TRUE, FALSE)
      )]
    
    label_fp = label_fp[!is.na(label_fp)]
    
    label_sp = 
      gene_expr_tissue %>% 
      filter(gene==temp_genes[j]) %>%
      select(-gene) %>%
      unlist()
    
    label_sp =  names(label_sp)[label_sp]
    
    label_tp = 
      filt_res_starnet %>% 
      filter(gene==temp_genes[j], rsID==res_independentSNPs_depict$rsID[i]) %>%
      select(-gene, -rsID) %>%
      unlist()
    
    label_tp =  names(label_tp)[label_tp]
    
    temp_label[j] =
      paste0(
        temp_genes[j],
        "(",
        
        paste(
          c(
            
            label_fp,
            
            label_sp,
            
            label_tp
          ),
          collapse = ", "),
        
        ")"
      ) 
    
  }
  
  res_independentSNPs_depict$gene_annotation[i] = paste(temp_label, collapse = "; ")
  
}

res_independentSNPs_depict %<>% 
  mutate(
    location=ifelse(location %in% c("exonic"), "nonsynonymous_SNV", location)
  )

red_res_independentSNPs_depict = 
  res_independentSNPs_depict %>%
  select(chr, pos.hg19, pos.hg38, OA, EA, EAF, location, rsID, Locus, gene_annotation)


red_res_independentSNPs_depict %>%
  arrange(chr, pos.hg19) %>%
  fwrite(
    paste0(output_folder,"/annotation.txt"),
    sep="\t",
    na = "."
  )

