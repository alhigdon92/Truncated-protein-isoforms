library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(vroom)


peak_data <- vroom("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/20210312_peaksums_ATGstarts_peaksumgreaterthanorequalto6.tsv.gz")
data <- read.table("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/20210310_Fullgenes_bytimepoint_forpermutationtest_script.txt",header=FALSE,sep="\t",col.names = paste0("V",seq_len(14740)),fill=TRUE)

colnames(data) = c("chrom","start","stop","tid","score","strand","timepoint",0:14732)

find_pvalue = function(gene_data,n_samples,orf_peak_data) {
    #iterate through all rows of data file (genes x timepoint)
    new_df = tibble()
    for (row in 1:nrow(gene_data)) {
        #identify gene,timepoint
        gene_id = gene_data[row,]$tid
        gene_timepoint = gene_data[row,]$timepoint

        #find empirical distribution for gene,timepoint
        gene_data_cleaned = na.omit(as.numeric((gene_data[row,11:14740])))
        #remove values less than 1
        gene_data_cleaned = gene_data_cleaned[gene_data_cleaned > 1]
        
        if (length(gene_data_cleaned) != 0) {
            sampled_sums = replicate(n_samples, sum(sample(gene_data_cleaned, 3, replace=TRUE)))
            
            #filter peak_data to find all orfs from same gene and timepoint
            gene_timepoint_orfs = orf_peak_data %>% filter(tid==gene_id & timepoint==gene_timepoint)
        
            #for each orf, compare orf peak sum to distribution
            gene_timepoint_orfs_pval = gene_timepoint_orfs %>% rowwise() %>% mutate(pvalue = 1-ecdf(sampled_sums)(peak_sum))
        } else {
            gene_timepoint_orfs_pval = gene_timepoint_orfs %>% mutate(pvalue = NA)
        }

        new_df = bind_rows(new_df,gene_timepoint_orfs_pval)
        
    }
    return(new_df)
}

ATG_starts_pvalue <- find_pvalue(data,10000,peak_data)

vroom_write(ATG_starts_pvalue, "20210312_pvalues_ATGstarts.tsv.gz")
