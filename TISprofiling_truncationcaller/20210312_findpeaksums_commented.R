library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(vroom)


#import data 
data = read.table("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/20210310_Fullgenes_bytimepoint_forpermutationtest_script.txt",header=FALSE,sep="\t",col.names = paste0("V",seq_len(14740)),fill=TRUE)

#filter to get ATG starts only
starts = fread("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/20210321_inframe_starts.txt",header=FALSE,sep="\t")
starts = starts[,1:11]
colnames(starts) = c("chrom","start","stop","orfname","score","strand","tid","codon","orftype","orflength","dist_from_annotated")
colnames(data) = c("chrom","start","stop","tid","score","strand","timepoint",0:14732)

ATG_starts = starts %>% filter(codon=="ATG")

#function find peak sums
find_peak = function(orf_df,count_df) {
    new_df = tibble()
    for (row in 1:nrow(orf_df)) {
        orf_tid = orf_df[row,]$tid
        index = orf_df[row,]$dist_from_annotated*3
        count_data = count_df %>% filter(tid==orf_tid)
        
        for (timepoint in 1:nrow(count_data)) {
            peak_data = count_data[timepoint,c(7,(index+8):(index+10))]
            peak_data = peak_data %>% rename("0" = as.character(index),"1" = as.character(index+1),"2" = as.character(index+2))
            orf_data = bind_cols(orf_df[row,],peak_data)
            orf_data = orf_data %>% mutate(peak_sum = `0`+`1`+`2`)

            new_df = bind_rows(new_df,orf_data)
        }
    }
    return(new_df)
}

#call peak sums 
peak_data <- find_peak(ATG_starts,data)

vroom_write(peak_data, "20210312_peaksums_ATGstarts.tsv.gz")
