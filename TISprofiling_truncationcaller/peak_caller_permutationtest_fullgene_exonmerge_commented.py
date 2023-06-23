import re
import os
import sys
from Bio import SeqIO
import pysam
import numpy as np

def read_fasta(fastaFile):
    """creates fastaDict from fasta, where key=geneID and value=sequence"""
    fastaDict = {}
    for line in fastaFile:
        if ">" in line:
            gene_id = line.split(">")[1].strip()
        else:
            fastaDict[gene_id]=line.strip()
    return fastaDict  

def read_bed(bedfile):
    """bed file parser that returns a dictionary where key is 
    gene_id, value is [chrom, start, stop, strand, exon_num, exon_len, exon_loc]"""
    
    bed_dict = {}
    for line in bedfile:
        split = line.split("\t")
        
        gene_id = split[3]
        chrom,start,stop,strand = split[0],int(split[1]),int(split[2]),split[5]
        exon_num = int(split[9])
        exon_len = [int(x) for x in split[10].split(',')[:-1]]
        exon_loc = [int(x) for x in split[11].strip().split(',')[:-1]]
        
        exon_list = []
        for i in range(len(exon_loc)):

            exon_genomic_start = start + exon_loc[i]
            exon_genomic_stop = exon_genomic_start + exon_len[i]
                
            exon_info = [exon_genomic_start,exon_genomic_stop]
            exon_list.append(exon_info)
        
        if strand == "-":
            exon_list = exon_list[::-1]
        
        bed_dict[gene_id] = [chrom, start, stop, strand, exon_num, 
                             exon_len, exon_loc, exon_list]
    return bed_dict

                
gene_dict = read_bed(open("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/MvO_v1_annotation_formattedforpeakcaller_noYDL248W_noYML132W.bed"))

#read in bam files for all timepoints
samfile_vegexp=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/veg_exp_LTM_merged_unique.bam","rb")
samfile_vegsat=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/veg_sat_LTM_merged_unique.bam","rb")
samfile_0h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/0h_LTM_merged_unique.bam","rb")
samfile_15h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/1point5h_LTM_merged_unique.bam","rb")
samfile_3h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/3h_LTM_merged_unique.bam","rb")
samfile_45h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/4point5h_LTM_merged_unique.bam","rb")
samfile_6h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/6h_LTM_merged_unique.bam","rb")
samfile_8h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/8h_LTM_merged_unique.bam","rb")
samfile_10h=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/10h_LTM_merged_unique.bam","rb")
samfile_spore=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/22h_LTM_merged_unique.bam","rb")
samfile_mataa=pysam.AlignmentFile("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/MatAA_LTM_merged_unique.bam","rb")


bamlist = [samfile_vegexp,samfile_vegsat,samfile_0h,samfile_15h,samfile_3h,samfile_45h,
           samfile_6h,samfile_8h,samfile_10h,samfile_spore,samfile_mataa]

asites = open("/mnt/brarlab/andrea_higdon/June2016_paradox_copy/asites_LTMAug2016.txt")
countfile = open("20210310_Fullgenes_bytimepoint_forpermutationtest_script.txt","w+")
    
#create dictionary with a-site offset info for each read length; read length is key, offset is value
asite_dict = dict()
for line in asites:
    l = line.strip().split('\t')
    asite_dict[int(l[0])] = int(l[1])

timepoint_list = ["vegexp","vegsat","meiotic_0h","meiotic_1_5h","meiotic_3h","meiotic_4_5h",
                  "meiotic_6h","meiotic_8h","meiotic_10h","spore","mataa"]

for gene_id in gene_dict:
    chrom = gene_dict[gene_id][0]
    start = gene_dict[gene_id][1]
    stop = gene_dict[gene_id][2]
    strand = gene_dict[gene_id][3]
    #length = stop - start
    exon_list = gene_dict[gene_id][7]
    
    
    #create list for tracking read counts at each position
    for bam_index, bam in enumerate(bamlist):
        
        #create list to store counts for all exons
        exon_count_list = []
        
        for exon in exon_list:
            exon_start = exon[0]
            exon_stop = exon[1]
            length = exon_stop - exon_start
            position_list = [1]*length
    
            #identify a-site locations for each read      
            if strand == "+": 
                for read in bam.fetch(chrom,exon_start,exon_stop):
                    if read.flag == 0:
                        if read.alen not in asite_dict:
                            continue
                        asite_position = read.pos + asite_dict[read.alen] - exon_start
                        if asite_position >=0 and asite_position < (length-1):
                            position_list[asite_position] += 1

            if strand == "-":
                for read in bam.fetch(chrom,exon_start,exon_stop):
                    if read.flag == 16:
                        if read.alen not in asite_dict:
                            continue
                        asite_position = exon_stop - (read.pos + read.alen - asite_dict[read.alen])
                        if asite_position >=0 and asite_position < (length-1):
                            position_list[asite_position] += 1

            exon_count_list += position_list
        
        exon_count_list_str = [str(position) for position in exon_count_list]
                               
        printlist = [chrom, str(start), str(stop), gene_id, str(0), strand, timepoint_list[bam_index]] + exon_count_list_str
        
        countfile.write("\t".join(printlist))
        countfile.write("\n")
                    

samfile_vegexp.close()
samfile_vegsat.close()
samfile_0h.close()
samfile_15h.close()
samfile_3h.close()
samfile_45h.close()
samfile_6h.close()
samfile_8h.close()
samfile_10h.close()
samfile_spore.close()
samfile_mataa.close()

asites.close()
countfile.close()