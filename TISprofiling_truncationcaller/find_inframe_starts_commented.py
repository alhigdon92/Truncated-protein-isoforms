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

        bed_dict[gene_id] = [chrom, start, stop, strand, exon_num, exon_len, exon_loc]
    return bed_dict
    
def find_all_starts(fasta_dict, bed_dict):
    """Find all in-frame start codons within sequence, return dictionary 
    where key is gene_id, value is [position, codon]"""
    
    start_codons = ["ATG","CTG","GTG","TTG","AAG","ACG","AGG","ATC","ATT","ATA"]
    
    #create dictionary with keys for every gene_id
    start_dict = {}
    
    for gene_id in fasta_dict:
        start_dict[gene_id] = []
        
        #get gene sequence and strand information from fasta
        seq = fasta_dict[gene_id]
        strand = bed_dict[gene_id][3]
        
        #search entire gene sequence for in-frame codons in list of start codons, 
        #identify as either annotated start or possible truncation start
        for pos in range(0, len(seq)-2, 3):
            codon = seq[pos:pos+3]
            if codon in start_codons:
                if pos == 0:
                    orftype = "annotated"
                else:
                    orftype = "truncation"
                
                #create unique ID for every in-frame start codon and add to start_dict
                if strand == "+":
                    orflength = int((len(seq)-pos)/3)
                    unique_id = gene_id+'_'+str(orflength)+'aa_'+codon
                    dist_from_annotated = int(pos/3)
                    start_dict[gene_id].append([pos, codon, orftype, orflength, dist_from_annotated, unique_id])
                elif strand == "-":
                    rev_pos = len(seq)-pos-1
                    orflength = int((len(seq)-pos)/3)
                    dist_from_annotated = int(pos/3)
                    unique_id = gene_id+'_'+str(orflength)+'aa_'+codon
                    start_dict[gene_id].append([rev_pos, codon, orftype, orflength, dist_from_annotated, unique_id])
    
    return start_dict

def assign_coordinates(start_dict,bed_dict):
    """Use information in start_dict to assign genomic coordinates to each in-frame start, 
    incorporating splicing information"""
    
    #create dictionary to store coordinates for each in-frame start
    coord_dict = {}
    
    for gene_id in start_dict:
        
        coord_dict[gene_id] = []
        
        strand = bed_dict[gene_id][3]
        exon_num = bed_dict[gene_id][4]
        exon_loc = bed_dict[gene_id][6]
        exon_len = bed_dict[gene_id][5]
        chrom = bed_dict[gene_id][0]
        start = bed_dict[gene_id][1]
        stop = bed_dict[gene_id][2]
        
        exon_list = []
        
        #calculate genomic coordinate of start codons and add to coord_dict
        for i in range(len(exon_loc)):

            exon_genomic_coord = start + exon_loc[i]

            if i==0:
                exon_transcript_start = exon_loc[i]
                exon_transcript_stop = exon_len[i] - 1
            else:
                exon_transcript_start = exon_len[i-1] 
                exon_transcript_stop = exon_transcript_start + exon_len[i] - 1
                
            exon_info = [exon_transcript_start, exon_transcript_stop, exon_genomic_coord, exon_len, exon_loc]
            exon_list.append(exon_info)
         
        for start_codon in start_dict[gene_id]:
            pos = start_codon[0]
            codon = start_codon[1]
            orftype = start_codon[2]
            orflength = start_codon[3]
            dist_from_annotated = start_codon[4]
            unique_id = start_codon[5]
            
            for exon in exon_list:
                exon_start = exon[0]
                exon_stop = exon[1]
                genomic_coord = exon[2]
                
                if pos >= exon_start and pos <= exon_stop:
                    codon_genomic_coord = genomic_coord + (pos-exon_start+1)
                    if strand == "+":
                        coord_dict[gene_id].append([chrom, str(codon_genomic_coord-1), str(stop), 
                                                    unique_id, strand, "0", gene_id, codon, orftype, 
                                                    str(orflength), str(dist_from_annotated)])
                    elif strand == "-":  
                        coord_dict[gene_id].append([chrom, str(start), str(codon_genomic_coord), 
                                                    unique_id, strand, "0", gene_id, codon, orftype, 
                                                    str(orflength), str(dist_from_annotated)])
                else:
                    pass
 
    return coord_dict
                
gene_dict = read_bed(open("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/MvO_v1_annotation_formattedforpeakcaller_noYDL248W_noYML132W.bed"))

fasta_dict = read_fasta(open("/mnt/brarlab/andrea_higdon/20210201_LTMandCHX_relignedtoKeeneyLabGenome/MvO_v1_annotation_formattedforpeakcaller_annotatedgenes_spliced.fasta"))

start_dict = find_all_starts(fasta_dict,gene_dict)

coord_dict = assign_coordinates(start_dict,gene_dict)

printfile = open("20210321_inframe_starts.txt","w+")
for gene in coord_dict:
    for orf in coord_dict[gene]:
        orflist = orf + ["\n"]
        printfile.write("\t".join(orflist))
