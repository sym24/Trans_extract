#!/usr/bin/env python

# improvement can be done
# 1 .Use transcript ID as sequence header
# >[Transcript ID]
# 2. Combine genome and transcriptome sequence in one script

# [modified on Jan27]
# [Add funtion: ]
# Attract transcript_ids for each gene and
# stored in seperate file for each gene
# [Delete funtion: ]
# Delete check output exists option.
# Disabled buffer
# Disables prefix. It is default as gene name
import os
import sys
from pybedtools import BedTool
import pysam
import subprocess
import argparse
from multiprocessing import Pool

def get_gene_name(target_file):
    """Generate set contains gene name without duplication."""
    with open(target_file, 'r') as fi:
        gene_set = set(line.rstrip() for line in fi)
    print 'Total number of genes: ', len(gene_set)
    return (gene_set, len(gene_set))

#extract exon startp endp with same transcrit_id for each gene
def get_bounds(gene_set, gtf_file):
    """Go through gtf file extract the longest sequence for genes in gene_set"""
    print '...Extract gene boundaries...'
    bounds = {}
    gene2id = {}
    id2coord = {}
    for feature in BedTool(gtf_file):
        if feature[2] == 'exon' and feature.attrs.has_key('gene_id'):
             gene = feature.attrs['gene_id']
             if gene in gene_set:
                 trans_id = feature.attrs['transcript_id']
                 startp = int(feature.start)
                 endp = int(feature.stop)
                 if gene in gene2id:
                     if trans_id in gene2id[gene]:
                         id2coord[trans_id].append([feature[0], startp, endp, \
                                                    feature.attrs['exon_number']])
                     else:
                         gene2id[gene].append(trans_id)
                         id2coord[trans_id] = [[feature[0], startp, endp, \
                                                feature.attrs['exon_number']]]
                 else:
                     gene2id[gene] = [trans_id]
                     id2coord[trans_id] = [[feature[0], startp, endp, \
                                            feature.attrs['exon_number']]]
    find_num = 0
    miss_genes = []
    for gene in gene_set:
        if gene in gene2id:
            find_num += 1
        else:
            miss_genes.append(gene)
    print 'Number of genes successfully captured: %d'%(find_num)
    if miss_genes:                
        print 'Number of genes failed to captured: %d'%len(miss_genes)
        print 'Missing genes: %s'%', '.join(miss_genes)
    return gene2id, id2coord

def extract_seq(outdir, ref_file, gene2id, id2coord):
    """Extract genes from reference genome and write it to output file."""
    print '...Extract gene sequence...'
    outdic = {}
    for gene, id_list in gene2id.iteritems():
        outfile = os.path.join(outdir, gene + '.fa')
        outdic[gene] = outfile
        with open(outfile, 'w') as fi:    
            ref_fasta = pysam.Fastafile(ref_file)
            for id in id_list:
                last = len(id2coord[id]) - 1
                fi.write("%s%s %s %s:%d-%d\n" % ('>',gene, id, id2coord[id][0][0], \
                                                 id2coord[id][0][1], \
                                                 id2coord[id][last][2]))
                # sample header: >ZCCHC8 NM_017612 chr12:122956145-122985620
                for coord in id2coord[id]:
                    seq = ref_fasta.fetch(coord[0], coord[1], coord[2])
                    fi.write("%s" % seq)
                fi.write('\n')
    return outdic

def make_filter((infile, gene, outdir)):
    """Make bloom filter out of targeted genes"""
    try:
        subprocess.check_output(['samtools', 'faidx', infile])
        cmd_pass = ['biobloommaker', '-p', gene, '-o', outdir, infile]
        subprocess.check_output(cmd_pass)
    except subprocess.CalledProcessError as e:
        sys.exit(e)    

def cmd_pass():
    parser = argparse.ArgumentParser(description='Extract target genes transcriptome and '
                                     'generate bloom filter for each gene')
    # Required arguments'
    parser.add_argument('-g', '--genes', dest='genes', metavar='file',
                        type=argparse.FileType('r'), required=True,
                        help='Absolute path of text file containing '
                        'target gene names. [Required]')
    parser.add_argument('-t', '--gtf', dest='gtf', metavar='gtf_file',
                        type=file, required=True,
                        help='Absolute path of gtf file.[Required]')
    parser.add_argument('-r', '--ref', dest='ref', metavar='ref_genome',
                        type=file, required=True,
                        help='Absolute path of the reference sequence file. [Requried]')

    # Optional variables
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='STR', type=str,
                        help='Output file path. Default: [%(default)s]',
                        default=os.getcwd())
    return parser

def main():
    parser = cmd_pass()
    args = parser.parse_args()

    # Make output file name
    target_fa = os.path.join(args.outdir, 'fasta_files', '<gene name>'+'.fa')
    dir_fa = os.path.join(args.outdir, 'fasta_files')
    dir_bf = os.path.join(args.outdir, 'bf_files')    
    # Make output directory if not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if not os.path.exists(dir_fa):
        os.makedirs(dir_fa)
    if not os.path.exists(dir_bf):
        os.makedirs(dir_bf)
    print 'Genome file: ', os.path.abspath(args.ref.name)
    print 'Gene list: ', os.path.abspath(args.genes.name)
    print 'GTF file: ', os.path.abspath(args.gtf.name)
    print 'Output FASTA file: ', target_fa
    # Extract gene name from input file 
    gene_set, gene_num = get_gene_name(args.genes.name)
    # Get gene coordinates from GTF file  
    gene2id, id2coord = get_bounds(gene_set, args.gtf.name)
    # Get gene sequences from reference sequence file
    indict = extract_seq(dir_fa, args.ref.name, gene2id, id2coord)
    input = [[indict[gene], gene, dir_bf] for gene in gene2id]
    p = Pool(5)
    p.map(make_filter, input)

if __name__ == '__main__':
    main()
