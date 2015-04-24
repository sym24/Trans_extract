Sequence extract for transciptome sequence
usage: trans_extract.py [-h] -g file -t gtf_file -r ref_genome [-o STR]

Extract target genes transcriptome and generate bloom filter for each gene

optional arguments:
  -h, --help            show this help message and exit
  -g file, --genes file
                        Absolute path of text file containing target gene
                        names. [Required]
  -t gtf_file, --gtf gtf_file
                        Absolute path of gtf file.[Required]
  -r ref_genome, --ref ref_genome
                        Absolute path of the reference sequence file.
                        [Requried]
  -o STR, --outdir STR  Output file path. Default:
                        [/projects/btl/ymingsun/seqExtract/single_exon]

