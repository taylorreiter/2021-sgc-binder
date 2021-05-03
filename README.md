[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/taylorreiter/2021-sgc-binder/HEAD)

# Bin completion with spacegraphcats

Below is a tutorial on bin completion. 
The necessary sequences are in the `data` folder. 
`SRR1976948_1_head.fastq.gz` is a subsampled version of the metagenome `SRR1976948`; it contains the first 100,000 reads of the original sequencing file:
The original name of the metagenome is Hu SB1. 
The microbial community was sampled from an Alaskan Oil reservior by [Hu et al 2016](https://mbio.asm.org/content/7/1/e01669-15). 

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1976948/SRR1976948_1.fastq.gz
zcat SRR1976948_1.fastq.gz | head -n 400000 > SRR1976948_1_head.fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1976948/SRR1976948_2.fastq.gz
zcat SRR1976948_2.fastq.gz | head -n 400000 > SRR1976948_2_head.fastq
```

## Quality control

Remove adapters and Very Low Quality sequences

```
fastp --in1 data/SRR1976948_1_head.fastq.gz \
  --in2 data/SRR1976948_2_head.fastq.gz \
  --out1 SRR1976948_1.trim.fastq.gz \
  --out2 SRR1976948_2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 4 \
  --length_required 31 --correction \
  --json SRR1976948.trim.json \
  --html SRR1976948.trim.html
```

K-mer trim reads

``` 
interleave-reads.py SRR1976948_1.trim.fastq.gz SRR1976948_2.trim.fastq.gz | \
        trim-low-abund.py --gzip -C 3 -Z 18 -M 2e9 -V - -o SRR1976948.abundtrim.fq.gz
```

## Bin completion with spacegraphcats

Download a binned genome. 
This metagenome was analyzed via a *de novo* binning approach by [Hu et al. 2016](https://mbio.asm.org/content/7/1/e01669-15).
The bins were placed in NCBI Genbank.

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/508/995/GCA_001508995.1_ASM150899v1/GCA_001508995.1_ASM150899v1_genomic.fna.gz
```

Create a spacegraphcats configuration file. 

```
catlas_base: 'SRR1976948'
input_sequences:
- SRR1976948.abundtrim.fq.gz
ksize: 31
radius: 1
search:
- data/GCA_001508995.1_ASM150899v1_genomic.fna.gz
searchquick: data/GCA_001508995.1_ASM150899v1_genomic.fna.gz
```

Run spacegraphcats!

```
python -m spacegraphcats run conf1.yml extract_contigs extract_reads --nolock 
```

Spacegraphcats returns many files. 
`SRR1976948_k31_r1_search_oh0` contains the search results:

```
GCA_001508995.1_ASM150899v1_genomic.fna.gz.cdbg_ids.contigs.fa.gz
GCA_001508995.1_ASM150899v1_genomic.fna.gz.cdbg_ids.reads.gz
GCA_001508995.1_ASM150899v1_genomic.fna.gz.cdbg_ids.txt.gz
GCA_001508995.1_ASM150899v1_genomic.fna.gz.contigs.sig
GCA_001508995.1_ASM150899v1_genomic.fna.gz.frontier.txt.gz
GCA_001508995.1_ASM150899v1_genomic.fna.gz.response.txt
results.csv
```

+ `*cdbg_ids.reads.gz` contains the reads or contigs (e.g. fastq or fasta, depending on input format) for the query neighborhood.

## Exploring the output of spacegraphcats

Let's first determine if there are any reads in the query neighborhood for sequences that are not in the original genome bin.
We'll do this by mapping the reads in the query neighborhood against the original genome, and seeing how many reads do not map.

```
bwa index data/GCA_001508995.1_ASM150899v1_genomic.fna.gz
bwa mem -t 1 data/GCA_001508995.1_ASM150899v1_genomic.fna.gz SRR1976948_k31_r1_search_oh0/GCA_001508995.1_ASM150899v1_genomic.fna.gz.cdbg_ids.reads.gz > SRR1976948_sgc_GCA_001508995.1.sam
samtools flagstat SRR1976948_sgc_GCA_001508995.1.sam > SRR1976948_sgc_GCA_001508995.1.flagstat
```

Samtools `flagstat` tells us that 91.3% of reads mapped, meaning 8.07% did not.

```
14033 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
320 + 0 supplementary
0 + 0 duplicates
12900 + 0 mapped (91.93% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Let's extract the unmapped reads to take a peak at what they are. 

```
samtools view -f 4 SRR1976948_sgc_GCA_001508995.1.sam > SRR1976948_sgc_GCA_001508995.1_unmapped.sam
```

Convert the unmapped reads to fastq format

```    
samtools fasta SRR1976948_sgc_GCA_001508995.1_unmapped.sam > SRR1976948_sgc_GCA_001508995.1_unmapped.fastq
``` 

And assemble the unmapped reads with megahit.
*NB* some unmapped reads may originate from contigs that were assembled but did not bin. 
This method would pull out those contigs.
However, for reads that still don't assemble, one may need to try an amino acid assembler or analyze them with read-level tools (e.g. mifaser).

```
megahit -r SRR1976948_sgc_GCA_001508995.1_unmapped.fastq --min-contig-len 142
```

Megahit assembled 115 contigs, all of which are less than 1000 base pairs in length and thus would not bin.
Using [NCBI blastx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome), we see that one of these contigs

```
>k141_62 flag=1 multi=1.0000 len=404
TTATTTTTACAAGGACTGCGGGGGTGGATTTATTTCCACTCCTTTCCTTTTAGGAACGCTGTTAAGCATCTCTTTTAATGTCTAACCTTGCTTTTTTAATACTGGTTCTAATATGCTCTTTAACAATAAGAGTATCTTTACTATTAGTTAGGTTTAAGATTGGCTCCCTGTGATAATCTGTAAACTATCAATTGGTTAAGGCTGACGTTTCCTTCTTTTGTTCTCTCTGTAAGACTACTTTTATGGAGAGATTTCGGAATCCTAATCCGAAGTTGTCCGGAATATTCGTTCGTTATGGGTTCCGGTATAGGGTGGCCTATCTCCAGACAGGTCTATCCTGAAAATTAATTTGCCTTAGGAAGGTTAACTGCTAAGCGCGCTCATGTTCGCCTCTTATTCTTGTA
```

Matches to a hypothetical protein in a different *Desulfotomaculum sp.* genome, as well as to a toxin-antitoxin system HicB family antitoxin in the closely related *Pelotomaculum propionicicum.
This sequence was not originally binned, but represents additional functional content in the metagenome. 

```
hypothetical protein [Desulfotomaculum sp.]
Sequence ID: HAU31430.1
Length: 60
Number of Matches: 1

Alignment statistics for match #1 Score	Expect	Method	Identities	Positives	Gaps	Frame
120 bits(300) 	9e-33 	Compositional matrix adjust. 	59/60(98%) 	60/60(100%) 	0/60(0%) 	-1

Query  326  LEIGHPIPEPITNEYSGQLRIRIPKSLHKSSLTERTKEGNVSLNQLIVYRLSQGANLKPN  147
            +EIGHPIPEPITNEYSGQLRIRIPKSLHKSSLTERTKEGNVSLNQLIVYRLSQGANLKPN
Sbjct  1    MEIGHPIPEPITNEYSGQLRIRIPKSLHKSSLTERTKEGNVSLNQLIVYRLSQGANLKPN  60
```

```
toxin-antitoxin system HicB family antitoxin [Pelotomaculum propionicicum]
Sequence ID: WP_134215993.1
Length: 122
Number of Matches: 1

Alignment statistics for match #1 Score	Expect	Method	Identities	Positives	Gaps	Frame
106 bits(265) 	1e-26 	Compositional matrix adjust. 	53/61(87%) 	56/61(91%) 	1/61(1%) 	-1

Query  332  TCLEIGHPIPEPITNEYSGQLRIRIPKSLHKSSLTERTKEGNVSLNQLIVYRLSQGANLK  153
            TCLEIGH IPEPITNEYSGQLRIRIPKSLH+ SLTER KE N+SLNQLIVY+LSQG NLK
Sbjct  61   TCLEIGHSIPEPITNEYSGQLRIRIPKSLHR-SLTERAKEENISLNQLIVYKLSQGGNLK  119

Query  152  P  150
            P
Sbjct  120  P  120
```
