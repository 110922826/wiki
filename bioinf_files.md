#Commonly used files-formats and conversions

##BED
A file-format to describe positional features in genomic sequences.

**Fields**

_The first three BED fields are required,  9 additional optional BED fields are possible._

1. **chrom** - The name of the chromosome (e.g. chr3, chrY)
2. **chromStart** - The starting position of the feature in the chromosome or scaffold. **_The first base in a chromosome is numbered 0_**.
3. **chromEnd** - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
4. **name** - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
5. **score** - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray.
6. **strand** - Defines the strand - either '+' or '-'.
7. **thickStart** - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
8. **thickEnd** - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
9. **itemRgb** - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
10. **blockCount** - The number of blocks (exons) in the BED line.
11. **blockSizes** - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
12. **blockStarts** - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 

**Example**
```bash
chr19   6159262 6159269 ESR2_si 6.47    -
chr19   6159262 6159272 THA_f1  6.53    -
chr19   6159262 6159275 RORG_f1 8.4     -
chr19   6159262 6159279 THA_f2  9.66    -
```

Note! Often I only use the first 6 fields and attach my own fields after this. This allows me to use tools like [bedtools](tutorial_bedtools.md) on my files. I also tend to gzip all bedfiles to save space.

**WWW**
* [http://genome.ucsc.edu/FAQ/FAQformat](http://genome.ucsc.edu/FAQ/FAQformat)


##GFF
A file-format to describe positional features in genomic sequences.

**Fields**
1. **seqname** - The name of the sequence. Must be a chromosome or scaffold.
2. **source** - The program that generated this feature.
3. **feature** - The name of this type of feature. Some examples of standard feature types are "CDS", "start-codon", "stop-codon", and "exon".
4. **start** - The starting position of the feature in the sequence. **_The first base is numbered 1_**.
5. **end** - The ending position of the feature (inclusive).
6. **score** - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
7. **strand** - Valid entries include '+', '-', or '.' (for don't know/don't care).
8. **frame** - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
9. **group** - All lines with the same group are linked together into a single item.

**_Example_**
```bash
browser position chr22:10000000-10025000
browser hide all
track name=regulatory description="TeleGene(tm) Regulatory Regions"
visibility=2
chr22  TeleGene enhancer  10000000  10001000  500 +  .  touch1
chr22  TeleGene promoter  10010000  10010100  900 +  .  touch1
chr22  TeleGene promoter  10020000  10025000  800 -  .  touch2
```

**WWW**
* [http://www.genome.ucsc.edu/FAQ/FAQformat.html#format3](http://www.genome.ucsc.edu/FAQ/FAQformat.html#format3)

##SAM
A SAM file (.sam) is a tab-delimited text file that contains sequence alignment data.

**WWW**
* [http://samtools.sourceforge.net/SAM1.pdf](http://samtools.sourceforge.net/SAM1.pdf)

##BAM
A BAM file (.bam) is the binary version of a SAM file. 

##WIG
**WWW**
* [http://www.genome.ucsc.edu/goldenPath/help/wiggle.html](http://www.genome.ucsc.edu/goldenPath/help/wiggle.html)

##SRA
**WWW**
* [https://www.ncbi.nlm.nih.gov/books/NBK47537/](https://www.ncbi.nlm.nih.gov/books/NBK47537/)

##FASTA
```
>Batf
attgggccactgtctggatgttaggactagatcccaagtcctctggacgagcagccagtgttcctaactgctgagcatctttccagtcccCAAACTTTAAAAAAAAATAGCCAAAACACTGGGATGAGAATGGGACTGTTCACTAAAGCGACAAGCGTGTTCGTGTCATGGGACACATGCCAAACACCGCCACTCCACAGCATGCTTTACAGAGTGGGGGACATTCCTATTCAAATTGTCATCCCTCATTGCAGTAGCCCTGAAAGGTAGTAGGTGGGGCAGGTTTTAGGTGTGTGTGTGTGGTGGGAGGCTGGGGGGGGGGGGAGTAGGGATAGGGAGTCACTGGATTCTCATCATATCGATCTGACTCCAAGTCCAGAGCCTCTGTCCCTCATTTTGACGTTCTCTCCTTGTTGGGTTATGACGTCACTGACTTCTGAAGACTCAAGGTCAACCTTGAACTCCTGAAGGAGCCCCTCCACAGAAGGTAGGTATCCTTTCTCTAGGGAGAATGAAAGTGTCTTGCTAAGAGGACAGGGGATAGATAATAAACTTCTTTAGCATCGA ```

Note! A file can hold more then one sequence.

##FASTQ
Very similar to fasta-files. Contains qualtiy scores for the sequence in ASCII-format. These files are often used for storing sequencing results.

```
@EAS54_6_R1_2_1_413_324
CCCTTCTTGTCTTCAGCGTTTCTCC
+
;;3;;;;;;;;;;;;7;;;;;;;88
```

**WWW**
* [http://maq.sourceforge.net/fastq.shtml](http://maq.sourceforge.net/fastq.shtml)

##MOTIF: TRANSFAC
```
 //
AC  Arnt
XX
ID  MA0004.1
XX
NA  Arnt
XX
P0      A       C       G       T
01      4.000000        16.000000       0.000000        0.000000        C
02      19.000000       0.000000        1.000000        0.000000        A
03      0.000000        20.000000       0.000000        0.000000        C
04      0.000000        0.000000        20.000000       0.000000        G
05      0.000000        0.000000        0.000000        20.000000       T
06      0.000000        0.000000        20.000000       0.000000        G
XX
//
```
## MOTIF: JASPAR
```
>MA0004.1 Arnt
A  [ 4 19  0  0  0  0 ]
C  [16  0 20  0  0  0 ]
G  [ 0  1  0 20  0 20 ]
T  [ 0  0  0  0 20  0 ]
```

## MOTIF: MEME
``` 
MOTIF MA0004.1
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0.0e+00
 0.200000  0.800000  0.000000  0.000000
 0.950000  0.000000  0.050000  0.000000
 0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000
 0.000000  0.000000  0.000000  1.000000
 0.000000  0.000000  1.000000  0.000000
```

##MOTIF: fasta-like
```
>MA0004.1 Arnt
4       16      0       0
19      0       1       0
0       20      0       0
0       0       20      0
0       0       0       20
0       0       20      0
```

##File conversions

###GFF to BED
```bash
> zcat file.gff.gz | egrep '^chr' | awk '{OFS="\t"; print $1,$4-1,$5,$1":"$4-1"-"$5"_"$3,$6,$7}'  | gzip > file.bed.gz
```

###BAM to BED
```bash
> bedtools bamtobed -i file.bam | gzip > file.bed.gz
```

###BED to FASTA
```bash
# needs a genome to extract the sequence form, e.g. human genome fast-file hg19.fa
# strand-specific
> zcat file.bed.gz | bedtools getfasta -s -fi hg19.fa -bed stdin -fo stdout | gzip > file.fa.gz

# not strand-specific, only plus strand
> zcat file.bed.gz | bedtools getfasta -fi hg19.fa -bed stdin -fo stdout | gzip > file.fa.gz
```

###SRA to BED
```bash
# This will create a bed file, but only from aligned reads:
> sam-dump public/sra/SRR770745.sra | samtools view -Sb - | bedtools bamtobed -i stdin | gzip > SRR770745.bed.gz 
```

**_FILE: bioinf_files.md - Sebastian Schmeier - Last update: 2013/12/12_**
