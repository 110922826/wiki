#Analysing ChIP-seq data

## Raw data access

Downloading a raw read archive from the [sequence read archive](https://www.ncbi.nlm.nih.gov/sra).
For example for this experiment: [https://www.ncbi.nlm.nih.gov/sra/?term=SRR529957](https://www.ncbi.nlm.nih.gov/sra/?term=SRR529957), we can download the data with:
```bash
> wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR529/SRR529957/SRR529957.sra
```

## Raw data processing
To extract the acctual sequences and the sequencing scores from the read archive we use:
```bash
> fastq-dump SRR529957.sra
# This will create file SRR529957.fastq
```

## Mapping to a reference genome

Before the mapping procedure we need to find out the read size that we want to map. Here, our reads are all 50bp long single-end reads, thus the choice for the aligners and options.

### Using [bwa](http://bio-bwa.sourceforge.net/)

#### Preparing the reference index
Assuming ```mm10.fa``` contains the complete genome sequence you want to use as a reference from mouse:
```bash
> bwa index mm10.fa
```

#### Mapping to the reference index
Two-command process. First you produce a gapped/ungapped alignment and then you generate the final alignment either single or paired ended.

```bash
> bwa aln /data/EXTERNAL/genomes/bwa-indexes/mm10.fa ../SRR529957.fastq > SRR529957.sai
```

Note: You might want to use the option ```-q```  , which will quality trim the reads before aligning them, e.g. ```-q 20```  (where 20 is the phred-score -> 20 = 99% base call accuracy.)

Generate the final alignment, with single ended reads:

``` bash
> bwa samse /data/EXTERNAL/genomes/bwa-indexes/mm10.fa SRR529957.sai ../SRR529957.fastq > SRR529957.sam
```

or paried-end reads:

```bash
> bwa sampe ref.fa aln_sa1.sai aln_sa2.sai read1.fq read2.fq > aln-pe.sam
```


The .sam format contains the alignments for the reads to the reference.

### Using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 

#### Preparing the reference index
 Assuming ```mm10.fa``` contains the complete genome sequence you want to use as a reference from mouse:
 ```bash
> bowtie2-build mm10.fa mm10
```

We can also download pre-processed index files:
```bash
> wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
> unzip mm10.zip
```

#### Mapping to the reference index
```bash
> bowtie2 -x /data/EXTERNAL/genomes/bowtie2-indexes/mm10 -U ../SRR529957.fastq -S SRR529957.sam
```

## Post-mapping processing
We create a binary form of the alignment file and delete the sam-file using [SAMtools](http://samtools.sourceforge.net/).
```bash
> samtools view -bS SRR529957.sam > SRR529957.bam
> rm SRR529957.sam
```

Now we create a sorted bam-file and delete the original bam-file to save some space.
```bash
> samtools sort SRR529957.bam SRR529957.sorted
> rm SRR529957.bam
```

The sorted bam-file is input into the peak call software.

##Peak calling

### Using [macs2](https://github.com/taoliu/MACS/tree/master/MACS2)

We are calling peaks wich are called with an adjusted p-Value (qvalue) of 0.05.
```bash
> macs2 callpeak -t SRR529957.sorted.bam -f BAM -g mm -n SRR529957.macs2 -B -q 0.05 2> macs2.stderr &
```

This will result in a few files. For example, the SRR529957.macs2_summits.bed file which contains the regions and a score for each peak, which is the ```-log(qvalue)```.

```bash
>head SRR529957.macs2_summits.bed
chr1    4785811 4785812 SRR529957.macs2_peak_1    6.55327
chr1    4858664 4858665 SRR529957.macs2_peak_2    34.49584
chr1    5082988 5082989 SRR529957.macs2_peak_3    10.34553
chr1    6214509 6214510 SRR529957.macs2_peak_4    3.94005
chr1    6406612 6406613 SRR529957.macs2_peak_5    7.54560
chr1    6454425 6454426 SRR529957.macs2_peak_6    12.69790
chr1    6467192 6467193 SRR529957.macs2_peak_7    5.99226
chr1    7088737 7088738 SRR529957.macs2_peak_8    4.01822
chr1    7100720 7100721 SRR529957.macs2_peak_9    13.81626
chr1    7147198 7147199 SRR529957.macs2_peak_10   5.99226
```

This file we can intersect for example using [bedtools](bioinf_bedtools.md) with regions of our interest, to  find some ChIP-seq peaks of the DNA-binding protein under investigation.

 


**_FILE: bioinf_chipseq.md - Sebastian Schmeier - Last update: 2014/05/27_**
