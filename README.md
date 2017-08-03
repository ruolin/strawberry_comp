Simulation and evaluation framework used in Strawberry
==================================
This project contains a step-by-step tutorial for simulating RNA-Seq samples and benchmarking popular bioinformatics programs for genome-guied assembly and quantifcation using the simulated data. This tutorial will generate the result that is published in Strawberry paper. 

Human RNA-Seq data
=====================================

Prerequisites 
======================
* Human genome annotation.

  I currently use the genecode GRCh37.p13 version (https://www.gencodegenes.org/releases/19.html)
* Human genome.

  I currently use the GRCh37/hg19 genome. 
* GEUVADIS RNA-Seq raw reads data (http://www.geuvadis.org/web/geuvadis/rnaseq-project).

     I use 6 samples from this data set, which are: 
    
    `ERR188021, ERR188052, ERR188088, ERR188204, ERR188317, ERR188453`
* A RAN-Seq splice alinger, such as HISAT 2 (http://www.ccb.jhu.edu/software/hisat/index.shtml).
* Samtools (http://samtools.sourceforge.net/)
* Polyester simulator (https://bioconductor.org/packages/release/bioc/html/polyester.html)
* Programs that are being evaluted in the Strawberry's paper.
    * Strawberry (https://github.com/ruolin/strawberry/releases/tag/0.9.0)
    * Cufflinks (http://cole-trapnell-lab.github.io/cufflinks/)
    * StringTie (https://ccb.jhu.edu/software/stringtie/)
    
    After installing these programs, add the path to the executables to your current PATH environment variable. 
    In linux and a bash shell, this can be done by type 
    
    `export PATH=/path/to/your/intallation/bin/:$PATH`
    
Workflow
=====================================

1. git clone https://github.com/ruolin/strawberry_comp.git
2. Aligns the 6 samples using HISAT 2 or any other splice aligner. Then convert the results to BAM and sort them using Samtools.
3. Use Cufflinks to quantify known transcript abundances of the six samples.
    For example, assume you have the alignment file (BAM) and gene annotation file (gff3) under your current directory, then you do 
    
    `cufflinks -G /path/to/gencode.v19.annotation.gff3 ERR188021.sorted.bam -o ERR188021_cuff.count -p 8`.
    
    Cufflinks will generate a `transcripts.gtf` file. Rename this file with a meaningful name such as `ERR188021.gtf`. 
4. Generate the transcript expression matrix for Polyester.
  
    After you got all the gtf files produced by cufflinks, the following script can be used to generate a transcript FPKM matrix. This matirx is a row-major matrix, where each row represents a transcript and each columns represents a sample. 
  
    ` scripts/gff-print.py -f ERR188021_cuff.count/ERR188021.gtf ERR188052_cuff.count/ERR188052.gtf ERR188088_cuff.count/ERR188088.gtfERR188204_cuff.count/ERR188204.gtf ERR188317_cuff.count/ERR188317.gtf ERR188453_cuff.count/ERR188453.gtf --min_fpkm 1 --min_iso_num 2 --sub 3000 > data/tx_exp.mat`
    This step will subsamples 3000 transcripts from genes which contain at least two isoforms. 
    Each of the transcript has a minimum FPKM value at 1.0 for all samples. 
    
5. Simulate using Polyester.
  `R < polyster.R --no-save`
  Note that this scripts assmue you have made the gtf database produced by R package `AnnotationDbi`, e.g., gencode.v19.annotation.sqlite. Otherwise you have to generate this file by doing 
  `txdb <- makeTxDbFromGFF("/path/to/genes.gtf")` in R and this might take a while.
  
6. Align the simulated reads, convert to BAM and sort them.

   ` hisat2 -x /path/to/hg19.fa -1 simulated_reads/sample_01_1.fasta -2 simulated_reads/sample_01_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_01.sam`
  
   `samtools view -bS simulated_reads/sample_01.sam -o simulated_reads/sample_01.bam`
  
   ` samtools sort simulated_reads/sample_01.bam simulated_reads/sample_01.sorted`
  
   Remember to repeat this step for other simulated samples. 

7. Generate the truth gtf for evaluation. 
   ` awk '{print $2}' data/tx_exp.mat > data/tx_names.lst`
  
   `scripts/gff-select.py /path/to/gencode.v19.annotation.gff3 -l data/tx_names.lst > data/truth.gff`
  
8. Run all three programs on the simulated data.

   `./run_cuff.sh`
  
   `./run_straw.sh`
  
   `./run_string.sh`
  
9. Run Cuffcompare

   `./cuffcompare.sh/`
  
10. Run the evaluation script. This will generate 5 pdf plots. 

    `R < eval.R --no-save`
