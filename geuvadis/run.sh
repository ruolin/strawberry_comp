gff-print.py -f ~/Research/geuvadis/alpine/cufflinks/ERR188021_cuff.count/ERR188021.gtf ~/Research/geuvadis/alpine/cufflinks/ERR188052_cuff.count/ERR188052.gtf ~/Research/geuvadis/alpine/cufflinks/ERR188088_cuff.count/ERR188088.gtf ~/Research/geuvadis/alpine/cufflinks/ERR188204_cuff.count/ERR188204.gtf ~/Research/geuvadis/alpine/cufflinks/ERR188317_cuff.count/ERR188317.gtf ~/Research/geuvadis/alpine/cufflinks/ERR188453_cuff.count/ERR188453.gtf --min_fpkm 0.1 --min_iso_num 2 --sub 300 > data/tx_exp.mat

R < polyester.R --no-save

hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_01_1.fasta -2 simulated_reads/sample_01_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_01.sam
hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_02_1.fasta -2 simulated_reads/sample_02_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_02.sam
hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_03_1.fasta -2 simulated_reads/sample_03_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_03.sam
hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_04_1.fasta -2 simulated_reads/sample_04_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_04.sam
hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_05_1.fasta -2 simulated_reads/sample_05_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_05.sam
hisat2 -x ~/Research/Genome/hg19/hg19.fa -1 simulated_reads/sample_06_1.fasta -2 simulated_reads/sample_06_2.fasta --no-mixed --no-discordant -f -k 1 > simulated_reads/sample_06.sam

samtools view -bS simulated_reads/sample_01.sam -o simulated_reads/sample_01.bam
samtools view -bS simulated_reads/sample_02.sam -o simulated_reads/sample_02.bam
samtools view -bS simulated_reads/sample_03.sam -o simulated_reads/sample_03.bam
samtools view -bS simulated_reads/sample_04.sam -o simulated_reads/sample_04.bam
samtools view -bS simulated_reads/sample_05.sam -o simulated_reads/sample_05.bam
samtools view -bS simulated_reads/sample_06.sam -o simulated_reads/sample_06.bam

samtools sort simulated_reads/sample_01.bam simulated_reads/sample_01.sorted
samtools sort simulated_reads/sample_02.bam simulated_reads/sample_02.sorted
samtools sort simulated_reads/sample_03.bam simulated_reads/sample_03.sorted
samtools sort simulated_reads/sample_04.bam simulated_reads/sample_04.sorted
samtools sort simulated_reads/sample_05.bam simulated_reads/sample_05.sorted
samtools sort simulated_reads/sample_06.bam simulated_reads/sample_06.sorted

samtools index simulated_reads/sample_01.sorted.bam
samtools index simulated_reads/sample_02.sorted.bam
samtools index simulated_reads/sample_03.sorted.bam
samtools index simulated_reads/sample_04.sorted.bam
samtools index simulated_reads/sample_05.sorted.bam
samtools index simulated_reads/sample_06.sorted.bam
