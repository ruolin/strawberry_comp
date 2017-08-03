cuffcompare -r data/truth.gff cufflinks/1/transcripts.gtf -o cufflinks/1/cufflinks
cuffcompare -r data/truth.gff cufflinks/2/transcripts.gtf -o cufflinks/2/cufflinks
cuffcompare -r data/truth.gff cufflinks/3/transcripts.gtf -o cufflinks/3/cufflinks
cuffcompare -r data/truth.gff cufflinks/4/transcripts.gtf -o cufflinks/4/cufflinks
cuffcompare -r data/truth.gff cufflinks/5/transcripts.gtf -o cufflinks/5/cufflinks
cuffcompare -r data/truth.gff cufflinks/6/transcripts.gtf -o cufflinks/6/cufflinks

cuffcompare -r data/truth.gff stringtie/1/transcripts.gtf -o stringtie/1/stringtie
cuffcompare -r data/truth.gff stringtie/2/transcripts.gtf -o stringtie/2/stringtie
cuffcompare -r data/truth.gff stringtie/3/transcripts.gtf -o stringtie/3/stringtie
cuffcompare -r data/truth.gff stringtie/4/transcripts.gtf -o stringtie/4/stringtie
cuffcompare -r data/truth.gff stringtie/5/transcripts.gtf -o stringtie/5/stringtie
cuffcompare -r data/truth.gff stringtie/6/transcripts.gtf -o stringtie/6/stringtie

cuffcompare -r data/truth.gff strawberry/1/assembled_transcripts.gtf -o strawberry/1/strawberry
cuffcompare -r data/truth.gff strawberry/2/assembled_transcripts.gtf -o strawberry/2/strawberry
cuffcompare -r data/truth.gff strawberry/3/assembled_transcripts.gtf -o strawberry/3/strawberry
cuffcompare -r data/truth.gff strawberry/4/assembled_transcripts.gtf -o strawberry/4/strawberry
cuffcompare -r data/truth.gff strawberry/5/assembled_transcripts.gtf -o strawberry/5/strawberry
cuffcompare -r data/truth.gff strawberry/6/assembled_transcripts.gtf -o strawberry/6/strawberry

scripts/parse-stats.py -f cufflinks/1/cufflinks.stats cufflinks/2/cufflinks.stats cufflinks/3/cufflinks.stats cufflinks/4/cufflinks.stats cufflinks/5/cufflinks.stats cufflinks/6/cufflinks.stats > cufflinks/assemb_cuff.stats
scripts/parse-stats.py -f strawberry/1/strawberry.stats strawberry/2/strawberry.stats strawberry/3/strawberry.stats strawberry/4/strawberry.stats strawberry/5/strawberry.stats strawberry/6/strawberry.stats > strawberry/assemb_straw.stats
scripts/parse-stats.py -f stringtie/1/stringtie.stats stringtie/2/stringtie.stats stringtie/3/stringtie.stats stringtie/4/stringtie.stats stringtie/5/stringtie.stats stringtie/6/stringtie.stats > stringtie/assemb_string.stats
