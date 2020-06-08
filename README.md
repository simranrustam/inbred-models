# research-theis
++ pwd
++ rev
++ tr / '\n'
++ cut -d / -f-2
++ rev
++ grep RP
+ PROJ=
+ TSV=.tsv
++ echo .tsv
++ cut -d . -f1
+ SRP=
++ cut -f2 .tsv
++ sed 1d
++ sort -u
cut: .tsv: No such file or directory
+ ls '*bam'
+ parallel samtools index '{}'
ls: cannot access '*bam': No such file or directory
Academic tradition requires you to cite works you base your article on.
When using programs that use GNU Parallel to process data for publication
please cite:

  O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
  ;login: The USENIX Magazine, February 2011:42-47.

This helps funding further development; AND IT WON'T COST YOU A CENT.
If you pay 10000 EUR you should feel free to use GNU Parallel without citing.

To silence this citation notice: run 'parallel --citation'.

+ featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e6.bed.saf -o .1e6.tsv 'SRX*.bam'

ERROR: invalid parameter: 'SRX*.bam'

+ sed 1d .1e6.tsv
+ cut -f1,7-
sed: can't read .1e6.tsv: No such file or directory
+ featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e5.bed.saf -o .1e5.tsv 'SRX*.bam'

ERROR: invalid parameter: 'SRX*.bam'

+ sed 1d .1e5.tsv
+ cut -f1,7-
sed: can't read .1e5.tsv: No such file or directory
+ featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e4.bed.saf -o .1e4.tsv 'SRX*.bam'

ERROR: invalid parameter: 'SRX*.bam'

+ sed 1d .1e4.tsv
+ cut -f1,7-
sed: can't read .1e4.tsv: No such file or directory
+ featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e3.bed.saf -o .1e3.tsv 'SRX*.bam'

ERROR: invalid parameter: 'SRX*.bam'

+ sed 1d .1e3.tsv
+ cut -f1,7-
sed: can't read .1e3.tsv: No such file or directory
+ featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/Mus_musculus.GRCm38.98.gtf.saf -o .genes.tsv 'SRX*.bam'

ERROR: invalid parameter: 'SRX*.bam'

+ sed 1d .genes.tsv
+ cut -f1,7-
