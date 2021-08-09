[![DOI](https://zenodo.org/badge/390445435.svg)](https://zenodo.org/badge/latestdoi/390445435)

# Estuarine molecular bycatch as a landscape-wide biomonitoring tool

#### Code and data for:

Mariani, S., Harper, L.R., Collins, R.A., Baillie, C., Wangensteen, O.S., McDevitt, A.D., Heddell-Cowie, M., Genner, M.J. (2021). Estuarine molecular bycatch as a landscape-wide biomonitoring tool. _Insert Journal_. [https://doi.org/xxx](https://doi.org/xxx).

---

1. Clone this project repository onto your system and install all required R package versions.

```bash
# download repo and R packages
git clone https://github.com/genner-lab/molecular-bycatch.git
cd molecular-bycatch
Rscript -e "renv::restore()"
```

2. Download the 'meta-fish-pipe' bioinformatics module. Be sure to read instructions and install all system software as instructed at [github.com/genner-lab/meta-fish-pipe](https://github.com/genner-lab/meta-fish-pipe).

```bash
# download and install pipeline
git clone https://github.com/genner-lab/meta-fish-pipe.git
cd meta-fish-pipe
Rscript -e "renv::restore()"
cd ..
```

3. Download the custom UK fish reference library from [github.com/genner-lab/meta-fish-lib](https://github.com/genner-lab/meta-fish-lib).

```bash
# download the custom uk fish reference library
scripts/download-reference-library.R
```

4. Download the NCBI RefSeq reference library from [github.com/genner-lab/refseq-reflib](https://github.com/genner-lab/refseq-reflib).

```bash
# download RefSeq reference library
git clone https://github.com/genner-lab/refseq-reflib.git
cd refseq-reflib
Rscript -e "renv::restore()"
mkdir temp references
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
scripts/download.sh
scripts/extract.R -p mifish-u
scripts/annotate.R -s 42 -p mifish-u
rm temp/duckdb
cd ..
```

5. Download the fastq data from Dryad (mifish-u) and SRA (tele02).

```bash
# make dir
mkdir -p temp/data
# download dryad data from doi: 10.5061/dryad.b8f6s44
scripts/rdryad.R
# get SRA data
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra69/SRZ/014350/SRR14350412/SeaDNA_Teleo02_01_S1_L001_R1_001.fastq.gz -P temp/data
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra69/SRZ/014350/SRR14350412/SeaDNA_Teleo02_01_S1_L001_R2_001.fastq.gz -P temp/data
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra64/SRZ/014350/SRR14350411/SeaDNA_Teleo02_02_S2_L001_R1_001.fastq.gz -P temp/data
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra64/SRZ/014350/SRR14350411/SeaDNA_Teleo02_02_S2_L001_R2_001.fastq.gz -P temp/data
# check md5sums
md5sum -c assets/fastq.md5
```


6. Copy sample files into the meta-fish-pipe module.

```bash
# copy across files to bioinformatics pipeline
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv
cp assets/sequencing-master.csv meta-fish-pipe/assets/sequencing-master.csv
cp refseq-reflib/references/refseq207-annotated-mifish-u.csv meta-fish-pipe/assets/refseq207-annotated-mifish-u.csv
cp assets/assemble-results.R meta-fish-pipe/scripts/assemble-results.R
```

7. Generate session information and set up the bioinformatics pipeline.

```bash
# generate session stats
cd meta-fish-pipe
scripts/session-info.sh -r assets/refseq207-annotated-mifish-u.csv -c assets/meta-fish-lib-v244.csv
```

8. Prepare the directories for each library.

```bash 
# prepare each library
scripts/prepare-libraries.sh -p mifish-u -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib2
```

9. Copy across symbolic links to the fastq files.

```bash
# make symlinks MIFISH-U LIB1
ln -s -r ../temp/data/12S-mifishu-R1.fastq.gz temp/processing/mifish-u-lib1/fastq/R1.fastq.gz
ln -s -r ../temp/data/12S-mifishu-R2.fastq.gz temp/processing/mifish-u-lib1/fastq/R2.fastq.gz
# make symlinks TELE02 LIB1
ln -s -r ../temp/data/SeaDNA_Teleo02_01_S1_L001_R1_001.fastq.gz temp/processing/tele02-lib1/fastq/R1.fastq.gz
ln -s -r ../temp/data/SeaDNA_Teleo02_01_S1_L001_R2_001.fastq.gz temp/processing/tele02-lib1/fastq/R2.fastq.gz
# make symlinks TELE02 LIB2
ln -s -r ../temp/data/SeaDNA_Teleo02_02_S2_L001_R1_001.fastq.gz temp/processing/tele02-lib2/fastq/R1.fastq.gz
ln -s -r ../temp/data/SeaDNA_Teleo02_02_S2_L001_R2_001.fastq.gz temp/processing/tele02-lib2/fastq/R2.fastq.gz
```

10. Generate barcode files for each library.

```bash
# generate barcodes files
scripts/generate-barcodes.R -p mifish-u -l lib1 -f 21 -r 27 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib2 -f 18 -r 20 -m assets/sequencing-master.csv
```

11. Demultiplex each library with cutadapt.

```bash
# demultiplex with cutadapt
scripts/demultiplex.sh -p mifish-u -l lib1 -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG -t 8 -m 21
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib2 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
```

12. Denoise each library with dada2.

```bash
# denoise with dada2
scripts/dada2.R -p mifish-u -l lib1
scripts/dada2.R -p tele02 -l lib1
scripts/dada2.R -p tele02 -l lib2
```

13. Generate pipeline output stats.

```bash
# generate bioinfomatics stats
scripts/generate-stats.sh -p mifish-u -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib2 -t 8
```

14. Run the taxonomic assignment steps and combine results from each library 

```bash
# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p mifish-u
# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```

15. Copy results to 'molecular-bycatch' directory.

```bash
# copy dir
cd ..
cp -r meta-fish-pipe/results results
```
