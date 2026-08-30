# Genome Assembly, Evaluation and Annotation Workflow

This is a sanitized English command template for long-read genome assembly, chromosome anchoring, polishing and genome annotation. Replace every placeholder path before execution.

No personal usernames, private server paths, project identifiers or embedded note images are included.

## 1. Initial assembly

```bash
mkdir -p hifi && cd hifi
nohup hifiasm -t 190 --telo-m TTTAGGG -o assembly \
  hifi_1.fq.gz hifi_2.fq.gz &

mkdir -p ../hifi_hic && cd ../hifi_hic
nohup hifiasm -t 190 --telo-m TTTAGGG \
  --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz -o assembly \
  hifi_1.fq.gz hifi_2.fq.gz &

mkdir -p ../hifi_ont && cd ../hifi_ont
nohup hifiasm -t 190 --telo-m TTTAGGG \
  --ul ont_reads.fq.gz -o assembly hifi_*.fq.gz &

mkdir -p ../hifi_hic_ont && cd ../hifi_hic_ont
nohup hifiasm -t 190 --telo-m TTTAGGG \
  --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz \
  --ul ont_reads.fq.gz -o assembly hifi_*.fq.gz &
```

For ONT-only assembly, use `hifiasm --ont`. Convert GFA output to FASTA as follows:

```bash
gfatools gfa2fa assembly.gfa > assembly.fa
seqkit stat -a assembly.fa
samtools faidx assembly.fa
```

## 2. Select a primary assembly

```bash
mkdir -p merge && cd merge
ln -s /path/to/hifi/assembly.fa hifi.fa
ln -s /path/to/hifi_hic/assembly.fa hifi_hic.fa
ln -s /path/to/hifi_ont/assembly.fa hifi_ont.fa
ln -s /path/to/hifi_hic_ont/assembly.fa hifi_hic_ont.fa
seqkit stat -a *.fa > assembly_statistics.txt
```

Select the primary assembly using continuity and completeness metrics. The following sections use `hifi_hic_ont.fa` as an example.

## 3. Remove redundant haplotigs

```bash
mkdir -p purge && cd purge
ln -s ../hifi_hic_ont.fa raw.fa
split_fa raw.fa > raw.fa.split
seqkit seq -M 1000000 raw.fa.split > raw.fa.split.small

minimap2 -x asm5 -DP -t 100 raw.fa.split raw.fa.split.small > split.self.paf
minimap2 -I 20G -t 100 -x map-hifi raw.fa \
  /path/to/hifi_reads.fq.gz > hifi.paf
pbcstat hifi.paf
rm -f hifi.paf
calcuts PB.stat > cutoffs 2> calcuts.log
purge_dups -2 -T cutoffs -c PB.base.cov split.self.paf > dups.bed 2> purge_dups.log
get_seqs -e dups.bed raw.fa

mkdir -p a_85 && cd a_85
purge_dups -a 85 -2 -T ../cutoffs -c ../PB.base.cov \
  ../split.self.paf > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../raw.fa
```

## 4. Assembly completeness assessment

```bash
conda activate busco
busco -i /path/to/raw.fa -o busco_raw -m genome \
  -l /path/to/embryophyta_odb10 -c 150
busco -i /path/to/purged.fa -o busco_purge -m genome \
  -l /path/to/embryophyta_odb10 -c 150

compleasm run -a /path/to/raw.fa --threads 150 \
  --autolineage -o compleasm_raw
compleasm run -a /path/to/purged.fa --threads 150 \
  --autolineage -o compleasm_purge
compleasm run -a /path/to/raw.fa --threads 150 \
  -l poales_odb12 -o compleasm_raw_poales_odb12
```

## 5. Remove organelle and other contamination

Use organelle references appropriate for the target species.

```bash
mkdir -p remove_organelle && cd remove_organelle
ln -s /path/to/purged.fa raw.fa
minimap2 -x asm5 -t 100 /path/to/mitochondrion.fasta raw.fa > mitochondrion.paf
minimap2 -x asm5 -t 100 /path/to/chloroplast.fasta raw.fa > chloroplast.paf

for paf in *.paf; do
  python /path/to/filter_polluting_contigs.py "$paf" > "${paf%.paf}.txt"
done
cat *.txt | seqkit grep -v -f - raw.fa > result.fa
```

Repeat the same procedure with documented archaeal, bacterial and viral reference databases to remove other contaminants.

## 6. Hi-C chromosome anchoring and scaffolding

```bash
mkdir -p chromap && cd chromap
ln -s /path/to/remove_organelle/result.fa raw.fa
conda activate chromap
samtools faidx raw.fa
chromap -i -r raw.fa -o raw.index
chromap --preset hic -r raw.fa -x raw.index \
  --remove-pcr-duplicates -1 hic_R1.fq.gz -2 hic_R2.fq.gz \
  --SAM -o aligned.sam -t 50

yahs raw.fa aligned.bam
/path/to/yahs/juicer pre -a -o out_JBAT \
  yahs.out.bin yahs.out_scaffolds_final.agp raw.fa.fai > out_JBAT.log 2>&1
grep PRE_C_SIZE out_JBAT.log | awk '{print $2" "$3}' > asm.size
java -Xmx200G -jar /path/to/juicer_tools.jar pre \
  out_JBAT.txt out_JBAT.hic asm.size
juicer post -o out_JBAT out_JBAT.review.assembly \
  out_JBAT.liftover.agp raw.fa
```

For very large genomes, index parameters may need adjustment. Keep the thread count within the supported range of the installed Chromap version. Verify chromosome orientation and names after scaffolding.

## 7. Chromosome naming

Assign chromosome names using a documented CDS set from a suitable reference species.

```bash
mkdir -p gmap_rename && cd gmap_rename
gmap_build -D . -d DB ../out_JBAT.FINAL.fa
gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3 2> gmap.gff3.log
cut -f1 -d ';' gmap.gff3 | sed 's/ID=//g' | \
  awk -v OFS='\t' '{print $1,$9}' | sort | uniq -c | \
  sort -k1,1rn | head -n 7 > rename.list
awk '{print "s/>"$2"$/>chr"$3"/g"}' rename.list > sed.cmd
sed -f sed.cmd ../out_JBAT.FINAL.fa | seqkit head -n 7 | seqkit sort > tmp.chr.fa
awk '{print $2}' rename.list | seqkit grep -v -f - \
  ../out_JBAT.FINAL.fa | sed 's/>.*/>contig/g' | seqkit rename > tmp.contigs.fa
cat tmp.chr.fa tmp.contigs.fa > result.fa
```

Minimap2 with splice alignment can be used as an alternative to GMAP. Inspect possible whole-chromosome inversions before finalizing the renamed FASTA.

## 8. Gap closing

```bash
mkdir -p gapclose && cd gapclose
ln -s /path/to/chromosome_named/result.fa raw.fa
ragtag.py patch -o patch1 raw.fa /path/to/hifi.fa
ragtag.py patch -o patch2 patch1/ragtag.patch.fasta /path/to/hifi_hic.fa
ragtag.py patch -o patch3 patch2/ragtag.patch.fasta /path/to/ont.fa
ln -s patch3/ragtag.patch.fasta result.fa

samtools faidx result.fa
samtools faidx raw.fa
paste raw.fa.fai result.fa.fai | \
  awk '{print "s/>"$1"$/>chr"$6"/g"}' > sed.cmd
sed -i -f sed.cmd result.fa
```

Check for altered FASTA headers after patching. Gap closing may be omitted when only one assembly version is available; `tgsgapcloser` is an alternative for relatively small genomes.

## 9. Genome polishing

### NextPolish2

```bash
conda activate nextpolish
meryl count k=15 output merylDB raw.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -I 10G -t 50 -W repetitive_k15.txt -ax map-pb \
  raw.fa /path/to/hifi.fasta.gz | samtools sort -@ 50 \
  -o hifi.map.sorted.bam
samtools index hifi.map.sorted.bam

./yak/yak count -o k21.yak -k 21 -b 37 \
  <(zcat /path/to/survey_R*.clean.fastq.gz) \
  <(zcat /path/to/survey_R*.clean.fastq.gz)
./yak/yak count -o k31.yak -k 31 -b 37 \
  <(zcat /path/to/survey_R*.clean.fastq.gz) \
  <(zcat /path/to/survey_R*.clean.fastq.gz)
nextPolish2 -t 5 hifi.map.sorted.bam raw.fa k21.yak k31.yak > polished.fa
```

The process-substitution syntax requires Bash. For limited memory, split the assembly and alignment by chromosome.

### Racon

```bash
meryl count k=15 output merylDB raw.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -I 20G -t 50 -W repetitive_k15.txt -x map-pb \
  raw.fa /path/to/HiFi/*.gz -o hifi.paf
cat /path/to/HiFi/*.gz > hifi.fq.gz
racon -t 60 hifi.fq.gz hifi.paf raw.fa > racon.fa
```

## 10. Reference indexes and gene annotation

```bash
mkdir -p /path/to/results/annotation/{TE,RNA}
cd /path/to/results
mkdir -p index
ln -s /path/to/polished.fa genome.fa
bwa index -p index/genome genome.fa
hisat2-build -p 60 genome.fa index/genome
```

Trim paired-end RNA-seq reads, map them with HISAT2, convert SAM to sorted BAM, and create a two-column BAM list for EviAnn.

```bash
for r1 in /path/to/RNAseq/*R1.fq.gz; do
  trim_galore -q 20 -j 60 --phred33 --stringency 3 --length 20 \
    -e 0.1 --paired "$r1" "${r1/R1/R2}" --gzip -o "$(dirname "$r1")"
done

for r1 in /path/to/RNAseq/*1P.fq.gz; do
  sam=$(basename "${r1/_trim2*/.sam}")
  hisat2 -p 100 -x /path/to/index/genome -1 "$r1" \
    -2 "${r1/1P/2P}" -S "$sam" 2> "${sam%.sam}.log"
done

for sam in *.sam; do
  samtools view -bS -@ 20 "$sam" | \
    samtools sort -@ 20 -o "${sam%.sam}.bam"
done
for bam in "$PWD"/*.bam; do printf '%s\tbam\n' "$bam"; done > bam.txt

conda activate EviAnn
/path/to/EviAnn/bin/eviann.sh -t 150 \
  -g /path/to/genome.fa \
  -p /path/to/homologous_species_proteins.fasta \
  -r /path/to/bam.txt \
  -s /path/to/uniprot_sprot.fasta
```

## 11. Repeat, tRNA and ncRNA annotation

```bash
conda activate HiTE
python /path/to/HiTE/main.py --genome /path/to/genome.fa \
  --outdir /path/to/annotation/TE --thread 120 --plant 1 \
  --curated_lib /path/to/curated_repeat_library.fa \
  --recover 1 --annotate 1 --intact_anno 1 > HiTE.log 2>&1

conda install -c bioconda trnascan-se

wget -c https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmsearch --cpu 20 --tblout genome.ncRNA.tblout \
  /path/to/Rfam.cm /path/to/genome.fa > genome.cmsearch
```

## Reproducibility and safety checklist

- Replace every placeholder path before execution.
- Record software versions, reference assemblies, database releases and parameters.
- Confirm intermediate file extensions match their actual formats.
- Validate chromosome names and orientations after scaffolding, patching and polishing.
- Do not publish raw reads, credentials, internal hostnames, private filesystem paths or personal identifiers.
