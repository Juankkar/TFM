# Trabajo de Fin de Máster

## Universidad Internacional de Valencia (VIU)

### Alumno: Juan Carlos García Estupiñán

Tema Cáncer de colon

## Comandos para el preprocesado de los BAM

* BAM index

``` 
samtools index input.bam
```

* Chromosome Filtering

```
samtools view -b input.bam chr{5} > output.bam
```

* Sorted filter BAM file

```
samtools sort -n output.bam -o output_sorted.bam
```

* BAM file to a fastq file

```
samtools fastq -@ 8 IV_method1.bam \
    -1 ../../raw/IV_method1_1.fastq.gz \
    -2 ../../raw/IV_method1_2.fastq.gz \
    -0 /dev/null -s /dev/null -n
```

## Descargar el genoma de referencia en data/reference

```
wget -O genome.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.15.fa.gz

gzip -d genome.fa.gz
```
