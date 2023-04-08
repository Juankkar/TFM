# Trabajo de Fin de Máster

## Universidad Internacional de Valencia (VIU)

### Alumno: Juan Carlos García Estupiñán

Tema Cáncer de colorrectal.

## DATA sacada de un artículo en un repositorio público: 
* ### Identificador: [PRJEB7926](https://www.ebi.ac.uk/ena/browser/view/PRJEB7926)
* ### Artículo: [*GREM1* and POLE variants in hereditary colorectal cancer syndromes](https://onlinelibrary.wiley.com/doi/10.1002/gcc.22314)

## Comandos para el preprocesado de los BAM

* BAM index

``` 
samtools index input.bam
```

* Chromosome Filtering

```
samtools view -b input.bam chr{n} > output.bam
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
