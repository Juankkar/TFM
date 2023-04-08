# Trabajo de Fin de Máster

## Universidad Internacional de Valencia (VIU)

### Alumno: Juan Carlos García Estupiñán

Tema Cáncer de colorrectal.

## DATA sacada de un artículo en un repositorio público: 
* ### Identificador: [PRJEB7926](https://www.ebi.ac.uk/ena/browser/view/PRJEB7926)
* ### Artículo: [*GREM1* and POLE variants in hereditary colorectal cancer syndromes](https://onlinelibrary.wiley.com/doi/10.1002/gcc.22314)

## **Directorios**

### **[Código](code) empleado en la carpeta**

* A continuación enumero los scripts y su función:
    * [enviroments/TFM.yml](code/enviroments/TFM.yml): directorio secundario con el archivo *YML* con el ambiente conda para reproducir el repositorio.
    * [1dl_rawdata.sh](code/1dl_rawdata.bash): script que descarga los bams del repositorio público, está vinculado para su funcionamiento a los [metadatos (report.tsv)](metadata/report.tsv).
    * [2join_samfiles.sh](code/2join_samfiles.sh): unimos los archivos *SAM* que han sido obtenidos al pasar los *fastq.gz* a este formato.
    * [3sam_to_bam.sh](code/3sam_to_bam.sh): convertimos el archivo SAM unido anterior y lo transformamos en un *BAM*, además de indexarlo.
    * [4delete_duplicates](code/4delete_duplicates.sh): eliminamos los duplicados.
    * [5extracting_variants.sh](code/5extracting_variants.sh): extraemos los varientes en un archivo *VCF* tanto SNPs como INDELs, además hacemos filtros para obtener para separar ambos tipos de variantes.
    * [6samples.ipynb](code/6samples.ipynb): jupyter notebook que extrae información de la tabla complementaria del artículo para tener más información, guardamos ticha tabla procesada en los [metadatos (table3.csv)](metadata/table3.csv)

### **[Data](data) usada**

* Describibos los directorios adicionales de los que se hace uso:
    * [original_bam](data/original_bam/): en este repositorio se guardan los *BAM* del artículo que se descargan con el script 1.
        * En este además hay otro directorio [filtering](data/original_bam/filtering/), donde guardamos los bam que hemos filtrado según el cromosoma de interés. El preprocesado el siguiente (https://www.metagenomics.wiki/tools/samtools/converting-bam-to-fastq) 


* Comandos para el preprocesado de los BAM: 
1. BAM index

``` 
samtools index input.bam
```


2. Chromosome Filtering


```
samtools view -b input.bam chr{n} > output.bam
```

3. Sorted filter BAM file

```
samtools sort -n output.bam -o output_sorted.bam
```

4. BAM file to a fastq file

```
samtools fastq -@ 8 IV_method1.bam \
    -1 ../../raw/IV_method1_1.fastq.gz \
    -2 ../../raw/IV_method1_2.fastq.gz \
    -0 /dev/null -s /dev/null -n
```

* [raw](data/raw/): aquí es donde se han guardado las lecturas forward y reverse "crudas", en nuestro caso preprocesadas de los archivos *BAM* en el paso anterior.
* [processed](data/processed/): donde procesamos los datos de fastq medinate ```fastp```.
* [reference](data/reference/): donde guardamos el genoma de referencia, para probarlo todo he elegido estudiar el cromosoma 15 de una de las muestras
    * Comando para descargar el genoma: 

```
wget -O genome.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.15.fa.gz

gzip -d genome.fa.gz
```

### **[Metadatos](metadata)**

* [report.tsv](metadata/report.tsv): archivo con información y links de descarga de los datos proporcionados por el repositorio público. Se hace uso del campo 8 para descargar los datos en el script [1dl_rawdata.sh](code/1dl_rawdata.bash).
* [table3.csv](metadata/table3.csv): tabla que proporcionan los autores para una información más detallada.

### **[Resultados](results)**

* Aquí se van creando y guardando los rsultados que se van obteniendo al ejecurtar los scirpts

### **Archivos adicionales de control de flujo de trabajo**

* [Snakefile](Snakefile): archivo de snakemake para ejecutar el código y mantener el flujo de trabajo.

* [config.yaml](config.yaml): archivo que condicional al anterior.