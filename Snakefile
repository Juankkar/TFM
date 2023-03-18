rule dl_rawdata:
    input:
        "data/report/report.txt"
    output: 
        "data/raw/{sample}1.fastq.gz",
        "data/raw/{sample}2.fastq.gz"
    shell:
        ## de report.txt descargamos los datos a partir de los ftp
        ## de las lecturas
        """reads1=$(cat {input} \
        | cut -f7 \
        | cut -d ";" -f 1 \
        | tail -n 1)
        
        reads2=$(cat {input} \
        | cut -f7 \
        | cut -d ";" -f 2 \
        | tail -n 1)
        
        reads1_array=$(echo $reads1)
        reads2_array=$(echo $reads2)
        
        for r1 in $reads1_array
        do
            for r2 in $reads2_array
            do
                wget -O {output[0]} $r1
                wget -O {output[1]} $r2
            done
        done"""