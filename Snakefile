rule dl_rawdata:
    input:
        script = "code/1dl_rawdata.sh"
    output: 
        "data/report/filereport_read_run_PRJEB7926_tsv.txt"
    params:
        file = "filereport_read_run_PRJEB7926_tsv.txt"
    shell:
        """
        {input.script} {params.file}
        """