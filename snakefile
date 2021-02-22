samples = ["REH", "GM12078"]

rule all:
    input:
         {'.'.join(['_'.join(samples), 'rds'])}


rule cellranger_atac_count:
    input: 
        ref_genome="/crex/data/Chromium/cellranger-ATAC-data/1.2.0/rackham/refdata-cellranger-atac-GRCh38-1.2.0",
        fastq="/crex/proj/uppstore2017165/nobackup/scATAC/fastq" 
    output:
        directory("{sample}")
    params: 
        sample="{sample}"
    shell: 
        """
        module load bioinfo-tools
        module load cellranger-ATAC/1.2.0

        cellranger-atac count --id={params.sample} --reference={input.ref_genome} --fastqs={input.fastq} --sample={params.sample}
        """

rule aggr_csv:
    input: 
        expand("{sample}", sample=samples)
    output:
        "aggr.csv"
    shell:
        """
        FILE=aggr.csv

        echo library_id,fragments,cells >> $FILE

        for sample in {input}
        do
            echo $sample,$(realpath $sample/outs/fragments.tsv.gz),$(realpath $sample/outs/singlecell.csv) >> $FILE
        done
        """

rule cellranger_atac_aggr:
    input: 
        aggr_csv = "aggr.csv", 
        ref_genome="/crex/data/Chromium/cellranger-ATAC-data/1.2.0/rackham/refdata-cellranger-atac-GRCh38-1.2.0",
    output: 
        directory({'_'.join(samples)})
    params:
        id = '_'.join(samples)
    shell: 
        """
        module load bioinfo-tools
        module load cellranger-ATAC/1.2.0 
        
        cellranger-atac aggr --id={params.id} --csv={input.aggr_csv} --reference={input.ref_genome} --normalize=depth --nosecondary 
        """

rule prepro_seurat:
    input: 
        {'_'.join(samples)}
    output: 
        {'.'.join(['_'.join(samples), 'rds'])}
    shell:
        """
        module load bioinfo-tools
        module load R/4.0.0 
        module load R_packages/4.0.0 

        Rscript --vanilla atac_preproc_seurat.R {input} {output}
        """
