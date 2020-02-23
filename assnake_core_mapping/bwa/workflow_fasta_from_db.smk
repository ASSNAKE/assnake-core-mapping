import os

index_dir = config['bwa_index_dir']
fna_db_dir= config['fna_db_dir']

rule create_seq_set_index_bwa:
    input:
        ref       = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fna')
    output:
        ref_index = os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/index.sa')
    params:
        prefix    = os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/index')
    log:            os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/log.txt')
    benchmark:      os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/benchmark.txt')
    conda: 'env_0.7.17.yaml'
    shell: ('''echo -e "INFO: Creating BWA index for {input.ref}"; \n
         bwa index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1; \n
         echo -e "INFO: Finished creating BWA index for {input.ref}\n"''')
        

rule map_on_ref_bwa:
    input:
        r1 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R2.fastq.gz',
        ref_index = os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/index.sa')
    output:
        sam = '{fs_prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.sam'
    params:
        ind_prefix = os.path.join(index_dir, 'bwa/{path}/{seq_set_id}/index')
    log:      '{fs_prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/log.txt'
    benchmark:'{fs_prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/benchmark.txt'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    threads: 8
    conda: 'env_0.7.17.yaml'
    shell: ('(bwa mem -M -t {threads} {params.ind_prefix} {input.r1} {input.r2} > {output.sam}) >{log} 2>&1')





        


        

