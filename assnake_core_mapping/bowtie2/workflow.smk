import os

index_dir = config['bwa_index_dir']
fna_db_dir= config['fna_db_dir']

rule create_seq_set_index_bowtie2:
    input:
        ref       = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fna')
    output:
        ref_index = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index.1.bt2')
    params:
        prefix    = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index'),
        index_dir = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/'),
        seed      = 42
    log:            os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/log.txt')
    benchmark:      os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/benchmark.txt')
    threads: 8
    conda: 'env.yaml'
    shell: ('''mkdir -p {params.index_dir}; export PERL5LIB='';\n
    bowtie2-build --threads {threads} --seed {params.seed} {input.ref} {params.prefix} > {log} 2>&1''')
        

rule map_on_ref_bowtie2:
    input:
        r1 = '{fs_prefix}/{df}/reads/{preproc}/{df_sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}/{df_sample}_R2.fastq.gz',
        ref_index = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index.1.bt2'),
        params =  os.path.join(config['assnake_db'], 'params/bowtie2/{params}.yaml'),
    output:
        sam = '{fs_prefix}/{df}/mapped/bowtie2__{params}__{version}/{path}/{seq_set_id}/{df_sample}/{preproc}/{df_sample}.sam'
    params:
        ind_prefix = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index'),
        r_s = '{fs_prefix}/{df}/reads/{preproc}/{df_sample}_S.fastq.gz'
    log:      '{fs_prefix}/{df}/mapped/bowtie2__{params}__{version}/{path}/{seq_set_id}/{df_sample}/{preproc}/log.txt'
    benchmark:'{fs_prefix}/{df}/mapped/bowtie2__{params}__{version}/{path}/{seq_set_id}/{df_sample}/{preproc}/benchmark.txt'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    threads: 20
    conda: 'env.yaml'
    wrapper: "file://"+os.path.join(config['assnake-core-mapping'], 'bowtie2/wrapper.py')
# bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x metaphlan_databases/mpa_v295_CHOCOPhlAn_201901 -U metagenome.fastq