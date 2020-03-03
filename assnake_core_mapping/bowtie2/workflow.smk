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
        r1 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R2.fastq.gz',
        ref_index = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index.1.bt2')
    output:
        sam = '{fs_prefix}/{df}/mapped/bowtie2__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.sam'
    params:
        ind_prefix = os.path.join(index_dir, 'bowtie2/{path}/{seq_set_id}/index')
    log:      '{fs_prefix}/{df}/mapped/bowtie2__{params}/{path}/{seq_set_id}/{sample}/{preproc}/log.txt'
    benchmark:'{fs_prefix}/{df}/mapped/bowtie2__{params}/{path}/{seq_set_id}/{sample}/{preproc}/benchmark.txt'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    threads: 8
    conda: 'env.yaml'
    shell: '''export PERL5LIB='';\n
    bowtie2 -x {params.ind_prefix} -1 {input.r1} -2 {input.r2} -S {output.sam} -p {threads} -D 20 -R 5 -N 1 -L 20 -i S,1,2.50 > {log} 2>&1'''