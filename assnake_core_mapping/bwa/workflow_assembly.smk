def get_ref_fasta(wildcards):
    ddf = assnake.api.dataset.Dataset(wildcards.ass_df)
    # print(ddf)
    fs_prefix = assnake.api.dataset.Dataset(wildcards.ass_df).fs_prefix
    # final_contigs_wc: '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/final_contigs__{mod}.fa',
    return wc_config['final_contigs_wc'].format(
            fs_prefix = fs_prefix, 
            df = wildcards.ass_df,
            sample_set = wildcards.sample_set,
            mod = wildcards.mod,
            params = 'def',
            assembler = wildcards.assembler,
            assembler_version = wildcards.assembler_version,
            assembler_params = wildcards.assembler_params
            )

rule create_assembly_index_bwa:
    input:
        ref = get_ref_fasta
    output:
        ref_index = os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{assembler_params}/{ass_df}/{sample_set}/final_contigs__{mod}/index.sa')
    params:
        prefix    = os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{assembler_params}/{ass_df}/{sample_set}/final_contigs__{mod}/index')
    log:            os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{assembler_params}/{ass_df}/{sample_set}/final_contigs__{mod}/log.txt')
    benchmark:      os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{assembler_params}/{ass_df}/{sample_set}/final_contigs__{mod}/benchmark.txt')
    conda: 'env_0.7.17.yaml'
    shell: ('''echo -e "INFO: Creating BWA index for {input.ref}"; \n
         bwa index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1; \n
         echo -e "INFO: Finished creating BWA index for {input.ref}\n"''')


rule map_on_assembly_bwa:
    input:
        r1 = wc_config['fastq_gz_R1_wc'],
        r2 = wc_config['fastq_gz_R2_wc'],
        ref_fasta = os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{params}/{ass_df}/{sample_set}/final_contigs__{mod}/index.sa')
    output:
        sam   = '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{assembler}__{assembler_version}__{params}/{ass_df}/{sample_set}/final_contigs__{mod}/{sample}/{preproc}/{sample}.sam'
    params:
        ind_prefix = os.path.join(config['bwa_index_dir'], 'assembly/{assembler}__{assembler_version}__{params}/{ass_df}/{sample_set}/final_contigs__{mod}/index')
    log:      '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{assembler}__{assembler_version}__{params}/{ass_df}/{sample_set}/final_contigs__{mod}/{sample}/{preproc}/log.txt'
    benchmark:'{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{assembler}__{assembler_version}__{params}/{ass_df}/{sample_set}/final_contigs__{mod}/{sample}/{preproc}/benchmark.txt'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    conda: 'env_0.7.17.yaml'
    threads: 8
    shell: ('(bwa mem -M -t {threads} {params.ind_prefix} {input.r1} {input.r2} | /srv/common/bin/samtools view -SF 4 -h > {output.sam}) >{log} 2>&1')
