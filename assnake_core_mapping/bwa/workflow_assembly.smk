index_dir = config['bwa_index_dir']

def get_ref_fasta(wildcards):
    ddf = assnake.api.dataset.Dataset(wildcards.ass_df)
    # print(ddf)
    fs_prefix = assnake.api.dataset.Dataset(wildcards.ass_df).fs_prefix
    final_contigs_wc = '{fs_prefix}/{df}/assembly/{sample_set}/{assembler}__{assembler_version}__{params}/final_contigs__{mod}.fa'

    return final_contigs_wc.format(
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
        ref_index = os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/index.sa')
    params:
        prefix    = os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/index')
    log:            os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/log.txt')
    benchmark:      os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/benchmark.txt')
    conda: 'env_0.7.17.yaml'
    shell: ('''bwa index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1;''')


rule map_on_assembly_bwa:
    input:
        r1 = wc_config['fastq_gz_R1_wc'],
        r2 = wc_config['fastq_gz_R2_wc'],
        ref_fasta = os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/index.sa')
    output:
        sam   = '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/{sample}/{preproc}/{sample}.bam'
    params:
        ind_prefix = os.path.join(index_dir, 'assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/index')
    log:        '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/{sample}/{preproc}/log.txt'
    benchmark:  '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/assembly/{ass_df}/{sample_set}/{assembler}__{assembler_version}__{assembler_params}/final_contigs__{mod}/{sample}/{preproc}/benchmark.txt'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    conda: 'env_0.7.17.yaml'
    threads: 12
    shell: ('(bwa mem -M -t {threads} {params.ind_prefix} {input.r1} {input.r2} | samtools sort -@{threads} -o {output.sam} - ) >{log} 2>&1')
