fna_db_dir= config['fna_db_dir']

rule init_bam_mapped:
    input:
        sam = '{fs_prefix}/{df}/mapped/{mapper}__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.sam'
    output:
        bam = '{fs_prefix}/{df}/mapped/{mapper}__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.bam'
    params:
        tmp = '{fs_prefix}/{df}/mapped/{mapper}__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.tmp.bam'
    conda: './bwa/env_0.7.17.yaml'
    shell: 
        '''echo "{input.sam}";\n
        samtools view -bS {input.sam} -o {params.tmp};\n
        samtools sort {params.tmp} -o {output.bam};\n
        samtools index {output.bam};\n
        rm {params.tmp}'''

rule samtools_flagstat:
    input:  bam      = '{fs_prefix}/{df}/mapped/{mapper}__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.bam'
    output: flagstat = '{fs_prefix}/{df}/mapped/{mapper}__{params}/{path}/{seq_set_id}/{sample}/{preproc}/{sample}.flagstat.txt'
    conda: './bwa/env_0.7.17.yaml'
    shell: 'samtools flagstat {input.bam} > {output.flagstat}'


nucl_dir= config['fna_db_dir']
rule coverage:
    input:
        bam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.bam',
        ref = nucl_dir+"{type}/{fasta}.fa"
    output:
        cov = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.cov'
    run:
        shell('/srv/common/bin/genomeCoverageBed -bga -ibam {input.bam} > {output.cov}')

def get_reference_fasta(wildcards):
    """
    Reconstructs reference fasta location
    """
    path = wildcards.path.replace('___', '/')
    
    fasta_wc_loc = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    fasta = fasta_wc_loc.format(path=path, seq_set_id=wildcards.seq_set_id)

    return fasta


             

rule leave_aligned_bwa:
    input: 
        sam = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/{id_sample}.sam'
    output: 
        filtered = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/{id_sample}.aligned.sam'
    run:
        shell('samtools view -S -F 4 {input.sam} > {output.filtered}')
        

rule filter_mapping:
    input: 
        sam = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/mapped/{id_sample}.sam'
    output: 
        filtered = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/filtered/{id_sample}.sam'
    benchmark: 'time/filter-map/{df}/{preproc}/{type}/{id_seq_set}/{id_sample}.time'
    run:
        shell('{config[python.bin]} bin/scripts/mgSNP_sam-filter.py -i {input.sam} -o {output.filtered}.pre -l 25 -m 90')
        shell('grep -v "XA:" {output.filtered}.pre > {output.filtered}')
        shell('rm {output.filtered}.pre')

