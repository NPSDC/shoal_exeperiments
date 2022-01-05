from os.path import join
import yaml
configfile: "config.json"

with open('cluster_partition.conf', 'r') as f:
    partitions = yaml.safe_load(f)

def get_partition(rule_name):
    return partitions[rule_name] if rule_name in partitions else partitions['__default__']


out_reads = config['out_reads']
reads_dir = config['reads_path']
REPLICATES = config['replicates']

raw_reads_gz = join(reads_dir, 'out_{{rep}}', 'sample_0{{cond}}_{end}_shuffled.fa.gz')
L_reads_gz = raw_reads_gz.format(end=1)
R_reads_gz = raw_reads_gz.format(end=2)

raw_reads_fasta = join(reads_dir, 'out_{{rep}}', 'sample_0{{cond}}_{end}_shuffled.fasta')
L_reads_fa = raw_reads_fasta.format(end=1)
R_reads_fa = raw_reads_fasta.format(end=2)

output_dir_fmt = join(out_reads, '{rep}_{cond}')
SEED = 2021
#print(extend(output_dir_fmt, rep=REPLICATES, cond=[1,2]))

train_fasta_fmt = join(config["out_reads"], '{{rep}}_{{cond}}', 'train', \
    '{{fold}}', 'sample_0{{cond}}_{end}_shuffled.fasta')
train_fasta_l = train_fasta_fmt.format(end=1)
train_fasta_r = train_fasta_fmt.format(end=2)

test_fasta_fmt = join(config["out_reads"], '{{rep}}_{{cond}}', 'test', \
    '{{fold}}', 'sample_0{{cond}}_{end}_shuffled.fasta')
test_fasta_l = test_fasta_fmt.format(end=1)
test_fasta_r = test_fasta_fmt.format(end=2)

FOLDS = list(range(1, config['k']+1))
vbpriors = [1e4, 1e3, 1e2, 1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
shoal_params = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]
train_sal_out = join(config["out_sal_train"], '{rep}_{cond}', 'vb={vb}', \
    '{fold}', 'quant.sf')
train_sal_out_dir = join(config["out_sal_train"], '{rep}_{cond}', 'vb={vb}', \
    '{fold}')
train_sal_out_vb_dir = join(config["out_sal_train"], '{rep}_{cond}', 'vb=1', \
    '{fold}')
test_sal_out_dir = join(config["out_sal_test"], '{rep}_{cond}', \
    '{fold}')

comm_path=join("aux_info", "eq_classes.txt.gz")
perp_file_sal=join(config["perp_sal_path"], "vb={vb}", '{rep}_{cond}', '{fold}', 'perplexity.yml')
perp_file_shoal=join(config["perp_shoal_path"], "vb=1e0", "c={w}", '{rep}_{cond}', '{fold}', 'perplexity.yml')
train_sal_ecs=join(train_sal_out_dir, comm_path)
test_sal_ecs=join(test_sal_out_dir, comm_path)
train_shoal_ecs=join(train_sal_out_vb_dir, comm_path)

shoal_prior = join(config["out_shoal"], "vbprior=1e0", "{fold}", "prior.tsv")
shoal_quant_ind = join(config["out_shoal"], "vbprior=1e0", "{fold}", "{rep}_{cond}_c={w}_adapt.sf")
shoal_params = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]
tried_priors = [0, 0.00001, 0.1, 100, 100000]

sal_best_quant_dir = join(config["sal_best_quant"], "vbprior=1e0", "{rep}_{cond}")
sal_best_quant = join(config["sal_best_quant"], "vbprior=1e0", "{rep}_{cond}", "quant.sf")
shoal_best_quant = join(config["shoal_best_quant"], "vbprior=1e0", "c={w}", "{rep}_{cond}", "adapt_quant.sf")
shoal_best_prior = join(config["shoal_best_quant"], "vbprior=1e0", "c={w}", "prior.tsv")


rule shoal_quant_bestc:
    input:
        expand(shoal_best_quant, w=tried_priors, rep=REPLICATES, cond=[1,2])
    params:
        partition = get_partition('shoal_quant_bestc')

rule _shoal_quant_bestc:
    input:
        prior=shoal_best_prior,
        quants_dir=sal_best_quant_dir
    output:
        shoal_best_quant
    params:
        shoal_path=config["shoal_path"],
        partition = get_partition('_shoal_quant_bestc'),
        
    shell:
        '''
        {params.shoal_path} requant-v2\
            --prior {input.prior}\
            --sample {input.quants_dir}\
            --output {output}\
            --weight {wildcards.w}\
            -u 1\
            -f true
        '''

rule run_shoal_best_prior:
    input:expand(shoal_best_prior,w=tried_priors),
    params:partition = get_partition('run_shoal_best_prior')

rule _run_shoal_best_prior:
    input:
        expand(sal_best_quant, cond=[1,2], rep=REPLICATES)
    output:
        shoal_best_prior
    params:
        partition = get_partition('_run_shoal_best_prior'),
        shoal_path=config["shoal_path"]
    shell:
        '''
        {params.shoal_path} create-prior\
            --samples {input} \
            --output {output}
        '''

rule sal_quant_bestvb:
    input:
        expand(sal_best_quant, rep=REPLICATES, cond=[1,2])
    params:
        partition = get_partition('sal_quant_bestvb')

rule _sal_quant_bestvb:
    input:
        inp_fastq1 = L_reads_gz,
        inp_fastq2 = R_reads_gz,
        index = join(config['sal_ind_path'], "pos.bin")
    output:
        out_sal = sal_best_quant
    resources: cpus=10, mem=32000
    params:
        out_dir = sal_best_quant_dir,
        partition = get_partition('_sal_quant_bestvb'),
        sal_path = config['sal_path'],
        sal_ind = config['sal_ind_path']
    shell:
        '''
        {params.sal_path} quant -i {params.sal_ind}\
             -l A -p {resources.cpus} --gcBias\
                -d\
                --vbPrior 1\
                --output {params.out_dir} -1 {input.inp_fastq1} -2 {input.inp_fastq2}    
        '''


rule run_perplexity_shoal:
    input:expand(perp_file_shoal, fold=FOLDS, w=shoal_params, cond=[1,2],rep=REPLICATES)
    params:
        partition = get_partition('perplexity_shoal')

rule _run_perplexity_shoal:
    input:
        train_quant=shoal_quant_ind,
        train_ecs=train_shoal_ecs,
        test_ecs=test_sal_ecs
    output:
        perp_file=perp_file_shoal
    params:
        perp_path=config["perplexity_path"],
        partition = get_partition('run_perplexity')
    shell:
        '''
        {params.perp_path} eval \
            --output {output} \
            --quants {input.train_quant} \
            --train_ecs {input.train_ecs} \
            --test_ecs {input.test_ecs} \
        '''

rule run_shoal_requant:
    input:
        expand(shoal_quant_ind, cond=[1,2], rep=REPLICATES, fold=FOLDS, w=shoal_params)
    params:
        partition = get_partition('run_shoal_requant')

rule _run_shoal_requant:
    input:
        prior=shoal_prior,
        quants_dir=train_sal_out_vb_dir
    output:
        shoal_quant_ind
    params:
        shoal_path=config["shoal_path"],
        partition = get_partition('_run_shoal_requant'),
        
    shell:
        '''
        {params.shoal_path} requant-v2\
            --prior {input.prior}\
            --sample {input.quants_dir}\
            --output {output}\
            --weight {wildcards.w}\
            -u 1\
            -f true
        '''

rule run_shoal_prior:
    input:expand(shoal_prior, fold=FOLDS)

rule _run_shoal_prior:
    input:
        lambda w:expand(train_sal_out, cond=[1,2], rep=REPLICATES, vb="1", fold=w.fold)
    output:
        shoal_prior
    params:
        partition = get_partition('quant'),
        shoal_path=config["shoal_path"]
    shell:
        '''
        {params.shoal_path} create-prior\
            --samples {input} \
            --output {output}
        '''

# rule _run_shoal_prior:

rule perplexity_salmon:
    input:expand(perp_file_sal, fold=FOLDS, vb=vbpriors, cond=[1,2], rep=REPLICATES)
    params:
        partition = get_partition('perplexity_salmon')

rule _perplexity_salmon:
    input:
        train_quant=train_sal_out,
        train_ecs=train_sal_ecs,
        test_ecs=test_sal_ecs
    output:
        perp_file=perp_file_sal
    params:
        perp_path=config["perplexity_path"],
        partition = get_partition('_perplexity_salmon')
    shell:
        '''
        {params.perp_path} eval\
            --output {output}\
            --quants {input.train_quant}\
            --train_ecs {input.train_ecs}\
            --test_ecs {input.test_ecs}\
        '''

rule run_salmon_test:
    input:
        expand(test_sal_ecs, rep = REPLICATES, cond = [1,2], fold = FOLDS)
    params:
        partition = get_partition('run_salmon_test')

rule _run_salmon_test:
    input:
        inp_fastq1 = test_fasta_l,
        inp_fastq2 = test_fasta_r,
        index = join(config['sal_ind_path'], "pos.bin")
    output:
        out_sal = test_sal_ecs
    resources: cpus=10, mem=32000
    params:
        out_dir = test_sal_out_dir,
        partition = get_partition('_run_salmon_test'),
        sal_path = config['sal_path'],
        sal_ind = config['sal_ind_path']
    shell:
        '''
        {params.sal_path} quant -i {params.sal_ind}\
             -l A -p {resources.cpus} --gcBias\
                -d\
                --skipQuant\
                --output {params.out_dir} -1 {input.inp_fastq1} -2 {input.inp_fastq2}    
        '''

rule run_salmon_train:
    input:
        expand(train_sal_out, rep = REPLICATES, cond = [1,2], fold = FOLDS, vb = vbpriors)
    params:
        partition = get_partition('run_salmon_train')

rule _run_salmon_train:
    input:
        inp_fastq1 = train_fasta_l,
        inp_fastq2 = train_fasta_r,
        index = join(config['sal_ind_path'], "pos.bin")
    output:
        out_sal = train_sal_out
    resources: cpus=10, mem=32000
    params:
        out_dir = train_sal_out_dir,
        partition = get_partition('run_salmon_train'),
        sal_path = config['sal_path'],
        sal_ind = config['sal_ind_path']
    shell:
        '''
        {params.sal_path} quant -i {params.sal_ind}\
             -l A -p {resources.cpus} --gcBias\
                -d\
                --vbPrior {wildcards.vb}\
                --output {params.out_dir} -1 {input.inp_fastq1} -2 {input.inp_fastq2}    
        '''


rule all_folds:
    input:
        expand(output_dir_fmt, rep=REPLICATES, cond=[1,2])
    params:
         partition = get_partition('all_folds')

rule kfold:
    input:
        l=L_reads_gz,
        r=R_reads_gz,
        bin=config['kfold-script']
    output:
        directory(output_dir_fmt)
    params:
        k=config['k'],
        seed=SEED,
        l_fa=L_reads_fa,
        r_fa=R_reads_fa,
        partition = get_partition('kfold')
    shell:
        '''
            echo "Uncompressing left read"
            zcat {input.l} > {params.l_fa}
            echo "Uncompressing right read"
            zcat {input.r} > {params.r_fa}

            {input.bin} \
            -l {params.l_fa} \
            -r {params.r_fa} \
            -k {params.k} \
            -s {params.seed} \
            -o {output}

            echo "Deleting deCompressed Stuff"
            rm {params.l_fa} {params.r_fa}
        '''