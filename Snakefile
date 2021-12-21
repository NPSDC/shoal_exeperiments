from os.path import join
import yaml

### Try to run perplexity - done
### See the output on a single sample
### from here try to run shoal on multiple samples and see the output
configfile:
    "config.json"

comm_path=join("aux_info", "eq_classes.txt.gz")
K=[1,2,3,4,5]
shoal_params = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]
condition = ["A", "B"]
replicates = [1,2,3,4]

shoal_prior = join(config["shoal_out"], "vbprior=1e0", "A_B", "{k}","prior.tsv")
seqc_old_quant_dir = join(config["seqc_quant"], "BGI_FC1_{cond}_{rep}", "vbprior=1e0", "{k}")
seqc_old_quant_ind = join(seqc_old_quant_dir,"quant.sf")
shoal_quant_ind=join(config["shoal_out"], "vbprior=1e0", "A_B", "{k}", "BGI_FC1_{cond}_{rep}_c={w}_adapt.sf")
perp_file=join(config["perp_out_path"], "vbprior=1e0", "c={w}", 'BGI_FC1_{cond}_{rep}', '{k}', 'perplexity.yml')
train_ecs=join(config['seqc_quant'], 'BGI_FC1_{cond}_{rep}', "vbprior=1e0", '{k}', comm_path)
test_ecs=join(config['seqc_test'], 'BGI_FC1_{cond}_{rep}', '{k}', comm_path)
gene_quant_ind=join(config['shoal_gene_out'], "vbprior=1e0", "A_B", "{k}", 'BGI_FC1_{cond}_{rep}_c={w}_gene_quant.sf')
sal_gene_quant_ind=join(config["sal_gene_out"], "BGI_FC1_{cond}_{rep}", "vbprior=1e0", "{k}", "gene_quant.sf")

with open('cluster_partition.conf', 'r') as f:
    partitions = yaml.safe_load(f)

def get_partition(rule_name):
    return partitions[rule_name] if rule_name in partitions else partitions['__default__']

rule salgene:
    input:
        expand(sal_gene_quant_ind, cond = condition, rep = replicates, k = K)
    params:
        partition = get_partition('salgene')

rule _salgene:
    input:
        seqc_old_quant_ind
    output:
        sal_gene_quant_ind
    params:
        partition = get_partition('_salgene')
    shell:
        '''
            source ~/.bashrc
            conda activate R4.1
            Rscript tx2gene.R \
                {input} \
                {output}
        '''

rule tx2gene:
    input:
        expand(gene_quant_ind, cond = condition, rep = replicates, w = [10,100,1000], k = K)
    params:
        partition = get_partition('tx2gene')

rule _tx2gene:
    input:
        shoal_quant_ind
    output:
        gene_quant_ind
    params:
        partition = get_partition('_tx2gene')
    shell:
        '''
            source ~/.bashrc
            conda activate R4.1
            Rscript tx2gene.R \
                {input} \
                {output}
        '''


rule perplexity:
    input:expand(perp_file, k=K, w=shoal_params, cond=condition,rep=replicates)
    params:
        partition = get_partition('run_perplexity')

rule run_perplexity:
    input:
        train_quant=shoal_quant_ind,
        train_ecs=train_ecs,
        test_ecs=test_ecs
    output:
        perp_file=perp_file
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



rule shoal_requant:
    input:
        expand(shoal_quant_ind, cond=condition, rep=replicates,k=K, w=shoal_params)
    params:
        partition = get_partition('shoal_requant')

rule run_shoal_requant:
    input:
        prior=shoal_prior,
        quants_dir=seqc_old_quant_dir,
    output:
        shoal_quant_ind
    params:
        shoal_path=config["shoal_path"],
        partition = get_partition('run_shoal_requant'),
        

    shell:
        '''
        {params.shoal_path} requant\
            --prior {input.prior}\
            --sample {input.quants_dir}\
            --output {output}\
            --weight {wildcards.w}
        '''

rule shoal_prior:
    input:
        expand(shoal_prior, k=K)
    params:
        partition = get_partition('shoal_prior')

rule quant:
    input:
        expand(seqc_old_quant_ind, cond=condition, rep=replicates,k=K)
    params:
        partition = get_partition('quant')

rule run_shoal_prior:
    input:
        lambda w:expand(seqc_old_quant_ind, cond=condition, rep=replicates,k=w.k)
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


        


