rule run_manta:
    input:
        normal=lambda w: f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"results/{w.run}/{w.sample}/bam/{w.sample}.bam",
        refg=config["refs"]["genome_human"],
        regions_bed=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz",
        regions_tbi=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz.tbi",
    output:
        workspace_dir=directory("work/manta/{run}/{sample}/workspace"),
        workflow="work/manta/{run}/{sample}/runWorkflow.py",
        workflow_pickle="work/manta/{run}/{sample}/runWorkflow.py.config.pickle",
        vcf_somaticSV="work/manta/{run}/{sample}/results/variants/somaticSV.vcf.gz",
        tbi_somaticSV="work/manta/{run}/{sample}/results/variants/somaticSV.vcf.gz.tbi",
    params:
        res_dir="work/manta/{run}/{sample}",
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["manta_max_gb"],
    conda:
        "../envs/manta.yaml"
    log:
        "work/logs/Manta_{run}_{sample}.log",
    shell:
        """
        configManta.py --exome \
        --referenceFasta {input.refg} \
        --normalBam {input.normal} \
        --tumorBam {input.tumor} \
        --runDir {params.res_dir} \
        --callRegions {input.regions_bed} \
        > {log} 2>&1

        {output.workflow} \
        --mode local \
        --jobs {resources.threads} --memGb {resources.mem_gb} \
        >> {log} 2>&1
        """


rule filter_manta_variants:
    input:
        vcf="work/manta/{run}/{sample}/results/variants/somaticSV.vcf.gz",
        tbi="work/manta/{run}/{sample}/results/variants/somaticSV.vcf.gz.tbi",
    output:
        vcf="work/manta/{run}/{sample}/{sample}.SV.filtered.vcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} | bcftools sort -Ov -o {output.vcf}
        """


rule annotate_sv:
    input:
        vcf="work/manta/{run}/{sample}/{sample}.SV.filtered.vcf",
    output:
        tsv="work/manta/{run}/{sample}/{sample}.SV.annotated.tsv",
    params:
        annotations=config["refs"]["annotsv_annotations"]["path"],
        genome_build=config["refs"]["annotsv_annotations"]["genome_build"],
        output_prefix="work/manta/{run}/{sample}/{sample}.SV.annotated"
    threads: 1
    conda:
        "../envs/annotsv.yaml"
    log:
        "work/logs/AnnotSV_{run}_{sample}.log",
    shell:
        """
        variant_count=$(bcftools view -H {input.vcf} | wc -l)
        
        if [ "$variant_count" -eq 0 ]; then
            touch {output.tsv}
            echo "No variants to annotate" > {log}
        else
            AnnotSV \
            -SVinputFile {input.vcf} \
            -outputFile {params.output_prefix} \
            -genomeBuild {params.genome_build} \
            -annotationsDir {params.annotations} \
            > {log} 2>&1
        fi
        """