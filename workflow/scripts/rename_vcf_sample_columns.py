def rename_samples(input_vcf, output_vcf, normal_name, tumor_name):
    with open(input_vcf) as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                # Find indexes of NORMAL and TUMOR columns
                try:
                    normal_idx = fields.index('NORMAL')
                    tumor_idx = fields.index('TUMOR')
                    # Replace sample names
                    fields[normal_idx] = normal_name
                    fields[tumor_idx] = tumor_name
                    line = '\t'.join(fields) + '\n'
                except ValueError as e:
                    raise ValueError("Could not find NORMAL and/or TUMOR columns in VCF") from e
            outfile.write(line)


if __name__ == "__main__":
    rename_samples(snakemake.input[0],
                  snakemake.output[0],
                  f"{snakemake.params.normal}",
                  f"{snakemake.wildcards.sample}")
