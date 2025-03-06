mkdir -p refs/xengsort

xengsort index --index refs/xengsort/xengsort-k25 \
            -G "/media/data/NGS/refs/broad/Homo_sapiens_assembly38.fasta" \
            -H "/media/data/NGS/refs/mm39/mm39.fa" \
            -k 25 -n 4500000000 \
            --fill 0.88 \
            -W 16