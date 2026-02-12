FROM nanopore_pipe:v1
WORKDIR /work
COPY . /work

# python pipeline.py -i in/basecall_out.bam  -b in/barcode_reference_sample_2024-10-21.xlsx  -w in/well_reference_sample_2024-10-21.xlsx  -o out/
