import os
from sbc_ngs.pathway import PathwayAligner, parse_barcode_excel
from sbc_ngs.process_result import exact_clone_well_info
import subprocess
import argparse
import logging
import sys
import pandas as pd
from sva2 import NanoporeMappingPipeline
from sva2.validate import MutationDetector
import json

parser = argparse.ArgumentParser(description='Pipeline for Nanopore sequencing data')
parser.add_argument('-i', '--input', help='Input bam file', required=True)
parser.add_argument('-b', '--barcode', help='Barcode excel file', required=True)
parser.add_argument('-w', '--well', help='Well info excel file', required=True)
parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', default=64, type=int)
parser.add_argument('--quiet', help='Quiet mode', action='store_true')

nt_base = set('ATCG')

if __name__=='__main__':
    args = parser.parse_args()
    input_bam = args.input
    out_dir = args.out_dir
    threads = int(args.threads)
    quiet = args.quiet
    in_dir = os.path.dirname(input_bam)
    fastq = input_bam.replace('.bam', '.fastq')
    fastq_dir = os.path.join(in_dir, 'fastq')

    if not os.path.exists(fastq_dir): 
        os.mkdir(fastq_dir)
        subprocess.check_call(f"samtools fastq {input_bam} > {fastq}", shell=True)
        subprocess.check_call(f"split -l 20000 {fastq} {os.path.join(fastq_dir, 'chunk_')} --additional-suffix=.fastq", shell=True)
        os.remove(fastq)
    # os.remove(input_bam)

    barcodes_df = parse_barcode_excel(args.barcode)
    demultiplexed_dir = os.path.join(out_dir, 'demultiplexed')
    sva_result_dir = os.path.join(out_dir, 'sva_result')

    if not os.path.exists(sva_result_dir):
        os.mkdir(sva_result_dir)

    df = exact_clone_well_info(args.well).merge(barcodes_df, on='well')
    ref_df = pd.read_excel(args.barcode, sheet_name='reference')
    ref_df.sequence = ref_df.sequence.str.upper().apply(lambda x: ''.join([c for c in x if c in nt_base]))

    refs = ref_df.set_index('known_seq_id')['sequence'].str.upper().to_dict()
    other_refs = list(refs.values())

    if not os.path.exists(demultiplexed_dir) or len(os.listdir(demultiplexed_dir)) < barcodes_df.shape[0] * 0.9 or len(os.listdir(sva_result_dir)>=len(df)):
        print('Demultiplexing...')
        pathway = PathwayAligner(out_dir=demultiplexed_dir, in_dir=fastq_dir, barcodes_df=barcodes_df)
        pathway.demultiplex(tolerance=16, num_threads=threads, indiv_strand=False)

    # for i, row in df.iterrows():
    #     dir_ = f"{row['forward']}_{row['reverse']}"
    #     dir_ = os.path.join(demultiplexed_dir, dir_)
        
    #     if os.path.exists(dir_):
    #         dest = os.path.join(demultiplexed_dir, row['clone_id'])
    #         subprocess.check_call(f"mv {dir_} {dest}", shell=True)
    

    logger = logging.getLogger('sva')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.FileHandler(os.path.join(sva_result_dir, 'sva.log')))
    formater = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger.handlers[0].setFormatter(formater)

    if not quiet:
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.handlers[1].setFormatter(formater)

    mutation_detector = MutationDetector(model_path="model.pth", input_dim=25, embedding_dim=96)

    summary = []

    for i, row in df.iterrows():
        seq_id = str(row['known_seq_id'])
        clone_id = row['clone_id']

        clone_fastq = os.path.join(demultiplexed_dir, f"{row['forward']}_{row['reverse']}", 'sum_reads.fastq')
        reference = refs[seq_id]
        prefix = os.path.join(sva_result_dir, clone_id, 'sva')

        if os.path.exists(prefix+'_mutated.csv'):
            result = {"contig": pd.read_csv(prefix+'_summary.csv'), 
                      "mutations": pd.read_csv(prefix+'_mutated.csv') if os.path.getsize(prefix+'_mutated.csv') > 10 else pd.DataFrame(),
                      "plus_contig": pd.read_csv(prefix+'_plus_contig.csv'),
                      "insertion": pd.read_csv(prefix+'_insertion.csv') if os.path.getsize(prefix+'_insertion.csv') > 10 else pd.DataFrame(),
                      "status": json.load(open(prefix+'_status.json'))}
        elif not os.path.exists(clone_fastq):
            row['variants'] = []
            row['DP'] = []
            row['Risk'] = []
            summary.append(row)
            print(f"Clone {clone_id} does not exist.")
            continue

        else:
            mp = NanoporeMappingPipeline(clone_fastq, reference, detect_threshold=0.1,
                                         other_refs=other_refs,
                                         threads=threads, logger=logger, prefix=prefix, Q_threshold=10)
            
            result = mp.run_multiprocess()

        depth = result['contig'][[base + '_Count' for base in 'ACGT']].sum(axis=1)
        mutations = result['mutations']
            
        
        max_depth = depth.max()
        min_depth = depth.min()
        mean_depth = depth.mean()
        identity = (depth >= max(mean_depth * 0.05, 3)).mean()
        match_rate = result['status']['mapped_reads'] / (result['status']['valid_reads'] + 1e-6)
        valid_reads = result['status']['valid_reads']

        mutation_detector.init_with_result(result)

        XF = []
        DP = []
        Risk = []
        
        if len(mutations) > 0:
            i = 0

            max_mutation_ratio = mutations['Ratio'].max()
            
            continous_deletion = lambda x, y: x['Mutation'] == 'Gap' and y['Mutation'] == 'Gap' and x['Position'] + 1 == y['Position']
            while i < len(mutations):
                mutation = mutations.iloc[i].to_dict()

                # Continous deletion
                if i < len(mutations) - 1 and continous_deletion(mutation, mutations.iloc[i+1]):
                    start = i
                    while i < len(mutations) - 1 and continous_deletion(mutations.iloc[i], mutations.iloc[i+1]):
                        i += 1
                        mutation['Reference'] += mutations.iloc[i]['Reference']
                        mutation['Ratio'] += mutations.iloc[i]['Ratio']
                        mutation['Plus_Ratio'] += mutations.iloc[i]['Plus_Ratio']

                    mutation['Ratio'] /= i - start + 1
                    mutation['Plus_Ratio'] /= i - start + 1

                # homopolymer deletion
                if mutation['Mutation'] == 'Gap' and len(mutation['Reference']) == 1:
                    front = reference[:mutation['Position']-1][::-1]
                    back = reference[mutation['Position']:]
                    homo = 1
                    for base in front:
                        if base == mutation['Reference']:
                            homo += 1
                        else:
                            break
                    for base in back:
                        if base == mutation['Reference']:
                            homo += 1
                        else:
                            break
                    if homo >= 3:
                        mutation['Mutation'] += f"(homo:{homo})"
                        
                i += 1
                if mutation['Ratio'] > 0.85:
                    XF.append(mutation)
                else: # >0.5 -> DP: 2026-02-10
                    if mutation['Ratio'] > 0.5 or len(mutation['Reference'])>2 or (mutation_detector(mutation) and abs(mutation['Plus_Ratio']-0.5)<0.45):
                        DP.append(mutation)
                    else:
                        Risk.append(mutation)

        row['identity'] = identity
        row['max_depth'] = max_depth
        row['mean_depth'] = mean_depth
        row['match_rate'] = match_rate
        row['valid_reads'] = valid_reads
        row['variants'] = XF
        row['DP'] = DP
        row['Risk'] = Risk

        summary.append(row)
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(sva_result_dir, 'summary.csv'), index=False)

    pass
