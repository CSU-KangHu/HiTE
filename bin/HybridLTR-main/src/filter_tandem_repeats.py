#!/usr/bin/env python3

import argparse
import os
import subprocess
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean up tandem repeat sequence and sequencing gap sequence in the fasta file")
    parser.add_argument('-f', required=True, help="Input fasta file")
    parser.add_argument('-misschar', default='n',
                        help="Define the letter representing unknown sequences; case insensitive; default: n")
    parser.add_argument('-Nscreen', type=int, choices=[0, 1], default=1,
                        help="Enable (1) or disable (0) the -nc parameter; default: 1")
    parser.add_argument('-nc', type=int, default=0,
                        help="Ambiguous sequence len cutoff; discard the entire sequence if > this number; default: 0")
    parser.add_argument('-nr', type=float, default=1,
                        help="Ambiguous sequence percentage cutoff; discard the entire sequence if > this number; default: 1")
    parser.add_argument('-minlen', type=int, default=100,
                        help="Minimum sequence length filter after clean up; default: 100 (bp)")
    parser.add_argument('-cleanN', type=int, choices=[0, 1], default=0,
                        help="Retain (0) or remove (1) the -misschar target in output sequence; default: 0")
    parser.add_argument('-trf', type=int, choices=[0, 1], default=1,
                        help="Enable (1) or disable (0) tandem repeat finder (trf); default: 1")
    parser.add_argument('-trf_path', default='', help="Path to the trf program")
    return parser.parse_args()


def run_trf(trf_path, file, align_score, max_seed):
    cmd = f"{trf_path} {file} 2 7 7 80 10 {align_score} {max_seed} -ngs -h -l 6"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout


def main():
    args = parse_args()

    if not os.path.isfile(args.f):
        print("Input fasta file not found!")
        exit(1)

    if not args.trf_path:
        result = subprocess.run("which trf", shell=True, capture_output=True, text=True)
        args.trf_path = result.stdout.strip()

    if not args.trf_path or not os.path.isfile(args.trf_path):
        print("Tandem Repeat Finder not found!")
        exit(1)

    tandem_sequences = {}
    if args.trf == 1:
        tandem_output = run_trf(args.trf_path, args.f, 1000, 2000)
        for line in tandem_output.split('\n'):
            if line.startswith('@'):
                seq_id = line[1:].strip()
                tandem_sequences[seq_id] = seq_id

    with open(f"{args.f}.cleanup", 'w') as info_file, open(args.f, 'r') as fasta_file:
        sequences = SeqIO.parse(fasta_file, 'fasta')
        for record in sequences:
            seq = str(record.seq)
            seq_id = record.id
            seq_len = len(seq)
            mark = False

            count = seq.lower().count(args.misschar.lower())
            count_rate = count / seq_len

            # Missing control
            if count_rate >= args.nr:
                info_file.write(f"{seq_id}\t{count_rate} missing\n")
                mark = True
            elif count >= args.nc and args.Nscreen == 1 and args.nc > 0:
                info_file.write(f"{seq_id}\t{count} sequence gap\n")
                mark = True

            # Tandem sequence control
            if seq_id in tandem_sequences:
                info_file.write(f"{seq_id}\ttandem sequence\n")
                mark = True

            # Remove missing seq and length control
            if args.cleanN == 1:
                seq = seq.replace(args.misschar.lower(), '').replace(args.misschar.upper(), '')
                new_len = len(seq)
                if new_len < args.minlen:
                    info_file.write(f"{seq_id}\tOnly {new_len} bp left after cleanup\n")
                    mark = True

            if not mark:
                print(f">{seq_id}\n{seq}")


if __name__ == "__main__":
    main()
