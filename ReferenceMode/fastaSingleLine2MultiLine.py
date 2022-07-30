import argparse

from Util import read_fasta


def print_seqs(header, sequence, length, outfile):
    print('>' + header, file=outfile)
    while len(sequence) > 0:
        print(sequence[:length], file=outfile)
        sequence = sequence[length:]

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run kmerRepFinder...')
    parser.add_argument('-f', metavar='fasta path',
                        help='input fasta path')

    args = parser.parse_args()

    singleLineFasta = args.f

    tmp_output_dir = '/home/hukang/KRF_output/dmel/CRD.2022-07-09.16-8-30'
    singleLineFasta = tmp_output_dir + '/test.fa'

    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-07-10.21-5-55'
    # singleLineFasta = tmp_output_dir + '/TE.merge.fa'
    contigNames, contigs = read_fasta(singleLineFasta)

    multiLineFasta = tmp_output_dir + '/test.ml.fa'

    outfile = open(multiLineFasta, 'w')  # open outfile for writing
    for name in contigNames:
        print_seqs(name, contigs[name], 50, outfile)