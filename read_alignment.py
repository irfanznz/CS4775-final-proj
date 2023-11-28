import argparse


def read_alignment(alignment_file, format):
    """
    Reads a file with specified format containing an alignment and returns a list of sequences
    """
    sequences = []
    if format == 'fasta':
        with open(alignment_file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    sequences.append('')
                else:
                    sequences[-1] += line.strip()
    if format == 'align':
        with open(alignment_file, 'r') as f:
            for line in f:
                sequences.append(line.strip())
    return sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads a file with an alignment and returns a list of sequences')
    parser.add_argument('alignment_file', type=str,
                        help='Path to the alignment file')
    parser.add_argument('format', type=str,
                        help='Format of the alignment file (fasta or align)')
    args = parser.parse_args()
    sequences = read_alignment(args.alignment_file, args.format)
    print(sequences)
