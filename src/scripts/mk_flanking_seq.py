import argparse, sys


def mk_mat(args, seqs):
    with open(args.matrix, 'r') as inf, open(args.outFile, 'w') as of:
        header = inf.readline().strip() + '\t' + 'flanking sequence'
        print(header, file=of)
        for line in inf:
            pos = line.strip().split('\t')[-2]
            new_line = line.strip() + '\t' + seqs.get(pos)
            print(new_line, file=of)

def get_seqs(args):
    seqs = {}
    with open(args.faFile, 'r') as inf:
        while True:
            line = inf.readline().strip()
            print(line)
            if line == '':
                break
            pos = line.split(':')[1].split('-')[0]
            print(pos)
            seq_line = inf.readline().strip()        
            print(seq_line)
            seqs[pos] = seq_line
    return seqs

def main(args):
    seqs = get_seqs(args)
    mk_mat(args,seqs)
    

if __name__ == "__main__":
    desc = 'Add flanking regions from fa to matrix.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('faFile', 'matrix', 'outFile')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
