# Solubility-Weighted Index,
# Bhandari, B.K., Gardner, P.P. and Lim, C.S.,(2020),
# doi: 10.1093/bioinformatics/btaa578

import argparse
import re
import os
import sys

weights = {'A': 0.8356471476582918,
           'C': 0.5208088354857734,
           'E': 0.9876987431418378,
           'D': 0.9079044671339564,
           'G': 0.7997168496420723,
           'F': 0.5849790194237692,
           'I': 0.6784124413866582,
           'H': 0.8947913996466419,
           'K': 0.9267104557513497,
           'M': 0.6296623675420369,
           'L': 0.6554221515081433,
           'N': 0.8597433107431216,
           'Q': 0.789434648348208,
           'P': 0.8235328714705341,
           'S': 0.7440908318492778,
           'R': 0.7712466317693457,
           'T': 0.8096922697856334,
           'W': 0.6374678690957594,
           'V': 0.7357837119163659,
           'Y': 0.6112801822947587}

# Constants from logistic fitting
# prob = 1 / (1 + exp(-(a * x + b)));

A = 81.0581
B = -62.7775


def fasta_reader(file):
    '''Converts .fasta to a pandas dataframe with accession as index
    and sequence in a column 'sequence'
    '''
    print('Reading file...', end='\n')
    valid = re.compile('^[ACEDGFIHKMLNQPSRTWVY]+$')
    fasta_df = pd.read_csv(file, sep='>', lineterminator='>',
                           header=None, encoding='utf-8')
    fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', 1,
                                                                expand=True)
    fasta_df['Accession'] = fasta_df['Accession']
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True).\
        astype(str).str.upper().replace('U', 'C')
    fasta_df = fasta_df[fasta_df['Sequence'].apply(
        lambda x: True if valid.search(x) else False)]
    total_seq = fasta_df.shape[0]
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df = fasta_df[fasta_df.Sequence != '']
    fasta_df = fasta_df[fasta_df.Sequence != 'NONE']
    final_df = fasta_df.dropna()
    remained_seq = final_df.shape[0]
    if total_seq != remained_seq:
        print("{} sequences were removed due to inconsistencies in"
              "provided file.".format(total_seq-remained_seq))
    return final_df


def check_arg(args=None):
    '''arguments.
    '''
    parser = argparse.ArgumentParser(prog='SWI calculator',
                                     description='Solubility calculator using SWI. doi: 10.1093/bioinformatics/btaa578',
                                     epilog='(c) authors')
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s ' + '1',
                        help="Show program's version number and exit.")
    parser.add_argument('-f', '--fasta',
                        type=fasta_reader,
                        help='Input sequence (fasta)',
                        required=True)
    results = parser.parse_args(args)
    return (results.fasta)


def main():
    print('Computing Solubility-Weighted Index...', end='\n')
    df['SWI'] = df['Sequence'].apply(
        lambda x: np.mean([weights[i] for i in x]))
    print('Computing Probability of solubility...', end='\n')
    df['Prob. of Solubility'] = 1/(1 + np.exp(-(A*df['SWI'] + B)))
    print('Exporting file...', end='\n')
    output_fname = os.path.join(os.getcwd(), 'swi_results.csv')
    df.to_csv(output_fname, index=None)
    print('Done!')


if __name__ == '__main__':
    try:
        import pandas as pd
    except ImportError:
        print('\nPandas not found!\nInstall it first by:\n\npython3 -m pip install --user pandas\n')
        exit()
    try:
        import numpy as np
    except ImportError:
        print('\nNumpy not found!\nInstall it first by:\n\npython3 -m pip install --user numpy\n')
        exit()
    df = check_arg(sys.argv[1:])

    main()
