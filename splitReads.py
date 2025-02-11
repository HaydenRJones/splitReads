# HTbfRS
# Hacked Together (but/barely functional) read splitter

# findMotif.py
# HJ v0.1 - 28/10/24
# v0.2 - 25/12/24
# - file writing signficantly improved

# Find, plot distribution, and split reads on specific motifs in some sequence file
# Useful for identifying and fixing chimeric sequences (maybe!)

# Imports

from Bio import SeqIO
from itertools import chain

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import regex
import argparse

# Funcs

def loadSeqs(file):
    
    seqList = []
    temp = list(SeqIO.parse(f'{file}', 'fastq'))
    for record in temp:
        seqList.append([record.id, 
                        str(record.seq), 
                        ''.join([chr(score + 33) for score in record.letter_annotations['phred_quality']]), 
                        len(record.seq)])
    
    df = pd.DataFrame(seqList, columns = ['ID', 'Seq', 'Qual', 'Length'])
    
    print(f'loaded {file} : {len(df)} seqs')
    
    return(df)#, seqList)


def findMatches(df, pattern):
    
    df[f'{pattern}_absolute_matches'] = [[] for _ in range(len(df))]
    df[f'{pattern}_relative_matches'] = [[] for _ in range(len(df))]
    df[f'{pattern}_count']            = 0
    
    for i, seq in df.iterrows():
        
        absolutePosList = []
        relativePosList = []
        matches = regex.finditer(pattern, df.loc[i, 'Seq'])
        #matches = pattern.finditer(df.loc[i, 'Seq']) OLD METHOD USING RE
        
        for match in matches:
            
            absolutePosList.append(match.start())
            relativePosList.append(match.start() / seq['Length'])

        df.at[i, f'{pattern}_absolute_matches'] = absolutePosList        
        df.at[i, f'{pattern}_relative_matches'] = relativePosList
        df.at[i, f'{pattern}_count']            = len(relativePosList)

    print('found matches...')
    
    return()


def makePlot(df):
    
    plt.figure()
    plot = sns.boxplot(data = df, 
                x = f'{pattern}_count',
                y = 'Length',
                log_scale=True)

    plot.set(xlabel = 'Number of regex matches',
          title = f'{pattern}')

    temp_matches = list(chain.from_iterable(df[f'{pattern}_relative_matches']))
    matches_df = pd.DataFrame(temp_matches, columns = ['matches'])

    plt.figure()
    plot = sns.histplot(data = matches_df,
                 x = 'matches')
    plot.set(xlabel = 'Relative position in read',
          title = f'{pattern}')
    plot.set_xlim(-0.01, 1.01)
    
    return()


def splitRead(seq, pattern):
    
    regions = [0] + seq.loc[f'{pattern}_absolute_matches'] + [int(seq.loc['Length'])]
    
    newSequences = []
    
    for i in range(0, seq[f'{pattern}_count'] + 1):
        
        newName = f'''{seq.loc['ID']}_{i}'''
        span    = [regions[i], regions[i + 1]]
        
        newSeq  = seq.loc['Seq'][span[0]:span[1]]
        newQual = seq.loc['Qual'][span[0]:span[1]]
        
        newSequences.append(f'@{newName}\n{newSeq}\n+\n{newQual}\n')
    
    return(newSequences)


# Main loop
    
if __name__ == '__main__':
   
    parser = argparse.ArgumentParser(
        prog = 'findMotif.py',
        description = 'Find and plot distribution of specific motifs in some sequence file')
    
    parser.add_argument('file', help = 'input file\nAssumes .fastq')
    args = parser.parse_args()
    
    
    # Fixed patterns for now
    patterns = [['5primeFlank',     'GGTGCTG'],
                ['5primeFlank_rev', 'CAGCACC']]#,
                #['bc01', 'AAGAAAGTTGTCGGTGTCTTTGTG'],
                #['bc02', 'TCGATTCCGTTTGTAGTCGTCTGT'],
                #['bc03', 'GAGTCTTGTGTCCCAGTTACCAGG'],
                #['bc04', 'TTCGGATTCTATCGTGTTTCCCTA'],
                #['bc05', 'CTTGTCCAGGGTTTGTGTAACCTT'],
                #['bc06', 'TTCTCGCAAAGGCAGAAAGTAGTC'],
                #['bc07', 'GTGTTACCGTGGGAATGAATCCTT'],
                #['bc08', 'TTCAGGGAACAAACCAAGTTACGT'],
                #['bc09', 'AACTAGGCACAGCGAGTCTTGGTT'],
                #['bc10', 'AAGCGTTGAAACCTTTGTCCTCTC'],
                #['bc11', 'GTTTCATCTATCGGAGGGAATGGA'],
                #['bc12', 'CAGGTAGAAAGAAGCAGAATCGGA'],
                #['bc01_rev', 'CACAAAGACACCGACAACTTTCTT'],
                #['bc02_rev', 'ACAGACGACTACAAACGGAATCGA'],
                #['bc03_rev', 'CCTGGTAACTGGGACACAAGACTC'],
                #['bc04_rev', 'TAGGGAAACACGATAGAATCCGAA'],
                #['bc05_rev', 'AAGGTTACACAAACCCTGGACAAG'],
                #['bc06_rev', 'GACTACTTTCTGCCTTTGCGAGAA'],
                #['bc07_rev', 'AAGGATTCATTCCCACGGTAACAC'],
                #['bc08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA'],
                #['bc09_rev', 'AACCAAGACTCGCTGTGCCTAGTT'],
                #['bc10_rev', 'GAGAGGACAAAGGTTTCAACGCTT'],
                #['bc11_rev', 'TCCATTCCCTCCGATAGATGAAAC'],
                #['bc12_rev', 'TCCGATTCTGCTTCTTTCTACCTG']]  
    
    pattern = f'({patterns[0][1]}){{e<=0}}' # Error stuff is something that needs to be explored more
    
    
    # Load sequences and find matches
    
    if args.file == '': 
        print('NO FILE!')
        quit()
    
    df = loadSeqs(f'{args.file}.fastq') #, seqList
    findMatches(df, pattern)
    
    #makePlot(df)
    
    # Write sequences
    newSeqs = []
    for i in range(0, len(df)):
        newSeqs.append(splitRead(df.loc[i], pattern))  
    print('split sequences...')
    
    with open(f'{args.file}_split.fastq', 'a') as f:
        f.writelines(list(chain.from_iterable(newSeqs)))
    
    print(f'Wrote file : ./{args.file}_split.fastq')
    
    print(f'''Results:
          --------------------------------
          Input Reads        : {len(df)}
          Average Length     : {np.int16(df['Length'].mean())}
          New Total Reads    : {len(list(chain.from_iterable(newSeqs)))}
          New Average Length : {np.int16(sum(df['Length']) / len(list(chain.from_iterable(newSeqs))))}
          --------------------------------''')
    
    #{int(statistics.mean(list(chain.from_iterable(newSeqLengths))))}\n
