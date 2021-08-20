import pandas as pd
import os
from collections import Counter
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.ExPASy import ScanProsite

database = os.path.join(os.path.dirname(os.getcwd()), 'database', 'eukaryota_db')

def blast(sequence, database, out):

    """ Compares query sequence to a sequence database using BLAST algorithm.

        Returns a dataframe of hits and pairwise alignment scores. """

    hits = {}

    # Run BLAST
    blastp = NcbiblastpCommandline(cmd='blastp', outfmt=5, query=sequence, db=database, evalue=0.001, out=out)
    stdout, stderr = blastp()

    # Read BLAST results
    blast_result_handle = open(out)
    blast_record = NCBIXML.read(blast_result_handle)

    # Parse BLAST results
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001:
                hits[alignment.title.split('|')[3]] = hsp.score

    return pd.DataFrame.from_dict(hits, columns= ['score'], orient='index')


class kNN:

    def __init__(self, k=10):
        self.k = k


    def fit(self, X, y):
        self.X_train = X
        self.y_train = y


    def predict(self, X, output):

        print('prediciting')
        print(output)
        pred_labels = [self._predict(x=x, output=output) for x in SeqIO.parse(X, 'fasta')]

        return pred_labels


    def _predict(self, x, output):

        # Create a temporary query seqeunce file
        tmp_seq = os.path.join(os.path.dirname(os.getcwd()), 'tmp_seq.txt')
        with open(tmp_seq, 'w') as f:
            f.write('>' + str(x.name).split('|')[1] + '\n')
            f.write(str(x.seq) + '\n')

        # Before running BLAST check if the sequence contains the GNAT fold (by scanning PROSITE)
        if self._scan_prosite(x) == 'GNAT':
            # Calculate similarity between sequences using BLAST
            #output = os.path.join(os.path.dirname(os.getcwd()), 'results', job_name, str(x.name.split('|')[1]) + '.xml')
            print('Running BLAST')

            blast_output = os.path.join(output, 'BLAST_results', str(x.name.split('|')[1]) + '.xml')

            scores = blast(sequence=tmp_seq, database=database, out=blast_output)
            print(f'BLAST results for {str(x.name.split("|")[1])} written to {blast_output}')

            # Get labels of k nearest neighbors
            print('Counting k neighbours')
            k_indices = scores.sort_values(ascending=False, axis=0, by=['score'])[:self.k]
            knn_labels = k_indices.merge(self.y_train, left_index=True, right_index=True, how='inner')

            # Majority vote
            majority = Counter(knn_labels['Cluster numbers'].to_list()).most_common(1)
            print(majority[0][0])

            os.remove(tmp_seq)

            return (str(x.name).split('|')[1], majority[0][0])

        else:

            os.remove(tmp_seq)

            print(f'Sequence {str(x.name).split("|")[1]} does not contain the GNAT fold.')


    def _scan_prosite(self, x):

        """ Scans prosite and checks whether the query sequence contains the GNAT fold """

        # Possible prosite signatures
        prosite_gnat_codes = ['PS51186', 'PS51730', 'PS51731', 'PS51729']

        # Scan Prosite
        handle = ScanProsite.scan(x.seq)

        # Read the results
        results = ScanProsite.read(handle)

        print(f'Checking whether {str(x.name).split("|")[1]} contains the GNAT fold')
        c = 0
        for result in results:
            if result['signature_ac'] in prosite_gnat_codes:
                c+=1

        if c > 0:
            return 'GNAT'
        else:
            return 'not GNAT'

