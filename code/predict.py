# Common imports
import pandas as pd
import os
from kNN import kNN
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Specify the file with query sequences.")
    parser.add_argument("output", help="Specify the output file")
    parser.add_argument("-k", help='Select k-neighbours to be considered', type=int)
    parser.add_argument("-e", "--e_value", help="Specify E-value for BLAST.", type=int)
    #parser.add_argument("-m", "--metric", help="Specify the measure of similarity between proteins. "
    #                                          "The default is alignment score.")
    #parser.add_argument("-l", "--labels", help="Select file with protein labels.")
    #parser.add_argument("-p", "--prosite", help="Provide a list of PROSITE signatures your group of proteins" \
    #                                           "should contain")
    parser.add_argument("-s", "--score", help="Specify minimum alignment score to be considered for "
                                              "evaluating neighborhood")

    args = parser.parse_args()

    #job = input('Enter job name: ')
    if not os.path.exists(args.output):
        os.makedirs(os.path.join(args.output, 'BLAST_results'))
    else:
        replace_dir = input('This folder already exists. Overwrite the existing folder? [y,n] ')
        if replace_dir == 'y'.lower():
            pass
        elif replace_dir == 'n'.lower():
            pass


    # Load labes and sort the df
    labels = pd.read_csv(os.path.join(os.path.dirname(os.getcwd()), 'data', 'clean_cluster_labels.csv'), sep=',')
    labels['index_col'] = labels['Accession number']
    labels.set_index('index_col', inplace=True)

    X, y = labels[['Accession number', 'Description']], labels['Cluster numbers']

    if args.k:
        clf = kNN(k=args.k)
    else:
        clf = kNN(k=10)

    # Fit
    clf.fit(X, y)

    # Predict labels
    if args.e_value:
        predictions = clf.predict(X=args.input, output=args.output, e_value=args.e_value)
        print(predictions)
    else:
        predictions = clf.predict(X=args.input, output=args.output, e_value=0.001)
        print(predictions)


if __name__ == '__main__':
    main()



# Test accuracy of the kNN model
# query_file = ''

# from sklearn.model_selection import train_test_split

# Split the set
#X_train, X_test, y_train, y_test = train_test_split(labels[['Accession number', 'Description']],
#                                                    labels['Cluster numbers'], test_size=0.1, random_state=123)

#clf = kNN(k=10)

# Fit
#clf.fit(X=X_train, y=y_train)

# Predict
#predictions = clf.predict(X=query_file)

# Accuracy of predictions
#s = 0
#for prediction in predictions:
#    for index, item in y_test.iteritems():
#       if prediction[0] == index and prediction[1] == item:
#            s+=1

#acc = s / len(y_test.to_list())
#print(acc)


