# Common imports
import pandas as pd
import os
from kNN import kNN

job = input('Enter job name: ')
os.mkdir(os.path.join(os.path.dirname(os.getcwd()), 'results', job))

# Load labes and sort the df
labels = pd.read_csv(os.path.join(os.path.dirname(os.getcwd()), 'data', 'clean_cluster_labels.csv'), sep=',')
labels['index_col'] = labels['Accession number']
labels.set_index('index_col', inplace=True)
# Load the file with query sequences
query_file = input('Path to query file: ')

X, y = labels[['Accession number', 'Description']], labels['Cluster numbers']

clf = kNN(k=10)

# Fit
clf.fit(X, y)

# Predict labels
predictions = clf.predict(X=query_file, job_name=job)
print(predictions)


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


