'''
Created on Oct 8, 2014

@author: agr9
'''

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel
import math
import csv
import os
import time
     

def main(protein_dict,q1,q2):
    if not protein_dict:
        angle_list = ['9999']
    elif len(protein_dict) < 2:
        angle_list = ['9999']
    else:
        train_set = [' '.join(protein_dict[x]) for x in protein_dict]
        proteins = [x for x in protein_dict]        
        tfidf_vectorizer = TfidfVectorizer()
        tfidf = tfidf_vectorizer.fit_transform(train_set)  #finds the tfidf score with normalization
#        print 'tfidf[0:1]', tfidf[0:1]
#        print 'tfidf[0:2]', tfidf[0:2]
#        print 'tfidf', tfidf
        cosine_similarities = linear_kernel(tfidf[0:1], tfidf).flatten()
#        print 'cosine_similarities', cosine_similarities
        related_docs_indices = cosine_similarities.argsort()[:-5:-1]

        degrees_list = []
        for a in (cosine_similarities[related_docs_indices].tolist()):

            try:
                angle_list = []
                angle_in_radians = math.acos(a)
                angle_in_degrees = math.degrees(angle_in_radians)
                degrees_list.append(angle_in_degrees) 
                angle_list.append(angle_in_degrees)   
            except ValueError:
                angle_list = ['9999']
                
        if len(degrees_list) > 1:
            return_list = [degrees_list[1]]
        else:
            return_list = degrees_list

        return return_list   
             
#        current_dir = os.getcwd()
#        write_path = current_dir + '/text_files/output_known_interactions.txt'
#        write_path = current_dir  + '/text_files/output_random_interactions.txt'
#        with open (write_path,'a') as f:
#            f.write(q1+'\t'+q2+'\t')
#            f.write(str(angle_list))
#            f.write(str('\t'))
#            f.write(str(protein_dict))
#            f.write(str('\n'))

    

if __name__=="__main__":
    MDM2= ['overexpression', 'expression', 'over', 'expression', 'inhibition', 'over', 'expression', 'association', 'association', 'interaction', 'binds', 'binds', 'interaction', 'affinity', 'affinity']
    TERT= ['expression', 'inducing', 'expression', 'detected', 'lacking', 'expression', 'inducing', 'expression', 'detected', 'lacking', 'expression']
    RANDOM= ['overexpression', 'adam', 'went', 'to', 'the', 'store', 'to', 'find', 'some', 'tomatoes', 'ended', 'up', 'buying', 'fish']
    train_set = [MDM2, TERT]
    protein_dict = {'MDM2':MDM2,'TERT':TERT}
    
    main(protein_dict, MDM2, TERT)
    