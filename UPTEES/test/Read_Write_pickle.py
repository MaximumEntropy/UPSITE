'''
Created on Jan 29, 2015

@author: adam
'''
import pickle
import os




# write
#with open(the_filename, 'wb') as f:
#    pickle.dump(my_list, f)
    
#read
#with open(the_filename, 'rb') as f:
#    my_list = pickle.load(f)
    
    
with open('/home/adam/workspace/TEES/text_files/Author_Lists/output/SOM_Faculty.tsv_tfidf', 'rb') as f2:
    q = pickle.load(f2)
    print 'q', q