'''
Created on Jan 28, 2015

@author: adam
'''

import os
import csv

file_location = '/home/adam/workspace/TEES/text_files/Author_Lists/SOM_Faculty_short.tsv'  #School of Medicine Faculty, cut short for testing
base_name = os.path.basename(file_location)
output_file = '/home/adam/workspace/TEES/text_files/Author_Lists/output/{0}' .format(base_name)

Author_frequency_dict = {'Yu Yanping': [('LSP1 gene', 0.2), ('YIF1A', 0.2), ('TBC1D22A', 0.2), ('MAML2', 0.2), ('RP1-1777G6', 0.2)], 'Simmons William': [], 'Jeyabalan Geetha': [('HMGB1', 0.222222), ('TLR4', 0.122222), ('IRF-1', 0.1), ('TNF-alpha', 0.044444), ('alpha-GalCer', 0.033333), ('A2AR', 0.033333), ('JNK', 0.033333), ('IL-13 neutralizing antibody', 0.022222), ('nitric oxide synthase', 0.022222), ('IL-6', 0.022222), ('IRF-1 mRNA', 0.022222), ('iNOS', 0.022222), ('IL-13', 0.022222), ('IFN-gamma', 0.022222), ('SH58261 diminished alpha', 0.011111), ('IL-1R', 0.011111), ('MAP kinase c-Jun NH(2)-terminal kinase', 0.011111), ('IRF-1 gene', 0.011111), ('arginase', 0.011111), ('calmodulin', 0.011111), ('High-mobility group box 1', 0.011111), ('Interferon (IFN) regulatory factor-1', 0.011111), ('ROS', 0.011111), ('glycolipid antigen alpha-galactosylceramide', 0.011111), ('nuclear protein high-mobility group box 1', 0.011111), ('adenosine A2AR', 0.011111), ('CD1d', 0.011111), ('CaMK', 0.011111), ('adenosine A2A receptor', 0.011111), ('High mobility group box 1', 0.011111), ('IFN-gamma and IFN-beta', 0.011111), ('IRF-1 protein', 0.011111), ('GalCer', 0.011111), ('IL-1beta', 0.011111), ('hepatocellular IRF-1', 0.011111), ('ICAM-1', 0.011111), ('NOS', 0.011111)]}

# with open(output_file, 'wb') as f:  # Just use 'w' mode in 3.x
#     w = csv.writer(f)
#     w.writerows(Author_frequency_dict.items())

# with open(output_file, 'wb') as f:  # Just use 'w' mode in 3.x
#     w = csv.DictWriter(f, Author_frequency_dict.keys())
#     w.writeheader()
#     w.writerow(Author_frequency_dict)
 
def saveDict(fn,dict_rap):
    f=open(fn, "wb")
    w = csv.writer(f)
    for key, val in dict_rap.items():
        w.writerow([key, val])
    f.close()
     
def readDict(fn):
    f=open(fn,'rb')
    dict_rap={}
    for key, val in csv.reader(f):
        dict_rap[key]=eval(val)
    f.close()
    return(dict_rap)

saveDict(output_file, Author_frequency_dict)
a = readDict(output_file)
print a

