'''
Created on Oct 27, 2014

@author: Adam
'''

import BANNER2
import operator
import sys
from sklearn.feature_extraction.text import TfidfVectorizer


def get_author_entity_dict(docs_dict, Author, current_author_num, tot_authors):
    number_docs = len(docs_dict)
    author_entity_dict = {}
    tfidf_lists = []
    for num, key in enumerate(docs_dict):
        print Author, "----------------------- Author {0}/{1}....Publication {2}/{3} -----------------------" .format(current_author_num, tot_authors, num + 1, number_docs)
        try:
            a = docs_dict[key][0]
        except IndexError:
            a = []
        if not isinstance(a, str):
            if not isinstance(a, unicode):
                pass
        else:
#   a = docs_dict['23231918'][0]
            b = a.encode('ascii','ignore')
            c = b.split('.')
            sent_list_final= []
            for sent in c:
                if len(sent) > 1:
                    sent_list_final.append(sent)
            entity_list = BANNER2.main(sent_list_final)
            tfidf_lists.append(entity_list)
            for x in entity_list:
                if x in author_entity_dict:
                    author_entity_dict[x] += 1
                if x not in author_entity_dict:
                    author_entity_dict[x] = 1
            
        
#   print 'author_entity_dict'
#   print author_entity_dict
#    author_entity_dict = {'TRPML1 KD': 4, 'Ca(2+) release': 1, 'EBP50': 1, 'PLCgamma': 1, 'lysosomal ion channel TRPML1': 1, 'transient receptor potential': 1, 'TRPC2': 9, 'transcription factor MTF-1': 1, 'ML1': 10, 'Gb3': 3, 'TRP': 16, 'Toc': 1, 'alpha-galactosidase': 1, 'proapoptotic protein Bax': 1, 'trp2 mutant': 1, 'phospholipase C': 1, 'CatB': 3, 'lysosomal ion homeostasis': 1, 'calmodulin': 1, 'IP(3) receptors': 1, 'LRRK2': 4, 'TRPML3': 1, 'lysosomal SNARE proteins': 1, 'mitochondrial Ca2+': 1, 'GPCR': 3, 'plasma membrane receptors': 1, 'TRPML1': 28, 'zinc transporter ZnT4': 1, 'InaD': 1, 'ROS': 1, 'VAMP7 KD': 1, 'apolipoprotein B hydrolysis in MLIV': 1, 'Ca(2+) release channels': 1, 'caveolin': 1, 'GPI': 1, 'lysosomal enzymes': 1, 'TRPC': 15, 'leucine-rich repeat kinase 2': 1, 'LRRK2 G2019S': 1, 'tyrosine kinase receptors': 1, 'caspase': 1, 'NEHRF': 1, 'TRP family': 3, 'TRP2': 1, 'TRPML1in zinc transport': 1, 'MCOLN1': 1, 'G protein-coupled receptors': 1, 'R1441C': 1, 'scaffolding proteins': 1, 'lysosomal protease cathepsin B': 1, 'TRPML1in': 1, 'RPE1': 1, 'G protein coupled receptors': 1, 'synaptotagmin VII': 1, 'ROS chelator': 1, 'VAMP7': 1, 'transient receptor potential mucolipin 1': 2, 'GFP': 2, 'KD': 2, 'ZnT4': 1, 'reactive oxygen species': 1, 'SYT7': 1, 'Zn(2+) transporters': 1, 'retinal pigmented epithelial 1': 1, 'Fe2': 2, 'TRPML2': 1}
#    author_entity_dict = {'TRPML1 KD': 4, 'Ca(2+) release': 1, 'EBP50': 1, 'PLCgamma': 1, 'lysosomal ion channel TRPML1': 1, 'transient receptor potential': 1, 'TRPC2': 9, 'transcription factor MTF-1': 1, 'ML1': 10, 'Gb3': 3, 'TRP': 16, 'Toc': 1, 'alpha-galactosidase': 1, 'proapoptotic protein Bax': 1, 'trp2 mutant': 1, 'phospholipase C': 1, 'CatB': 3, 'lysosomal ion homeostasis': 1, 'calmodulin': 1, 'IP(3) receptors': 1, 'LRRK2': 4, 'TRPML3': 1, 'lysosomal SNARE proteins': 1, 'mitochondrial Ca2+': 1, 'GPCR': 3, 'plasma membrane receptors': 1, 'TRPML1': 28, 'zinc transporter ZnT4': 1, 'InaD': 1, 'ROS': 1, 'VAMP7 KD': 1, 'apolipoprotein B hydrolysis in MLIV': 1, 'Ca(2+) release channels': 1, 'caveolin': 1, 'GPI': 1, 'lysosomal enzymes': 1, 'TRPC': 15, 'leucine-rich repeat kinase 2': 1, 'LRRK2 G2019S': 1, 'tyrosine kinase receptors': 1, 'caspase': 1, 'NEHRF': 1, 'TRP family': 3, 'TRP2': 1, 'TRPML1in zinc transport': 1, 'MCOLN1': 1, 'G protein-coupled receptors': 1, 'R1441C': 1, 'scaffolding proteins': 1, 'lysosomal protease cathepsin B': 1, 'TRPML1in': 1, 'RPE1': 1, 'G protein coupled receptors': 1, 'synaptotagmin VII': 1, 'ROS chelator': 1, 'VAMP7': 1, 'transient receptor potential mucolipin 1': 2, 'GFP': 2, 'KD': 2, 'ZnT4': 1, 'reactive oxygen species': 1, 'SYT7': 1, 'Zn(2+) transporters': 1, 'retinal pigmented epithelial 1': 1, 'Fe2': 2, 'TRPML2': 1}
    

    total_entity_count = 0
    for entity in author_entity_dict:
        entity_count = author_entity_dict[entity]
        total_entity_count += entity_count
    
    return author_entity_dict, total_entity_count , tfidf_lists

def normalize_dict(author_entity_dict):
    pass


def get_ent_frequency(author_entity_dict, total_entity_count):
    frequency_dict = {}
    for entity in author_entity_dict:
        frequency = float(author_entity_dict[entity]) / float(total_entity_count)
        frequency = round(frequency, 6)
        frequency_dict[entity] = frequency

    sorted_frequency_list = sorted(frequency_dict.items(), key=operator.itemgetter(1), reverse= True)
    for x in sorted_frequency_list:
        print x    
    return sorted_frequency_list

def get_tfidf(tfidf_lists):
    tfidf_vectorizer = TfidfVectorizer()
        
    tfidf_lists_rm_blank = []
    tot_num_entities = 0
    for l in tfidf_lists:
        tot_num_entities += len(l)
        if len(l) > 0:
            tfidf_lists_rm_blank.append(l)
    
    print tfidf_lists_rm_blank
    tfidf_vect_input = [' '.join(x) for x in tfidf_lists_rm_blank]
    print 'tfidf_vect_input', tfidf_vect_input
    
    tfidf = tfidf_vectorizer.fit_transform(tfidf_lists_rm_blank)
    print 'tfidf', tfidf
    
    
    

def main(docs_dict, Author, current_author_num, tot_authors):
    author_entity_dict, total_entity_count, tfidf_lists = get_author_entity_dict(docs_dict, Author, current_author_num, tot_authors)
#    normalized_author_entity_dict = normalize_dict(author_entity_dict)
    entity_frequency_list = get_ent_frequency(author_entity_dict, total_entity_count)
    
#    tfidf_lists = [[], ['AC genotype', 'GCLC gene', 'rs6458939', 'glutamate-cysteine ligase', 'glutathione S-transferase alpha 3', 'GSTA3', 'CC genotype', 'glutathione-S-transferase', 'GST', 'GCLC gene'], [], ['MGMT', 'MTHFS', 'CBS', 'MGMT', 'MTHFS', 'CBS', 'DNMT3L genes'], [], ['Bmp6', 'Bone morphogenetic protein 6', 'Bmp6', 'Bmp6', 'Bmp6', 'Bmp6 allele']]
#    tfidf_list = get_tfidf(tfidf_lists)
    
    
    return entity_frequency_list, tfidf_lists



if __name__=="__main__":
    docs_dict = {}
    Author = ''
    current_author_num = ''
    tot_authors = ''
    main(docs_dict, Author, current_author_num, tot_authors)
