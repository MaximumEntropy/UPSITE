import pmids
import classify
import os
import gzip
import queries
import Cosine_Sim
import csv
import time
import math
from os import listdir
from os.path import isfile, join
from random import randrange
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def run_tees_batch(q, id_list):
    current_dir = os.getcwd()
    gene_dir = current_dir + 'output/batch/genes/{0}' .format(q)
    already_downloaded_pmids = []
    if os.path.exists(gene_dir):
        onlyfiles = [ f for f in listdir(gene_dir) if isfile(join(gene_dir,f)) ]
        file_name_list = list(file_name for file_name in onlyfiles)
        for x in file_name_list:
            split_at_dash = x.split('-')
            for piece_of_file_name in split_at_dash:
                if piece_of_file_name[0].isdigit():
                    if piece_of_file_name in already_downloaded_pmids:
                        continue
                    else:
                        already_downloaded_pmids.append(str(piece_of_file_name)) 
 
    ignore_list_path = current_dir + '/text_files/batch/id_ignore_list.txt'
    pmid_ignore_list = []
    with open (ignore_list_path, 'r') as f:
        reader = csv.reader(f,delimiter='\t')
        if reader:
            for pmid_list in reader:
                pmid_ignore_list = pmid_list
        else:
            pmid_ignore_list = []
            
    pmid_run_list = [] 
    for pmid in id_list:
        if pmid in pmid_ignore_list:
            continue
        elif pmid in already_downloaded_pmids:
            continue
        else:
            pmid_run_list.append(pmid)
            
    #===========================================================================
    # float_num = float(len(pmid_run_list))/10
    # rounded_up_num = math.ceil(float_num)
    # list_of_pmid_lists = []
    # for x in range(int(rounded_up_num)):
    #     list_of_pmid_lists.append(pmid_list[10*x:(10*x)+9])       
    #===========================================================================
###################       list_of_pmid_lists.append(pmid_list[10*x:(10*x)+10])  

    float_num = float(len(pmid_run_list))/25
    rounded_up_num = math.ceil(float_num)
    list_of_pmid_lists = []
    for x in range(int(rounded_up_num)):
        list_of_pmid_lists.append(pmid_list[25*x:(25*x)+24])    
    
    file_path_list = [] 
    if not list_of_pmid_lists:
        return file_path_list
    else:
        run_single_pmids = [] 
        for plist in list_of_pmid_lists:
            file_path = current_dir + 'output/batch/genes/{0}/{1}' .format(q , '-'.join(plist))
            addition = '-pred.xml.gz'
            file_path_check = file_path + addition
            file_path_list.append(file_path)
            
            if os.path.exists(file_path_check):
                print '--------------------------------SKIPPING ALREADY DOWNLOADED ABSTRACTS {0}-------------------------------------------' .format(plist)
            else:               
                try:    
                    classify.classify('-'.join(plist),'GE11',file_path)
                except (ValueError, UnicodeEncodeError, AssertionError, IndexError) as e:
                    print 'file path list remove!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print 'error,', e
                    file_path_list.remove(file_path)
                    run_single_pmids.append(plist)
        if run_single_pmids:
            single_pmids_file_path_list = run_tees(q, plist)
            file_path_list += single_pmids_file_path_list
                     
        return file_path_list
    
def run_tees(q, id_list):
    current_dir = os.getcwd()
    ignore_list_path = current_dir + '/text_files/id_ignore_list.txt'
    pmid_ignore_list = []
    with open (ignore_list_path, 'r') as f:
        reader = csv.reader(f,delimiter='\t')
        for pmid_list in reader:
            pmid_ignore_list = pmid_list

    file_path_list = []    
    for pmid in id_list:
        if pmid in pmid_ignore_list:
            pass
        else:
            file_path = current_dir + 'output/genes/{0}/{1}' .format(q ,pmid)
            addition = '-pred.xml.gz'
            file_path_check = file_path + addition
            file_path_list.append(file_path)
            if os.path.exists(file_path_check):
                print '--------------------------------SKIPPING ALREADY DOWNLOADED ABSTRACT {0}-------------------------------------------' .format(pmid)
                pass
            else:
                try:
                    classify.classify(pmid,'GE11',file_path)
                except (ValueError, UnicodeEncodeError, AssertionError, IndexError) as e:
                    file_path_list.remove(file_path)
                    with open(ignore_list_path, 'a') as f:
                        f.write(pmid + '\t')
                        f.close()
                    pass
    return file_path_list

def get_info_from_interaction_xml(file_paths):
    new_file_paths = []
    addition = '-pred.xml.gz'
    for f in file_paths:
        f+= addition
        new_file_paths.append(f)
    final_dict = {}
    for file_path in new_file_paths:
        infile = gzip.open(file_path, 'r')
        tree = ET.ElementTree(file=infile)

	entity_dict = {}
	trigger_dict = {}
        entity_trigger_dict = {}
        for elem in tree.iter(tag='entity'):
	    entity_trigger_dict[elem.attrib['id']]=elem.attrib['text']
            if 'source' in elem.attrib:
	        entity_dict[elem.attrib['id']] = elem.attrib['text']
	    if 'umConf' in elem.attrib:
		trigger_dict[elem.attrib['id']] = elem.attrib['text']
#        entity_trigger_dict = dict(entity_dict.items() + trigger_dict.items())
	for elem in tree.iter(tag='interaction'):
	    e1 = elem.attrib['e1']
            e2 = elem.attrib['e2']
            if (e1 in entity_trigger_dict) and (e2 in entity_trigger_dict):
                e1_text = entity_trigger_dict[e1]
	        e2_text = entity_trigger_dict[e2]
	    if (e1 in entity_dict) and (e2 in trigger_dict):
                if e1_text not in final_dict:
	            final_dict[e1_text]=[e2_text]
                else:
	            final_dict[e1_text].append(e2_text)  
	    elif (e2 in entity_dict) and (e1 in trigger_dict):
                if e2_text not in final_dict:
                    final_dict[e2_text]=[e1_text]
                else:
                    final_dict[e2_text].append(e1_text)
            elif (e1 in entity_dict) and (e2 in entity_dict):
                if e1_text not in final_dict:
                    final_dict[e1_text]=[e2_text]
                else:
                    final_dict[e1_text].append(e2_text)
                if e2_text not in final_dict:
                    final_dict[e2_text]=[e1_text]
                else:
                    final_dict[e2_text].append(e1_text)
         
    return final_dict	

def get_protein_dict(q1, q2, q1_dict, q2_dict):
    if q1 in q1_dict:
        print q1, 'in q1 dict: ', q1, q1_dict[q1]
        q1_words = q1_dict[q1]
    else: 
        print q1, 'not in q1 dict:' , q1 , q1_dict
        q1_words = []
        
    if q2 in q2_dict:   
        print q2, 'in q2 dict: ', q2, q2_dict[q2]
        q2_words = q2_dict[q2]    
    else: 
        print q2, 'not in q2 dict:' , q2_dict
        q2_words = []
        
    if q1_words and q2_words:
        protein_dict = {q1:q1_dict[q1], q2:q2_dict[q2]}
        print 'protein dict', protein_dict
        return protein_dict
    else:
        return {}
    
def get_all_words_dict(q1, q2, q1_dict, q2_dict):
    if q1_dict:
        all_q1_words = []
        for q in q1_dict:
            all_q1_words.extend(q1_dict[q])
    if q2_dict:
        all_q2_words = []
        for q in q2_dict:
            all_q2_words.extend(q2_dict[q])
    if q1_dict and q2_dict:
        all_words_dict = {q1:all_q1_words, q2:all_q2_words}
        return all_words_dict
    else:
        return {}
    
def normalize_dict(dict_in, query):
    dict_x = {}
    for k, v in dict_in.iteritems():
        v_lower = []
        for word in v:
            v_lower.append(word.lower())
        dict_x[k.lower()] = v_lower   

    return_dict = {}
    for entity in dict_x:
        try:
            if entity in query.q1_syns or entity == query.q1.lower():
                if query.q1.lower() in return_dict:
                    return_dict[query.q1.lower()] += dict_x[entity]
                else:
                    return_dict[query.q1.lower()] = dict_x[entity]
            if entity in query.q2_syns or entity == query.q2.lower():
                if query.q2.lower() in return_dict:
                    return_dict[query.q2.lower()] += dict_x[entity]
                else:
                    return_dict[query.q2.lower()] = dict_x[entity]
            else:
                pass
        except TypeError:
            pass
            
    return return_dict

def combine_dictionaries(query_dicts):
    combined_dict = {}
    for d in query_dicts:
        for k in d:
            if k not in combined_dict:
                combined_dict[k] = d[k]
            else:
                combined_dict[k] += d[k]
    return combined_dict

def print_pair_score_dict(angle_list, protein_dict, q1, q2, input_type, outputFileName):
    current_dir = os.getcwd()
    
    if input_type == 'known':
        write_path = current_dir + '/text_files/output_known_interactions.txt' 
                  
    elif input_type == 'unknown' or input_type =='random':
        write_path = current_dir  + '/text_files/output_random_interactions.txt' 
                   
    elif outputFileName:
        write_path = current_dir + '/text_files/' +str(outputFileName)
        
    else:       
        write_path = current_dir + '/text_files/' +str(input_type)  
                 
    with open (write_path,'a') as f:
        f.write(q1+'\t'+q2+'\t')
        f.write(str(angle_list))
        f.write(str('\t'))
        f.write(str(protein_dict))
        f.write(str('\n'))  
        
def get_already_batch_downloaded_list(q):        
    current_dir = os.getcwd()
    gene_dir = current_dir + 'output/batch/genes/{0}' .format(q)
    already_downloaded_pmids = []
    if os.path.exists(gene_dir):
        onlyfiles = [ f for f in listdir(gene_dir) if isfile(join(gene_dir,f)) ]
        file_name_list = list(file_name for file_name in onlyfiles)
        for x in file_name_list:
            split_at_dash = x.split('-')
            for piece_of_file_name in split_at_dash:
                if piece_of_file_name[0].isdigit():
                    if piece_of_file_name in already_downloaded_pmids:
                        continue
                    else:
                        already_downloaded_pmids.append(str(piece_of_file_name))  
                        
def get_already_downloaded_list(q):
    current_dir = os.getcwd()
    gene_dir = current_dir + 'output/genes/{0}' .format(q)
    already_downloaded_pmids = []
    if os.path.exists(gene_dir):
        onlyfiles = [ f for f in listdir(gene_dir) if isfile(join(gene_dir,f)) ]
        file_name_list = list(file_name for file_name in onlyfiles)
        for x in file_name_list:
            split_at_dash = x.split('-')
            for piece_of_file_name in split_at_dash:
                if piece_of_file_name[0].isdigit():
                    if piece_of_file_name in already_downloaded_pmids:
                        continue
                    else:
                        already_downloaded_pmids.append(str(piece_of_file_name)) 
    print 'already_downloaded_pmids', already_downloaded_pmids

    file_path_list = []    
    for pmid in already_downloaded_pmids:
        file_path = current_dir + 'output/genes/{0}/{1}' .format(q ,pmid)
        file_path_list.append(file_path)
    
    return file_path_list
        
        
def main(q1, q2, articles, batch, input_type, outputFileName, dictType):
    num_articles = int(articles)
    query = queries.main(q1,q2)
    
    if batch == "yes":
        q1_id_list, q2_id_list = pmids.main(query, num_articles)
        q1_file_paths= run_tees_batch(q1, q1_id_list)
        q2_file_paths= run_tees_batch(q2, q2_id_list)
    
    if batch == "no":
        q1_id_list, q2_id_list = pmids.main(query, num_articles)  
        q1_file_paths= run_tees(q1, q1_id_list)
        q2_file_paths= run_tees(q2, q2_id_list)

    q1_dict = get_info_from_interaction_xml(q1_file_paths)
    q2_dict = get_info_from_interaction_xml(q2_file_paths)
    
    if dictType == 'all' or dictType == 'both':
        all_words_dict = get_all_words_dict(q1, q2, q1_dict, q2_dict)
        normalized_all_words_dict = normalize_dict(all_words_dict, query)
        angle_list_all = Cosine_Sim.main(normalized_all_words_dict, q1, q2)
        print_pair_score_dict(angle_list_all, normalized_all_words_dict, q1, q2, input_type, outputFileName)        
    
    if dictType == 'protein' or dictType == 'both':
        query_dicts = [q1_dict, q2_dict]
        combined_dict = combine_dictionaries(query_dicts)
        normalized_protein_dict = normalize_dict(combined_dict, query)
        angle_list_protein = Cosine_Sim.main(normalized_protein_dict, q1, q2)
        print_pair_score_dict(angle_list_protein, normalized_protein_dict, q1, q2, input_type, outputFileName)

if __name__=="__main__":
    from optparse import OptionParser
    optparser = OptionParser(description="Get XML from PubMed")
    optparser.add_option("-q", "--q1", default='HMGB1', dest="q1", help="query1")
    optparser.add_option("-w", "--q2", default='KL', dest="q2", help="query2")
    optparser.add_option("-n", "--n", default=40, dest="articles", help="Number of Pubmed Papers to download per gene/protein")
    optparser.add_option("-o", "--outputFileName", default="", dest="outputFileName", help="output file name")
    optparser.add_option("-i", "--input", default='known', dest="input", help="single or 50examples or known or unknown")
    optparser.add_option("-b", "--batch", default="no", dest="batch", help="y for yes or n for no")
    optparser.add_option("-d", "--dict", default="all", dest="dictType", help="all or protein or both")
    (options, args) = optparser.parse_args()
    
    if options.input == 'single':
        main(options.q1, options.q2, options.articles, options.batch, options.input, options.outputFileName, options.dictType)
        
    elif options.input =='50examples':
        file_entry = r'/text_files/madhavi_example_protein_interactions.txt'
        current_dir = os.getcwd()
        dir_entry = current_dir + file_entry
        list_of_protein_pairs = []
        with open(dir_entry, 'r') as my_file:
            reader = csv.reader(my_file, delimiter='\t')
            for row in reader:
                list_of_protein_pairs.append(row)
        for protein_pair in list_of_protein_pairs:
            q1 = protein_pair[0]
            q2 = protein_pair[1]   
#            t = randrange(1,3)
#            time.sleep(t) 
            main(q1, q2, options.articles, options.batch, options.input, options.outputFileName, options.dictType)
            
    elif options.input =='known':
        file_entry = r'/text_files/known_interactions.tsv'
        current_dir = os.getcwd()
        dir_entry = current_dir + file_entry
        list_of_protein_pairs = []
        with open(dir_entry, 'r') as my_file:
            reader = csv.reader(my_file, delimiter='\t')
            for row in reader:
                list_of_protein_pairs.append(row)
        for protein_pair in list_of_protein_pairs:
            q1 = protein_pair[0]
            q2 = protein_pair[1]   
#            t = randrange(1,2)
#            time.sleep(t) 
            main(q1, q2, options.articles, options.batch, options.input, options.outputFileName, options.dictType)

    elif options.input =='unknown' or options.input =='random':
        file_entry = r'/text_files/random_interactions.tsv'
        current_dir = os.getcwd()
        dir_entry = current_dir + file_entry
        list_of_protein_pairs = []
        with open(dir_entry, 'r') as my_file:
            reader = csv.reader(my_file, delimiter='\t')
            for row in reader:
                list_of_protein_pairs.append(row)
        for protein_pair in list_of_protein_pairs:
            q1 = protein_pair[0]
            q2 = protein_pair[1]   
#            t = randrange(1,2)
#            time.sleep(t) 
            main(q1, q2, options.articles, options.batch, options.input, options.outputFileName, options.dictType)
        
        