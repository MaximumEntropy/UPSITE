# from pyquery import PyQuery as pq
import urllib2
import os
import math
import get_genes
import time
import csv
import sys
import pickle
from random import randrange
from xml.dom.minidom import parseString
from cookielib import domain_match
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def get_xml_minidom(url):
    xml = urllib2.urlopen(url)
    dom = parseString(xml.read())

    print dom
    
    return dom

def get_xml(url):
    if url:
        print 'url', url
        file = urllib2.urlopen(url)
        xml = file.read()
        file.close()
        return xml

def make_search_url(base_URL, articles, university, type_field):
    max_papers = "&retmax=%d" % articles
#    title_abstract_add = "[tiab]"
    search_url_add = "esearch.fcgi?db=pubmed&term="      #University%20of%20Pittsburgh%5BAffiliation%5D"
    url = base_URL + search_url_add + type_field + university + max_papers
    return url

def get_ID_list(xml):
    try:
        root = ET.fromstring(xml)
        ID_List_ofElements = root.findall("./IdList/Id")
        ids = []
        for element in ID_List_ofElements:
            singleID_string = ET.tostring(element, method='text')
            singleID_string_stripped = singleID_string.replace("\n", "")
            ids.append(singleID_string_stripped)
    except AttributeError:
        ids = []
        print("No Papers with both queries were found on PubMed")
    existing_papers = []  # Use this in the future to make database of existing IDs 
    papers_to_download = []
    for ind_id in ids:
        papers_to_download.append(ind_id)

    full_ID_List = {"existing_papers":existing_papers,
                                    "papers_to_download":papers_to_download}
    return full_ID_List

def make_fetch_url(base_URL, get_abstract_portion_URL, ids, articles):

    max_papers = "&retmax=%d" % articles
    fetch_id_string = ",".join(ids)
    fetch_url_add = "efetch.fcgi?db=pubmed&id=%s" % fetch_id_string
    full_url = base_URL + fetch_url_add + get_abstract_portion_URL + max_papers
    return full_url
#     else:
#         max_papers = "&retmax=%d" % articles
#         fetch_id_string = ",".join(ids["papers_to_download"])
#         fetch_url_add = "efetch.fcgi?db=pubmed&id=%s" % fetch_id_string
#         full_url = base_URL + fetch_url_add + get_abstract_portion_URL + max_papers
#         return None


def get_info_from_docs_xml(xml_list, ids):
    
    return_dict = {}
    for xml in xml_list:
#        root = ET.fromstring(xml)
        tree = ET.ElementTree(ET.fromstring(xml))
        root = tree.getroot()
        
        PMIDS = []
        for elem in tree.iter(tag = 'MedlineCitation'):
            for child in elem:
                if child.tag == 'PMID':
                    PMID = child.text
                    PMIDS.append(PMID)
        
        partial_return_dict_values = []
        for elem in tree.iter(tag='Article'):
            individual_paper_info = []
            for child in elem:
                if child.tag == 'Abstract':
                    for grandchild in child:
                        abstract_text = grandchild.text
           #             abstract_text = 'PLACE HOLDER'
                        individual_paper_info.append(abstract_text)
                if child.tag == 'ArticleTitle':
                    ArticleTitle = child.text
                if child.tag =='AuthorList':
                    Authors = []
                    for grandchild in child:
                        if grandchild.tag == 'Author':
                            Last_Name = ''
                            First_Name = ''
                            Affiliation = ''
                            for sub_branch in grandchild:
                                if sub_branch.tag == 'LastName':
                                    Last_Name = sub_branch.text
                                if sub_branch.tag == 'ForeName':
                                    First_Name = sub_branch.text
                                if sub_branch.tag == 'Affiliation':
                                    Affiliation = sub_branch.text
                            single_author = (First_Name, Last_Name, Affiliation)
                            Authors.append(single_author)
                    individual_paper_info.append(Authors)
            partial_return_dict_values.append(individual_paper_info)
        
        if len(PMIDS) != len(partial_return_dict_values):
            print 'Different number of PMIDS and results!!!!'
            
        for num, ID in enumerate(PMIDS):
            if ID in return_dict:
                pass
            else:
                return_dict[ID] = partial_return_dict_values[num]
 #           partial_return_dict[ID] = partial_return_dict_values[num]
        
    return return_dict
        
#     return_dict = {}
#     for xml in xml_list:
#         root = ET.fromstring(xml)    
#         def findall(whattofind):  # closure function -- http://en.wikipedia.org/wiki/Closure_%28computer_programming%29
#             listofelements = []
#             for b in root.findall(whattofind):  
#                 c = b.text
#                 if isinstance(c, unicode):
#                     c = c.encode('ascii', 'ignore')  # Note: ignores unicode, does not keep unicode letters
#                 listofelements.append(c)
#             return listofelements
#     
#         id_list = findall(".//ArticleId[@IdType='pubmed']")
# #    if id_list > 0:
# #        print str(len(id_list)) + " papers with co-occurrence found"
#         title_list = findall(".//ArticleTitle")
#         abstract_list = findall(".//AbstractText")
#         last_name_list = findall(".//LastName")
#         first_name_list = findall(".//ForeName")
#         affiliation_list = findall(".//Affiliation")
#         
#         if not return_dict:
#             return_dict = {"fetched_id_list" : id_list, "title_list":title_list, "abstract_list":abstract_list, "last_name_list":last_name_list,
#                        "first_name_list":first_name_list, "affiliation_list":affiliation_list}
#         elif return_dict:
#             return_dict['fetched_id_list'] += id_list
#             return_dict['title_list'] += title_list
#             return_dict['abstract_list'] += abstract_list
#             return_dict['last_name_list'] += last_name_list
#             return_dict['first_name_list'] += first_name_list
#             return_dict['affiliation_list'] += affiliation_list
#     return return_dict

def get_info_from_PubMed(articles, university, type_field):  # Creates URL to search PubMed
    base_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    get_abstract_portion_URL = "&rettype=abstract"
    
    search_url = make_search_url(base_URL, articles, university, type_field)
    
    if len(search_url) > 2000:
        return_dict = {}
        raise Exception('Length of search URL is greater than 2000')
        time.sleep(3)
        return return_dict
    
    id_xml_as_String = get_xml(search_url)
    full_ID_List = get_ID_list(id_xml_as_String)
    id_list = full_ID_List['papers_to_download']
    print '{0} Publications found! ' .format(len(id_list))
    info_from_PubMed = {}
    if full_ID_List["papers_to_download"]:
        size = 200
        id_lists = [id_list[i:i+size] for i in range(0, len(id_list), size)]    # splits total id lists into sets of 200 or less so url isn't too long 
        docs_xml_list = []
        for sub_list in id_lists:
            fetch_url = make_fetch_url(base_URL, get_abstract_portion_URL, sub_list, articles)
            docs_xml = get_xml(fetch_url)
            docs_xml_list.append(docs_xml)
#            time.sleep(randrange(3,6))               # commented out bc banner takes so long, probably dont need to sleep
        info_from_PubMed = get_info_from_docs_xml(docs_xml_list, full_ID_List)
    else:
        print 'Search returned no results'
        info_from_PubMed = {}
        
    return info_from_PubMed
        
        
def make_paper_objects(dict_of_info):
    """takes in dict of info, returns dictionary paper objects like
        ["paper_id#"]: paper object """
    ID_paper_obj_dict = {}
    if "fetched_id_list" in dict_of_info:
        fetched_id_list = dict_of_info["fetched_id_list"]
        title_list = dict_of_info["title_list"]
        abstract_list = dict_of_info["abstract_list"]
        author_list = dict_of_info["authors_list"]
    print author_list

def get_university(university):
    university = university.replace(' ', '%20')
    affiliation = '%5BAffiliation%5D'
    univ_search = university+affiliation
    return univ_search

def get_name_list(file_location):       # use if first, last separated by a space as in author_list
        try:         # if names split by spaces
            with open(file_location, 'r') as my_file:
                reader=csv.reader(my_file, delimiter = '\t')
                Author_strings = []
                Authors = []
                for row in reader:
                    for full_name in row:
                        name = full_name.split(' ')
                        Author_first =  name[0]
                        Author_last = name[1]
                        Author = [Author_first,Author_last]
                        Authors.append(Author)
                        Author_string = Author_last + '%2C%20' +Author_first + '%5BAuthor%5D'
                        Author_strings.append(Author_string)
        except IndexError:      # if names split by tabs
            with open(file_location, 'r') as my_file:
                reader=csv.reader(my_file, delimiter = '\t')
                Author_strings = []
                Authors = []
                row_num = 0
                for row in reader:
                    if row_num > 1534:   # row before first empty row
                        Author_last =  row[0]
                        Author_first = row[1]
                        Authors.append(row)
                        Author_string = Author_last + '%2C%20' + Author_first + '%5BAuthor%5D'
                        Author_strings.append(Author_string)
                        row_num += 1
                    else:
                        row_num += 1
                    
        return Author_strings, Authors
                

def main():
################search_type############################
#   select (1) for university+author, OR 
#   select (2) for university+gene   ####PROBABLY NOT GOING TO USE
#   select (3) if using a list of authors in one university

    search_type = '3'


    if int(search_type) == 1:
        articles = 10
        university = 'University of Pittsburgh'
        Author_first = 'Catalina'
        Author_last = 'Cleves Bayon'
        Author_string = Author_last + '%2C%20' +Author_first + '%5BAuthor%5D'
        Author = Author_first + ' ' + Author_last
        current_author_num = 1
        tot_authors = 1
        university = get_university(university)
        docs_dict  = get_info_from_PubMed(articles, university, Author_string)
        print 'docs dict', docs_dict
        entity_frequency_list, tfidf_lists = get_genes.main(docs_dict, Author, current_author_num, tot_authors)
        
    if int(search_type) == 2:
        gene = ''
        gene_field = '{0}%20AND%20' .format(gene)
        docs_xml  = get_info_from_PubMed(articles, university, gene_field , Author_first, Author_last)
        
    if int(search_type) == 3:
        articles = 40
#        file_location = '/home/adam/workspace/TEES/text_files/Author_Lists/Author_list'  #Madhavi suggested example
        file_location = '/home/adam/workspace/TEES/text_files/Author_Lists/SOM_Faculty.tsv'  #School of Medicine Faculty
#        file_location = '/home/adam/workspace/TEES/text_files/Author_Lists/SOM_Faculty_short.tsv'  #School of Medicine Faculty, cut short for testing
        university = 'University of Pittsburgh'
        university = get_university(university)
        Author_strings, Authors = get_name_list(file_location)
        Author_keys = []
        for Author in Authors:
            full_name = Author[0] + ' ' + Author[1]
            Author_keys.append(full_name)
        print 'author keys', Author_keys
        frequency_values = []
        for num, author_string in enumerate(Author_strings):
            indiv_dict = {}
            tot_authors = len(Author_strings)
            current_author_num = int(num) + 1
            Author_first_last = Authors[num]
            Author = Author_first_last[0] + ' ' + Author_first_last[1]
            docs_dict  = get_info_from_PubMed(articles, university, author_string)   
            entity_frequency_list, tfidf_lists = get_genes.main(docs_dict, Author, current_author_num, tot_authors)
            frequency_values.append(entity_frequency_list)
            indiv_dict[Author_keys[num]] = entity_frequency_list
            base_name = os.path.basename(file_location)
            output_file = '/home/adam/workspace/TEES/text_files/Author_Lists/output/{0}' .format(base_name)
            with open(output_file, 'a') as f:  # Just use 'w' mode in 3.x
                w = csv.writer(f)
                w.writerows(indiv_dict.items())
            tfidf_output_file = output_file + '_tfidf'
            with open(tfidf_output_file, 'a') as f2:
                pickle.dump(tfidf_lists, f2)
                
            
        Author_frequency_dict = dict(zip(Author_keys,frequency_values))
  #      Author_frequency_dict = {'Kirill Kiselyov': [('TRP', 0.148148), ('TRPC', 0.138889), ('TRPML1', 0.12963), ('ML1', 0.092593), ('TRPC2', 0.083333), ('GPCR', 0.027778), ('CatB', 0.027778), ('TRP family', 0.027778), ('TRPML1 KD', 0.009259), ('Ca(2+) release', 0.009259), ('EBP50', 0.009259), ('PLCgamma', 0.009259), ('plasma membrane receptors', 0.009259), ('trp2 mutant', 0.009259), ('proapoptotic protein Bax', 0.009259), ('phospholipase C', 0.009259), ('lysosomal SNARE proteins', 0.009259), ('lysosomal ion homeostasis', 0.009259), ('calmodulin', 0.009259), ('IP(3) receptors', 0.009259), ('mitochondrial Ca2+', 0.009259), ('G protein-coupled receptors', 0.009259), ('InaD', 0.009259), ('synaptotagmin VII', 0.009259), ('VAMP7 KD', 0.009259), ('Ca(2+) release channels', 0.009259), ('caveolin', 0.009259), ('lysosomal enzymes', 0.009259), ('tyrosine kinase receptors', 0.009259), ('caspase', 0.009259), ('TRP2', 0.009259), ('MCOLN1', 0.009259), ('scaffolding proteins', 0.009259), ('lysosomal protease cathepsin B', 0.009259), ('TRPML3', 0.009259), ('transient receptor potential', 0.009259), ('G protein coupled receptors', 0.009259), ('TRPML2', 0.009259), ('VAMP7', 0.009259), ('SYT7', 0.009259), ('Zn(2+) transporters', 0.009259), ('apolipoprotein B hydrolysis in MLIV', 0.009259), ('NEHRF', 0.009259)], 'Madhavi Ganapathiraju': [('ANKLE1', 0.25), ('ORAOV1', 0.25), ('TMEM45B', 0.25), ('human protein', 0.25)]}
        
        print Author_frequency_dict
        
#         base_name = os.path.basename(file_location)
#         output_file = '/home/adam/workspace/TEES/text_files/Author_Lists/output/{0}' .format(base_name)
#         with open(output_file, 'wb') as f:  # Just use 'w' mode in 3.x
#             w = csv.writer(f)
#             w.writerows(Author_frequency_dict.items())
            
    
if __name__=="__main__":
    main()




