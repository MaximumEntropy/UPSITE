'''
Created on Mar 1, 2014

@author: Adam
'''
# from pyquery import PyQuery as pq
import urllib2
import xml.etree.ElementTree as ET
import urllib

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def get_xml(url):
    if url:
        file = urllib.urlopen(url)
        xml = file.read()
        file.close()
        return xml

def make_search_url(base_URL, q, articles, addition):
    max_papers = "&retmax=%d" % articles
    title_abstract_add = "[tiab]"
    search_url_add = "esearch.fcgi?db=pubmed&term=(%s)" % addition
    url = base_URL + search_url_add + max_papers
    return url

def get_ID_list(xml, articles):
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
        print("No Papers with your query were found on pubmed")
    sliced_id_list = ids[:articles]

    return sliced_id_list

def make_fetch_url(base_URL, get_abstract_portion_URL, ids, articles):
    if ids["papers_to_download"]:
        max_papers = "&retmax=%d" % articles
        fetch_id_string = ",".join(ids["papers_to_download"])
        fetch_url_add = "efetch.fcgi?db=pubmed&id=%s" % fetch_id_string
        full_url = base_URL + fetch_url_add + get_abstract_portion_URL + max_papers
        return full_url
    else:
        max_papers = "&retmax=%d" % articles
        fetch_id_string = ",".join(ids["papers_to_download"])
        fetch_url_add = "efetch.fcgi?db=pubmed&id=%s" % fetch_id_string
        full_url = base_URL + fetch_url_add + get_abstract_portion_URL + max_papers
        return None


def get_info_from_docs_xml(xml, ids):
       
    root = ET.fromstring(xml)
    def findall(whattofind):  # closure function -- http://en.wikipedia.org/wiki/Closure_%28computer_programming%29
        listofelements = []
        for b in root.findall(whattofind):
            
            c = b.text
            if isinstance(c, unicode):
                c = c.encode('ascii', 'ignore')  # Note: ignores unicode, does not keep unicode letters
            listofelements.append(c)
        return listofelements
    
    id_list = findall(".//ArticleId[@IdType='pubmed']")
    title_list = findall(".//ArticleTitle")
    abstract_list = findall(".//AbstractText")
    authors_list = []

    return_dict = {"fetched_id_list" : id_list, "title_list":title_list, "abstract_list":abstract_list, "authors_list": authors_list}
    return return_dict

def get_info_from_PubMed(q, num_articles, addition):  # Creates URL to search PubMed
    base_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    get_abstract_portion_URL = "&rettype=abstract"
    search_url = make_search_url(base_URL,q, num_articles, addition)
        
    if len(search_url) > 2000:
        print '---------------------URL TOO LONG-----------------------------'
        return_dict = {}
        return return_dict
    id_xml_as_String = get_xml(search_url)
    full_ID_List = get_ID_list(id_xml_as_String, num_articles)
    return full_ID_List
        
 

def main(q, num_articles, q_search_string, evaluation_mode):
    if evaluation_mode == 'yes':
        id_list = []
    elif evaluation_mode == 'no':
        id_list = get_info_from_PubMed(q, num_articles, q_search_string)
    return id_list

if __name__=="__main__":
    from optparse import OptionParser
    optparser = OptionParser(description="Get XML from PubMed")
    optparser.add_option("-q", "--query1", default='TERT', dest="q1", help="query1")
    optparser.add_option("-w", "--query2", default='MDM2', dest="q2", help="query2")
    optparser.add_option("-n", "--numpapers", default=2, dest="articles", help="Number of Pubmed Papers")
    optparser.add_option("-o", "--output", default="/home/ubuntu/Documents/upsite/protein.txt", dest="output", help="output directory")
    (options, args) = optparser.parse_args()

    main(options.q1, options.q2, options.articles)







