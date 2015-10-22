import pmids
#import tees_wrapper
#import pmids_xml
import classify
import tempfile
import os



def main(q1, q2, articles):
    id_list = pmids.main(q1,q2,articles)
#    pmid_xml = pmids_xml.main(q1,q2,articles)
#   tees_wrapper.main(pmid_xml)
    print id_list
    classify.classify('9668063','GE11','home/ubuntu/output/one_at_a_time/oneatatime')

#    for pmid in id_list: 
#	file_path = '/home/ubuntu/output/pmids/%s' % pmid
#	if os.path.isdir(file_path):
#	    continue
#	else:
#            classify.classify(pmid,'GE11',file_path)    


#    f = tempfile.NamedTemporaryFile()
#    try:
#        print 'temp:', temp
#        print 'temp.name:', temp.name
#    f.write(pmid_xml)
    print 'a' 






if __name__=="__main__":
    from optparse import OptionParser
    optparser = OptionParser(description="Get XML from PubMed")
    optparser.add_option("-q", "--query1", default='TERT', dest="q1", help="query1")
    optparser.add_option("-w", "--query2", default='MDM2', dest="q2", help="query2")
    optparser.add_option("-n", "--numpapers", default=50, dest="articles", help="Number of Pubmed Papers")
    optparser.add_option("-o", "--output", default="/home/ubuntu/Documents/upsite/protein.txt", dest="output", help="output directory")
    (options, args) = optparser.parse_args()

    main(options.q1, options.q2, options.articles)
