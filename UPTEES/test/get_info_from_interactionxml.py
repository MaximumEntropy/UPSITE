'''
Created on Feb 2, 2015

@author: adam
'''
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import gzip

def get_info_from_interaction_xml(file_paths):
    final_dict = {}
    for file_path in file_paths:
        try:
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
                if 'conf' in elem.attrib:
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
        except (IOError, ET.ParseError):
            continue
           
    return final_dict    

def main():
    file_paths = ['/home/adam/workspace/TEESoutput/testing/REL11/genes/mdm2/25629974-pred.xml.gz', '/home/adam/workspace/TEESoutput/testing/REL11/genes/mdm2/25629302-pred.xml.gz']
    final_dict = get_info_from_interaction_xml(file_paths)
    print 'final dict', final_dict
    
    
if __name__ == "__main__":
    main()