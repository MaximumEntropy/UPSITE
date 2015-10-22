'''
Created on Jan 20, 2015

@author: adam
'''
import csv
import os

current_dir = os.getcwd()
file_entry = '/home/adam/workspace/TEES/text_files/final_protein_word_dict_known60.txt'

protein_word_dict = {}
with open(file_entry, 'r') as my_file:
    reader = csv.reader(my_file, delimiter='\t')
    proteins = []
    lists = []
    for row in reader:
        proteins.append(row[:1][0])
        lists.append(row[1:])

a = zip(*[iter(proteins)] * 2)
b = zip(*[iter(lists)] * 2)
print a
