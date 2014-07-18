
# coding: utf-8

# In[174]:

import pandas as pd
from IPython.display import HTML
import matplotlib as pl
import pylab as pyl
import matplotlib.pyplot as pyp
import numpy as np
from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import matplotlib.pyplot as pyp
import StringIO
from rosetta import *
rosetta.init()
from fuzzywuzzy import *
from fuzzywuzzy import fuzz
from Levenshtein import distance
from best_matches import overlapping_nmer
from best_matches import best_nmers
from best_matches import find_motif_nmer
from best_matches import get_reference_sequence
from best_matches import reverse_dictionary
from best_matches import get_cysteines
# from best_matches import get_all_motif_positions
from best_matches import motif_frequency
from collections import Counter
import re


# In[175]:

# find_motif_nmer('H5 complete.txt', 'N[^P][ST]', 3, 5)


# In[176]:

ref_dict = overlapping_nmer(get_reference_sequence('H5 complete.txt'), 13)
ref_dict[21]


# In[177]:

def hamd(string1, string2):
    score = 0
    for n, char in enumerate(string1):
        if string2[n] != char:
            score += 1
    return score
hamd('ATGATC', 'ATGAAC')


# In[178]:

# # def reverse_dictionary(dictionary):
# #     new_dictionary = {}
# #     for key, value in dictionary.items():
#         new_dictionary[value] = key
#     return new_dictionary


# In[179]:

# reverse_dictionary(ref_dict)


# In[180]:

pos_freqs = motif_frequency(best_nmer_matches)


# In[180]:




# In[173]:

best_nmer_matches = best_nmers('H5 complete.txt','N[^P][ST]', 3, 5, 0, 1, 20)


# In[181]:

pos_freqs.sort()
pos_freqs


# In[182]:

long_ref_seq = get_reference_sequence('H5 complete.txt')
long_ref_dict = overlapping_nmer(long_ref_seq, 13)
long_reverse_ref_dict = reverse_dictionary(long_ref_dict)
half_buffer = 20


# In[183]:

aichi_ref = get_reference_sequence('A_Aichi_2_68.txt')
aichi_dict = overlapping_nmer(aichi_ref, 13)


# In[195]:

long_cys = get_cysteines('H5 complete.txt', 5, 2, 20)
aichi_ref = open('A_Aichi_2_68.txt')
for accession in long_cys.keys():
    for pos in long_cys[accession].keys():
        if long_cys[accession][pos] == []:
            long_cys[accession].pop(pos)
aichi = aichi_ref.readlines()[1:-1]
aichi_seq = ''
for line in aichi:
    aichi_seq += line[:-1]
aichi_num = list(enumerate(aichi_seq))
aichi_cys = []
for char in aichi_num:
    if char[1] == 'C':
        aichi_cys.append(char)
# aichi_cys
# aichi_seq.index('C')
long_cys


# In[207]:

positions = []
for accession in long_cys.values():
    for scored_match in accession.values():
        positions.append(scored_match[0][1])
cysteine_counts = Counter(positions)
# cysteine_counts
cysteine_freqs = []
for item in cysteine_counts:
    if float(cysteine_counts[item])/len(long_cys.keys()) >= 0.9:
        cysteine_freqs.append([item, float(cysteine_counts[item])/len(long_cys.keys())])
cysteine_freqs


# In[208]:

import matplotlib.pyplot as plt
plt.figure(0, figsize=(12,4))
plt.bar([row[0] for row in pos_freqs], [row[1] for row in pos_freqs], width=5, color='blue')
plt.ylabel('Possible NLG Frequency')
plt.xlabel('Position on Reference H5', labelpad = 50)
plt.title('H5 Cysteine Frequency (Longest Reference)')
for cys in cysteine_freqs:
    plt.annotate('C'+str(cys[0]+1), xy=(cys[0]+1, 0), xytext=(cys[0]+1,-0.2),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.25, width=1, headwidth=5))
plt.show()


# In[209]:

plt.figure(1, figsize=(12, 4))
plt.axes()
plt.title('H5 NLG Position Frequency (Aichi Reference)')
plt.xlabel('Position on Reference', labelpad=50)
plt.ylabel('Possible NLG Frequency')
plt.bar([row[0] for row in pos_freqs], [row[1] for row in pos_freqs], width=5, color='blue')
for even in range(0, len(aichi_cys), 2):
    plt.annotate('C'+str(aichi_cys[even][0]+1), xy=(aichi_cys[even][0]+1, 0.03), xytext=(aichi_cys[even][0]+1,-0.2),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.25, width=1, headwidth=5))
for odd in range(1, len(aichi_cys), 2):
    plt.annotate('C'+str(aichi_cys[odd][0]+1), xy=(aichi_cys[odd][0]+1, 0), xytext=(aichi_cys[odd][0]+1,-0.1),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.5, width=1, headwidth=5))


# In[ ]:



