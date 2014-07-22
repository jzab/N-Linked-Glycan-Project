
# coding: utf-8

# In[1]:

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
from best_matches import get_other_cysteines
from best_matches import best_other_nmers
# from best_matches import get_all_motif_positions
from best_matches import motif_frequency
from collections import Counter
import re
import csv


# In[2]:

ref_dict = overlapping_nmer(get_reference_sequence('H5 complete.txt'), 13)
best_nmer_matches = best_nmers('H5 complete.txt','N[^P][ST]', 3, 5, 0, 1, 20)
pos_freqs = motif_frequency(best_nmer_matches)
pos_freqs.sort()
# pos_freqs


# In[3]:

long_ref_seq = get_reference_sequence('H5 complete.txt')
long_ref_dict = overlapping_nmer(long_ref_seq, 13)
long_reverse_ref_dict = reverse_dictionary(long_ref_dict)
half_buffer = 20


# In[4]:

long_cys = get_cysteines('H5 complete.txt', 5, 2, 20)
for accession in long_cys.keys():
    for pos in long_cys[accession].keys():
        if long_cys[accession][pos] == []:
            long_cys[accession].pop(pos)
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
# cysteine_freqs


# In[5]:

import matplotlib.pyplot as plt
plt.figure(0, figsize=(12,4))
plt.bar([row[0] for row in pos_freqs], [row[1] for row in pos_freqs], width=5, color='blue')
plt.ylabel('Possible NLG Frequency')
plt.xlabel('Position on Reference H5', labelpad = 50)
plt.title('H5 Possible NLG Frequency (Longest Reference)')
for cys in cysteine_freqs:
    plt.annotate('C'+str(cys[0]+1), xy=(cys[0]+1, 0), xytext=(cys[0]+1,-0.2),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.25, width=1, headwidth=5))
plt.annotate('Figure 1: Mapping of Possible NLG Frequency on Longest Reference Sequence with 95% Conserved Cysteines as Reference Position Markers', xy=(0,0), xytext=(0,-0.4), size='small')
plt.show()


# In[6]:

aichi_ref = get_reference_sequence('A_Aichi_2_68.txt')
aichi_dict = overlapping_nmer(aichi_ref, 13)


# In[7]:

aichi_ref = open('A_Aichi_2_68.txt')
aichi = aichi_ref.readlines()[1:-1]
aichi_seq = ''
for line in aichi:
    aichi_seq += line[:-1]


# In[8]:

best_aichi_nmer_matches = best_other_nmers('H5 complete.txt', aichi_seq, 'N[^P][ST]', 3, 5, 0, 1, 20)
aichi_pos_freqs = motif_frequency(best_aichi_nmer_matches)
aichi_pos_freqs.sort()


# In[9]:

# aichi_num = list(enumerate(aichi_seq))
# aichi_cs = []
# for char in aichi_num:
#     if char[1] == 'C':
#         aichi_cs.append(char)
# # aichi_cs
# # aichi_seq.index('C')
# # long_cys


# In[10]:

aichi_cys = get_other_cysteines('H5 complete.txt', aichi_seq, 5, 4, 20)
for accession in aichi_cys.keys():
    for pos in aichi_cys[accession].keys():
        if aichi_cys[accession][pos] == []:
            aichi_cys[accession].pop(pos)
# aichi_cys
aichi_positions = []
for accession in aichi_cys.values():
    for scored_match in accession.values():
        aichi_positions.append(scored_match[0][1])
# aichi_positions
aichi_cysteine_counts = Counter(aichi_positions)
aichi_cysteine_counts
aichi_cysteine_freqs = []
for item in aichi_cysteine_counts:
    if float(aichi_cysteine_counts[item])/len(aichi_cys.keys()) >= 0.9:
        aichi_cysteine_freqs.append([item, float(aichi_cysteine_counts[item])/len(aichi_cys.keys())])
aichi_cysteine_freqs
# aichi_cys


# In[11]:

plt.figure(1, figsize=(12, 4))
plt.axes()
plt.title('H5 Possible NLG Position Frequency (Aichi Reference)')
plt.xlabel('Position on Reference', labelpad=50)
plt.ylabel('Possible NLG Frequency')
plt.bar([row[0] for row in aichi_pos_freqs], [row[1] for row in aichi_pos_freqs], width=5, color='blue')
for cys in aichi_cysteine_freqs:
    plt.annotate('C'+str(cys[0]+1), xy=(cys[0]+1, 0), xytext=(cys[0]+1,-0.2),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.25, width=1, headwidth=5))
plt.annotate('Figure 2: Mapping of Possible NLG Frequency on Aichi Reference Sequence with 95% Conserved Cysteines as Reference Position Markers', xy=(0,0), xytext=(0,-0.4), size='small')
plt.show()
# for even in range(0, len(aichi_cys), 2):
#     plt.annotate('C'+str(aichi_cys[even][0]+1), xy=(aichi_cys[even][0]+1, 0.03), xytext=(aichi_cys[even][0]+1,-0.2),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.25, width=1, headwidth=5))
# for odd in range(1, len(aichi_cys), 2):
#     plt.annotate('C'+str(aichi_cys[odd][0]+1), xy=(aichi_cys[odd][0]+1, 0), xytext=(aichi_cys[odd][0]+1,-0.1),rotation='vertical',size='xx-small',arrowprops=dict(facecolor='black', shrink=0.2, frac=0.5, width=1, headwidth=5))


# In[13]:

freq_file = open('H5_Possible_NLG_Frequencies_All.csv','wb')
wr = csv.writer(freq_file, dialect='excel')
wr.writerow(['Longest H5 Sequence'])
for freq in pos_freqs:
    wr.writerow([str(freq[0]), str(freq[1]), str(ref_dict[freq[0]])])
wr.writerow(['Aichi H3 Reference'])
for aichi_freq in aichi_pos_freqs:
    wr.writerow([str(aichi_freq[0]), str(freq[1]), str(ref_dict[freq[0]])])


# In[12]:



