"""
Author: Justin Zabilansky
Date: 7/17 July 2014

Purpose of this set of custom functions:

To find the best match for a small query string in a larger reference 
string and locate it position

"""


from fuzzywuzzy import fuzz
from Levenshtein import *
from Levenshtein import distance
import re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter



def get_reference_sequence(sequences):
    '''
    get_reference_sequence takes in a txt file composed of accession numbers and sequences and parses them in to a a list of sequence lengths, accession numbers, and sequences. It then returns the longest sequence in the txt file to be used as the reference sequence to which motifs will be mapped.

    Input:
    -   sequences: txt file in fasta format composed of accession numbers and their corresponding sequences

    Output:
    -   max(seqs)[2]: a string of the longest sequence found within the set of sequences given.

    '''
    handle = open(sequences, 'rU')
    seqs = []
    for record in SeqIO.parse(handle, 'fasta'):
        seqs.append([len(record.seq), record.id, str(record.seq)])
    return max(seqs)[2]



def find_motif_nmer(sequences, motif, motif_length, edge_length):
    '''
    find_motif_nmer finds a motif within multiple sequences and returns the start of the motif in the sequence, a nmer with a specified number of of residues bordering the motif on each side (edge_length), the accession number of the sequence the motif was found in, and the sequence

    Input:
    -   sequences: txt file of sequence records in fasta format
    -   motif: set of characters in agreement with regular expression syntax
    -   motif_length: number of characters in motif
    -   edge_length: number of characters to be included on edge of motif when
        creating an nmer that will be searched for in reference sequence

    Output:
    -   motifs: a list of positons of motif matches, nmers containing motifs, accession numbers of sequence in which motif was found, and full sequence records (Seq instance)
    '''
    handle = open(sequences, 'rU')
    accessions_seqs = []
    for record in SeqIO.parse(handle, 'fasta'):
        #
        record_accession = record.id
        accessions_seqs.append([record_accession, record.seq])
    
    motifs = {rec[0]: {match.start(): str(rec[1][match.start()-edge_length:match.start()+motif_length+edge_length]) for match in \
    re.finditer(motif, str(rec[1]))} for rec in accessions_seqs}

    # for record in accessions_seqs:
    #     motif_match = re.finditer(motif, str(record[1]))
    #     for match in motif_match:
    # motifs.append([match.start(), record[1][match.start()-edge_length:match.start()+motif_length+edge_length], 
                # record[0], record[1]])
    return motifs



def overlapping_nmer(string, n):
    """
    This method takes in a string and breaks it into overlapping n-mer chunks as a dictionary.
    
    Input:
    -   string: a string
    -   n: the length of th n-mer
    
    Output:
    -   string_dict: a dictionary that of the structure {'string1':position, 'string2':position+1 etc...}
    """
    string_dict = {}
    for i,letter in enumerate(string):
        if i < len(string)-(n-1):
           string_dict[i] = string[i:i+n]
    return string_dict

def reverse_dictionary(dictionary):
    new_dictionary = {}
    for key, value in dictionary.items():
        new_dictionary[value] = key
    return new_dictionary
    
def best_nmers(sequences, motif, motif_length, edge_length, padding=0,score_buffer=0, position_buffer=0):
    '''
    best_nmers takes in a txt file containing accession numbers, finds the nmers in which the motif of interest is contained, finds the longest sequence to be used as the reference sequence which it turns into a dictionary with keys being overlapping nmers and positions as values. The function then creates a list of lists for each motif nmer, reference nmer pair, scores the pairs using Levenshtein distance and returns a list of the pairs with the lowest Levenshtein scores that are within a range of socre values and a range of positions close to the original position of the motif nmer.

    Inputs:
    -   sequences: txt file containing the accession numbers and sequences to be analyzed
    -   motif: sequence of interest to be mapped to reference sequence in accordance with regular expressions system of creating sets 'N[^P][ST]'
    -   motif_length: number of characters a motif being searched for will contain
    -   edge_length: number of positions to be included on each side of the motif when creating nmers containing the motif
    -   padding: extra positions to be considered when looking for matches within reference sequence in order to accound for in-dels
    -   score_buffer: allowed deviation from minimum Levenshtein distance when determining best matches for nmer in reference sequence
    -   position_buffer: number of positions a reference nmer match is allowed to be away from the original position of the motif nmer in the non-reference sequence

    Outputs:
    -   min_list: a list of motif nmer, reference nmer pairs with the form     [score, accession number, motif nmer, motif position, reference nmer, reference nmer position]
    '''

    # querys = []
    # for motif_match in find_motif_nmer(sequences, motif, motif_length, edge_length):
    #     querys.append([str(motif_match[1]), motif_match[0], motif_match[2]])
    ref_seq = get_reference_sequence(sequences)
    ref_dict = overlapping_nmer(ref_seq, motif_length + 2*edge_length +padding)
    reverse_ref_dict = reverse_dictionary(ref_dict)
    half_buffer = int(0.5*position_buffer)

    motifs = find_motif_nmer(sequences, motif, motif_length, edge_length)


    best_nmer_matches = {}
    for accession in motifs:
        best_nmer_matches[accession] = {item[0]: min([score for score in \
            [(distance(item[1],refmer), reverse_ref_dict[refmer], refmer, item[1]) for refmer \
            in [ref_dict[key] for key in \
            ref_dict.keys()[item[0]-half_buffer:item[0]+half_buffer]]]]) \
        for item in motifs[accession].items()}

    return best_nmer_matches

# def better_nmers(sequences, motif, motif_length, edge_length, padding=0, score_buffer=0, position_buffer=0):

#     ref_seq = get_reference_sequence(sequences)
#     ref_dict = overlapping_nmer(ref_seq, motif_length+2*(edge_length+padding))

#     motifs = find_motif_nmer(sequences, motif, motif_length, edge_length)

#     def lowest_score(motifs, ref_dict):
#         lowest_scores = {}
#         for accession in motifs.keys():
#             for nmer in motifs[accession]:



def get_cysteines(sequences, edge_length, score_buffer, position_buffer):
    ref_seq = get_reference_sequence(sequences)
    ref_dict = overlapping_nmer(ref_seq, 2*edge_length+1)
    reverse_ref_dict = reverse_dictionary(ref_dict)
    half_buffer = int(position_buffer/2)
    cysteines = find_motif_nmer(sequences, 'C', 1, edge_length)
    best_cysteines = {}
    for accession in cysteines:
        best_cysteines[accession] = {item[0]: [score for score in \
        [(distance(item[1], refmer), reverse_ref_dict[refmer], refmer, item[1]) for refmer in [ref_dict[key] for key in ref_dict.keys()[item[0] - half_buffer: item[0] + half_buffer]] if refmer[5] == 'C'] if score[0] <= score_buffer] for item in cysteines[accession].items() if len(item[1]) == 2*edge_length+1} 
    
    for accession in best_cysteines.keys():
        for pos in best_cysteines[accession].keys():
            if best_cysteines[accession][pos] == []:
                best_cysteines[accession].pop(pos)

    return best_cysteines

# def find_receptor_binders(sequences, edge_length, score_buffer=1, position_buffer):
    # ref_seq = get_reference_sequence(sequences)
    # ref_dict = overlapping_nmer(ref_seq, 2*edge_length+1)
    # reverse_ref_dict = reverse_dictionary(ref_dict)
    # half_buffer = int(position_buffer/2)
    # receptor_dict = {190: 'H......[DE]....Y', 220: 'K...DQ'}



    # for query in querys:
    #     for key in ref_dictionary.keys()[int(query[1]-half_buffer):int(query[1]+half_buffer)]:
    #         nmer_scores.append([distance(str(ref_dictionary[key]), \
    #             str(query[0])), query[2], query[0], query[1], \
    #         ref_dictionary[key], key])
    
    # min_list = []
    # for nmer_score in nmer_scores:
    #     if (nmer_score[0] <= min(nmer_scores)[0] + score_buffer):
    #         min_list.append(nmer_score)

    # min_list.sort()
    # zeros_list = []
    # for min_dist in min_list:
    #     if min_dist[0] == 0:
    #         zeros_list.append(min_dist)
    #     elif min_dist[0] != 0:
    #         if min_dist[5] in 


    # return min_list




    # for query in querys:
    #     for ref_nmer in ref_dictionary:

    #         nmer_scores.append([distance(ref_nmer, str(query[0])), query[2], query[0], query[1], ref_nmer, ref_dictionary[ref_nmer] + edge_length])  

    # min_list = []
    # for score in nmer_scores:
    #     if (score[0] <= min(nmer_scores)[0] + score_buffer):
    #         if (abs(score[3] - score[5]) <= position_buffer): 
    #             min_list.append(score)
    # return min_list    



# def get_all_motif_positions(min_list):
#     '''
#     get_all_motif_positions takes in the output of best_nmers (min_list) and returns a dictionary with keys = accession numbers, and values being lists of tuples with form (motif match nmer, reference position)

#     Input:
#     -   min_list: output of best_nmers which is a list of lists in format [[Levenshtein Distance, accession number, motif match nmer, motif match position, reference nmer match, reference nmer position]]

#     Output:
#     -   all_motif_positions: dictionary in form {Accession Number 1:[(Motif Match Nmer 1, Reference Position 1), ...], ...}
#     '''
#     all_motif_positions = {} 
#     for mini_access in min_list:
#         if mini_access[1] not in all_motif_positions.keys():
#             all_motif_positions[mini_access[1]] = {}
    
#     for mini_match in  min_list:
#         if mini_match[3] not in all_motif_positions[mini_match[1]].keys():
#             all_motif_positions[mini_match[1]][mini_match[3]] = mini_match[2]

#     # for accession_num in all_motif_positions.keys():
#     #     for motif_position in all_motif_positions[accession_num]:



#     return all_motif_positions



def motif_frequency(best_nmer_matches):
    '''
    motif_frequency takes in a dictionary of Accessions with their corresponding motif Nmers and Positions and returns lists of the frequency with which the motif occurs at this position and the frequency with which the nmer occurs

    Inputs:
    -   all_motif_positions: dictionary of accession numbers with corresponding nmers and positions at which they occurs

    Outputs:
    -   pos_freqs: list of tuples of positions at which the motif occurs and the frequency at which it is present at this position within the set of accession numbers
    '''
    positions = []
    for accession in best_nmer_matches.values():
        for scored_match in accession.values():
            positions.append(scored_match[1])

    pos_counts = Counter(positions)

    pos_freqs = []
    for item in pos_counts.items():
        pos_freqs.append([item[0], float(item[1])/len(best_nmer_matches.values())])

    return pos_freqs
    # positions = []
    # nmers = []
    # for accession in all_motif_positions.keys():
    #     for position in all_motif_positions[accession].keys():
    #         nmers.append(all_motif_positions[accession][position])
    #         positions.append(position)
    # pos_counts = Counter(positions)
    # nmer_counts = Counter(nmers)
    
    # pos_freqs = []
    # for item in pos_counts.items():
    #     pos_freqs.append([item[0],float(item[1])/len(all_motif_positions)])

    # nmer_freqs = []
    # for item in nmer_counts.items():
    #     nmer_freqs.append([item[0],float(item[1])/len(all_motif_positions)])

    # return pos_freqs, nmer_freqs


    # min_list = [[query_nmer, reference_nmer, reference_position] for score, query_nmer, query_nmer_postion, reference_nmer, reference_position in nmer_scores if score <= (min(nmer_scores) + score_buffer) and abs(query_nmer_position - reference_position) < position_buffer]
    
# def compare_to_cysteines(best_nmers, sequences): 
#     cysteines = find_motif_nmer(sequences, 'C', 1, 5)
#     for accession in best_nmers.keys():

    




# def best_nmer(ref_str, query_list, padding=0):
#     """
#     best_nmer takes in a reference string (nmer from motifs) and finds the strings where it 
#     best matches
    
#     Input:
#     -   ref_str: string to be searched
#     -   query_str: string that is being searched for
#     -   padding: make nmers being searched larger than query string if desired
    
#     Output:
#     -   tuple with the following structure (list of best string matches, best
#         score, dictionary of reference nmers) 
#     """
    
#     nmers = []
#     ref_dictionary = overlapping_nmer(ref_str, len(query_str)+padding)
#     for ref_nmer in ref_dictionary:
#         nmers.append([distance(ref_nmer,query_str), ref_nmer])   
#     min_list = [item for score, item in nmers if [score, item] \
#     <= min(nmers)]
#     min_score = min(nmers)[0]
#     # print sorted(nmers)


    # return (min_list, min_score, ref_dictionary)
    
# nmers.append([fuzz.partial_ratio(ref_nmer, query_str), ref_nmer])
        # nmers.append([fuzz.ratio(ref_nmer, query_str), ref_nmer])
        # score = fuzz.partial_ratio(ref_nmer, query_str)
        # if score > min_score:
            # min_list = []
            # min_list.append(ref_nmer)
            # min_score = score
        # elif score == min_score:
            # min_list.append(ref_nmer)
    #min_score = max(min_list)[0]


# def find_match_position(ref_str, query_str, padding=0, shift=0, score_buffer=0, position_buffer=0):
#     """
#     Takes in a reference string and returns the position which matches the query string best
    
#     Input:
#     -   ref_str: string to be searched
#     -   query_str: string that is being searched for
#     -   padding: make nmers being matched longer than the query string
#     -   shift: add a systematic shift in position to find position of specific 
#         residue of interest in nmer
    
#     Output:
#     -   nmer_position: list of positions of residues of interest that match 
#         the query string best 
#     """
#     min_list, min_score, ref_dictionary = best_nmer(ref_str, query_str, \
#         padding=padding)
#     nmer_position = []
#     for nmer in min_list:   
#         nmer_position.append(ref_dictionary[nmer] + shift)
#     return nmer_position
        
        