import db_handling
from memory_tempfile import MemoryTempfile
from parallelization import Parallelization
import subprocess
import os
import operator
from ete3 import NCBITaxa, ncbi_taxonomy
import warnings
import numpy as np
import platform
import sys
import signal
from progress.bar import Bar


class Alignment:


    # !!CONSTRUIR CLASSE ADEQUADA!!
    @staticmethod
    def get_system_type():
        system = platform.system()
        if system.find('Windows') != -1:
            return 'w'
        elif system.find('Linux') != -1:
            return 'l'
        elif system.find('CYGWIN') != -1:
            return 'c'

    comparing_ranks = ['species', 'genus', 'order', 'family', 'phylum', 'class']

    def compare_sequences_same_entry(self, seq_list):
        params_a , params_b = self.populate_parameters_same_entry_comparisons(seq_list)
        Parallelization.parallelize_7(self.perform_needle_comparison, [params_a, params_b])

    @staticmethod
    def populate_parameters_same_entry_comparisons(seq_list):
        previous_entry_id = ''
        previous_seq_id = ''
        previous_entry_value = ''

        sequence_list_a = [(None,None)]
        sequence_list_b = [(None,None)]

        for seq in seq_list:
            current_entry_id = seq[1]
            current_seq_id = seq[0]
            current_entry_value = seq[2]

            if (current_entry_id == previous_entry_id and 
                current_seq_id != previous_seq_id):

                comparison_abstence = True
                i = 0
                while comparison_abstence and i < len(sequence_list_a):
                    if ((sequence_list_a[i][0] == previous_seq_id and
                        sequence_list_b[i][0] == current_seq_id) or
                        (sequence_list_b[i][0] == previous_seq_id and
                        sequence_list_a[i][0] == current_seq_id)):

                        comparison_abstence = False

                    i += 1
                if comparison_abstence:
                    sequence_list_a.append((previous_seq_id , previous_entry_value))
                    sequence_list_b.append((current_seq_id , current_entry_value))
                    print('added', current_seq_id, previous_seq_id)
            previous_entry_id = seq[1]
            previous_seq_id = seq[0]
            previous_entry_value = seq[2]

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a , sequence_list_b

    def perform_needle_comparison(self, seq_a, seq_b):
        protDB = db_handling.ProteinDatabase()
        tempfile = MemoryTempfile(preferred_paths=['/mnt/tmp'])
        f1 = tempfile.NamedTemporaryFile(delete = False)
        f1.write(seq_a[1].encode('utf-8'))
        f1.close()
        f2 = tempfile.NamedTemporaryFile(delete = False)
        f2.write(seq_b[1].encode('utf-8'))
        f2.close()
        f3 = tempfile.NamedTemporaryFile(delete = False)
        subprocess.run(['needle --asequence '+f1.name+' --bsequence '+f2.name+' --gapopen 10.0 --gapextend 0.5 --stdout '+f3.name+' < enter_param > /dev/null 2>&1'], stdout=subprocess.PIPE, shell=True)
        lines = [line.decode('utf-8') for line in f3.readlines()]
        f3.close()
        values = self.get_values_from_needle_align(lines)
        values.insert(0,seq_a[0])
        values.insert(0,seq_b[0])
        if values[2] == 0 and values[3] == 0 and values[4] == 0 and values[5] == 0 and values[6] == 0.0:
            os.unlink(f1.name)
            os.unlink(f2.name)
            os.unlink(f3.name)
            self.perform_needle_comparison(seq_a, seq_b)
            return
        protDB.insert_needle_alignment(values)
        os.unlink(f1.name)
        os.unlink(f2.name)
        os.unlink(f3.name)

            

    @staticmethod
    def get_values_from_needle_align(lines):
        alignment_length = 0
        alignment_identity = 0
        alignment_similarity = 0
        alignment_gaps = 0
        alignment_score = 0.
        for line in lines:
            if line.startswith('#'):
                if line.find('Length') != -1:
                    alignment_length = int(line.split(':')[1].strip())
                if line.find('Identity') != -1:
                    alignment_identity = int(line.split(':')[1].split('/')[0].strip())
                if line.find('Similarity') != -1:
                    alignment_similarity = int(line.split(':')[1].split('/')[0].strip())
                if line.find('Gaps') != -1:
                    alignment_gaps = int(line.split(':')[1].split('/')[0].strip())
                if line.find('Score') != -1:
                    alignment_score = float(line.split(':')[1].strip())
        return [alignment_length, alignment_identity, alignment_similarity, alignment_gaps, alignment_score]

    def delete_redundant_comparisons(self, seq_list):
        protDB = db_handling.ProteinDatabase()
        params_a , params_b = self.populate_parameters_same_entry_comparisons(seq_list)
        Parallelization.parallelize_7(protDB.delete_redundant_comparison, [params_a, params_b])


    def assing_equivalent_seq_ids(self, seq_list):
        protDB = db_handling.ProteinDatabase()
        params_a, params_b  = self.populate_parameters_equivalent_sequences(seq_list)
        Parallelization.parallelize_7(protDB.update_biological_sequence, [params_a, params_b])

    def assing_representant_entry_ids(self, seq_list):
        protDB = db_handling.ProteinDatabase()
        params_a, params_b  = self.populate_parameters_equivalent_sequences(seq_list)
        Parallelization.parallelize_7(protDB.update_biological_sequence, [params_a, params_b])

    @staticmethod
    def populate_parameters_equivalent_sequences(seq_list):
        previous_entry_id = ''
        previous_seq_id = ''

        sequence_list_a = [{}]
        sequence_list_b = [None]

        current_entry_equivalent_seq_ids = {}

        for seq in seq_list:
            current_entry_id = seq[2]
            current_seq_id = seq[0]
            current_equivalent_seq_id = seq[1]

            if (current_entry_id == previous_entry_id and 
                current_seq_id == previous_seq_id):
                current_entry_equivalent_seq_ids[current_seq_id].append(current_equivalent_seq_id)
            elif(current_entry_id == previous_entry_id and 
                current_seq_id != previous_seq_id):              

                current_entry_equivalent_seq_ids[current_seq_id] = [current_equivalent_seq_id]
            elif current_entry_id != previous_entry_id:
                biggest = 0
                chosen_id = ''
                for key, items in current_entry_equivalent_seq_ids.items():
                    if len(items) > biggest:
                        chosen_id = key
                if chosen_id:
                    for equivalent_id in current_entry_equivalent_seq_ids[chosen_id]:
                        sequence_list_a.append({'equivalent_to':chosen_id})
                        sequence_list_b.append(equivalent_id)

                current_entry_equivalent_seq_ids = {current_seq_id:[current_equivalent_seq_id]}
        
            previous_entry_id = current_entry_id
            previous_seq_id = current_seq_id

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a, sequence_list_b

    def ellect_representative_seqs_by_equivalence(self, seq_list):
        protDB = db_handling.ProteinDatabase()
        params_a , params_b = self.populate_parameters_representative_seqs_by_equivalence(seq_list)
        Parallelization.parallelize_7(protDB.update_biological_sequence, [params_a, params_b])

    @staticmethod
    def populate_parameters_representative_seqs_by_equivalence(seq_list):
        protDB = db_handling.ProteinDatabase()
        previous_entry_id = ''
        previous_seq_id = ''

        sequence_list_a = [{}]
        sequence_list_b = [None]

        entry_equivalent_to = {}

        for seq in seq_list:
            current_entry_id = seq[1]
            current_seq_id = seq[0]
            current_seq_equivalent_to = seq[3]

            if (current_entry_id == previous_entry_id and 
                current_seq_id != previous_seq_id):
                if current_seq_equivalent_to != None:
                    if current_seq_equivalent_to in entry_equivalent_to:
                        entry_equivalent_to[current_seq_equivalent_to] += 1
                    else:
                        entry_equivalent_to[current_seq_equivalent_to] = 1
            elif current_entry_id != previous_entry_id:
                biggest_hits = 0
                ellected_seq = ''
                for seq_id, hits in entry_equivalent_to.items():
                    if hits > biggest_hits:
                        biggest_hits = hits
                        ellected_seq = seq_id
                if ellected_seq != '':
                    if protDB.translation_matches(ellected_seq):
                        sequence_list_a.append({'ellected_seq':1})
                        sequence_list_b.append(ellected_seq)
                entry_equivalent_to = {}
                if current_seq_equivalent_to != None:
                    entry_equivalent_to[current_seq_equivalent_to] = 1

            previous_entry_id = seq[1]
            previous_seq_id = seq[0]

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a , sequence_list_b


    def ellect_representative_seqs_by_abstence(self, seq_list):
        protDB = db_handling.ProteinDatabase()
        params_a , params_b = self.populate_parameters_representative_seqs_by_abstence(seq_list)
        Parallelization.parallelize_7(protDB.update_biological_sequence, [params_a, params_b])

    @staticmethod
    def populate_parameters_representative_seqs_by_abstence(seq_list):
        protDB = db_handling.ProteinDatabase()
        previous_entry_id = ''
        previous_seq_id = ''

        sequence_list_a = [{}]
        sequence_list_b = [None]

        entry_equivalent_to = {}

        for seq in seq_list:
            current_entry_id = seq[1]
            current_seq_id = seq[0]
            current_seq_value = seq[2]
            current_seq_status = seq[3]

            ellected = False

            possible_repres = {}
            
            if (current_entry_id == previous_entry_id and 
                current_seq_id != previous_seq_id):
                if not ellected and current_seq_value:
                    if protDB.translation_matches(current_seq_id):
                        sequence_list_a.append({'ellected_seq':1})
                        sequence_list_b.append(current_seq_id)
                        ellected = True
                    else:
                        possible_repres[current_seq_id] = current_seq_status
            elif current_entry_id != previous_entry_id:
                if not ellected and current_seq_value:
                    if protDB.translation_matches(current_seq_id):
                        sequence_list_a.append({'ellected_seq':1})
                        sequence_list_b.append(current_seq_id)
                        ellected = True
                    else:
                        ellected_anyways = ''
                        altered_status_seq_id = ''
                        if possible_repres:
                            for seq_id, status in possible_repres.items():
                                if not status:
                                    ellected_anyways = seq_id
                                else:
                                    altered_status_seq_id = seq_id
                            if ellected_anyways != '':
                                sequence_list_a.append({'ellected_seq':1})
                                sequence_list_b.append(ellected_anyways)
                            elif altered_status_seq_id != '':
                                sequence_list_a.append({'ellected_seq':1})
                                sequence_list_b.append(altered_status_seq_id)
                        else:
                            sequence_list_a.append({'ellected_seq':1})
                            sequence_list_b.append(current_seq_id)

                        possible_repres = {}
                        ellected = False

            previous_entry_id = seq[1]
            previous_seq_id = seq[0]

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a , sequence_list_b

    @staticmethod
    def populate_parameters_same_taxon_comparisons(seq_list, latest_seq_a=None, latest_seq_b=None, partial=False, taxon_rank=None):
        protDB = db_handling.ProteinDatabase()
        previous_entry_taxon = ''
        previous_seq_id = ''

        sequence_list_a = [(None,None)]
        sequence_list_b = [(None,None)]

        found_pair = True

        if latest_seq_a and latest_seq_b:
            found_pair = False


        taxon_id_list = []
        n = 0
        for seq in seq_list:
            current_entry_taxon = seq[2]
            current_seq_id = seq[0]
            current_entry_value = seq[1]
            if taxon_rank and taxon_rank != 'species':
                current_entry_comparison_rank = seq[3]
 
            if (current_entry_taxon == previous_entry_taxon and 
                current_seq_id != previous_seq_id):
                if taxon_rank and taxon_rank != 'species':
                    taxon_id_list.append((current_seq_id , current_entry_value, current_entry_comparison_rank))
                else:
                    taxon_id_list.append((current_seq_id , current_entry_value))
            elif current_entry_taxon != previous_entry_taxon or n == len(seq_list) - 1:

                indices = np.tril_indices(n=len(taxon_id_list),k=-1)

                for index in range(len(indices[0])):
                    if latest_seq_a and latest_seq_b and latest_seq_a == taxon_id_list[indices[0][index]][0] and latest_seq_b == taxon_id_list[indices[1][index]][0]:
                        found_pair = True
                        continue

                    if latest_seq_a and latest_seq_b and not found_pair:
                        continue
                    

                    if taxon_rank and taxon_rank != 'species':
                        if taxon_id_list[indices[0][index]][2] == taxon_id_list[indices[1][index]][2]:
                            database_abstence = protDB.verify_comparison(taxon_id_list[indices[0][index]][0],taxon_id_list[indices[1][index]][0])
                        else:
                            database_abstence  = True
                    else:
                        database_abstence  = protDB.verify_comparison(taxon_id_list[indices[0][index]][0],taxon_id_list[indices[1][index]][0])
                    
                    if database_abstence and found_pair:
                        sequence_list_a.append((taxon_id_list[indices[0][index]][0],taxon_id_list[indices[0][index]][1]))
                        sequence_list_b.append((taxon_id_list[indices[1][index]][0],taxon_id_list[indices[1][index]][1]))
                        # print('added', taxon_id_list[indices[0][index]][0], taxon_id_list[indices[1][index]][0])
                               
                taxon_id_list = []


                if partial and previous_entry_taxon != '' and sequence_list_a != [(None,None)]:
                    sequence_list_a.pop(0)
                    sequence_list_b.pop(0)

                    return sequence_list_a, sequence_list_b, seq_list[(n-1):], previous_entry_taxon

            previous_entry_taxon = seq[2]
            previous_seq_id = seq[0]

            n+=1

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a , sequence_list_b , None, None

    def compare_sequences_same_taxon(self, seq_list):
        protDB = db_handling.ProteinDatabase()

        params_a , params_b, _ = self.populate_parameters_same_taxon_comparisons(seq_list)
        Parallelization.parallelize_7(self.perform_needle_comparison, [params_a, params_b])

        protDB.assign_identity_paramters_comparisons()

    def compare_sequences_taxon_rank(self, repres_list, rank='species'):
        protDB = db_handling.ProteinDatabase()

        sorted_repres_list = self.sort_collection_by_taxon_rank(repres_list,key=2,rank=rank)

        
        if sorted_repres_list != repres_list:
            params_a , params_b, sorted_repres_list = self.populate_parameters_same_taxon_comparisons(sorted_repres_list, partial=True)
            while params_a:
                print('paralellizing needleman alignments...')

                Parallelization.parallelize_7(self.perform_needle_comparison, [params_a, params_b])

                params_a , params_b, sorted_repres_list = self.populate_parameters_same_taxon_comparisons(sorted_repres_list, partial=True)
            
            protDB.assign_identity_paramters_comparisons()


    def populate_parameters_represented_by(self, align_list, identity_parameter=6, leaf=True, rank=None):
        protDB = db_handling.ProteinDatabase()
        previous_taxon_id = ''
        previous_entry_id_a = ''
        previous_seq_id_b = ''
        previous_seq_identity_parameter = ''

        ncbi = NCBITaxa()

        sequence_list_a = [{}]
        sequence_list_b = [None]
        
        compare_matrix = {}
        lines_and_cols_sum = {}

        for align in align_list:
            current_taxon_id = align[0]
            current_entry_id_a = align[1]
            current_entry_id_b = align[2]
            current_entry_identity_parameter = align[identity_parameter]

            if current_taxon_id == previous_taxon_id:
                if current_entry_id_a in compare_matrix.keys():
                    compare_matrix[current_entry_id_a][current_entry_id_b] = current_entry_identity_parameter
                    print('compare matrix('+str(current_taxon_id)+')<-'+str(current_entry_id_a)+','+str(current_entry_id_b)+'='+str(current_entry_identity_parameter))
                elif current_entry_id_b in compare_matrix.keys():
                    compare_matrix[current_entry_id_b][current_entry_id_a] = current_entry_identity_parameter
                    print('compare matrix('+str(current_taxon_id)+')<-'+str(current_entry_id_b)+','+str(current_entry_id_a)+'='+str(current_entry_identity_parameter))

            elif current_taxon_id != previous_taxon_id:
                new_chunk_of_list_a , new_chunk_of_list_b, chosen_id, just_inserted = self.ellect_id_from_compare_matrix(compare_matrix,lines_and_cols_sum,previous_taxon_id)
                previous_chosen_ids = []
                while not new_chunk_of_list_a and chosen_id != '' and lines_and_cols_sum[chosen_id] != 0:
                    lines_and_cols_sum.pop(chosen_id)
                    print('popped', chosen_id)
                    previous_chosen_ids.append(chosen_id)
                    new_chunk_of_list_a , new_chunk_of_list_b, chosen_id, just_inserted = self.ellect_id_from_compare_matrix(compare_matrix,lines_and_cols_sum,previous_taxon_id,previous_chosen_ids,just_inserted)
                
                for item in new_chunk_of_list_a:
                    sequence_list_a.append(item)
                for item in new_chunk_of_list_b:
                    sequence_list_b.append(item)

                
                protDB.renew_conn()

                if not leaf:
                    entry_ids = protDB.get_entry_ids_representative_of_taxon(current_taxon_id)
                    sorted_entry_ids = self.sort_collection_by_taxon_rank(entry_ids, key=1, rank=rank, rank_id=current_taxon_id)
                    entry_ids = sorted_entry_ids
                else:
                    entry_ids = protDB.get_entry_ids_from_taxon(current_taxon_id)

                print(entry_ids, '- entry ids')
                
                compare_matrix, lines_and_cols_sum = self.generate_compare_matrix_and_lines_cols(entry_ids, current_taxon_id)

            previous_taxon_id = current_taxon_id
            previous_entry_id_a = current_entry_id_a
            previous_entry_id_b = current_entry_id_b
            previous_entry_identity_parameter = current_entry_identity_parameter

        sequence_list_a.pop(0)
        sequence_list_b.pop(0)

        return sequence_list_a , sequence_list_b

    @staticmethod
    def generate_compare_matrix_and_lines_cols(entry_ids,taxon_id):
        print('generating compare matrix for '+str(taxon_id)+' ('+str(len(entry_ids))+')')
        compare_matrix = {}
        lines_and_cols_sum = {}
        for i in range(len(entry_ids)):
            lines_and_cols_sum[entry_ids[i][0]] = 0
            for j in range(len(entry_ids)):
                if i > j:
                    # print('matrix i='+str(i),entry_ids[i][0],'matrix j='+str(j),entry_ids[j][0])
                    if entry_ids[i][0] in compare_matrix.keys():
                        compare_matrix[entry_ids[i][0]][entry_ids[j][0]] = None
                    else:
                        compare_matrix[entry_ids[i][0]] = {entry_ids[j][0]:None}
        return compare_matrix, lines_and_cols_sum
        

    @staticmethod
    def ellect_id_from_compare_matrix(compare_matrix,lines_and_cols_sum,taxon_id, previous_chosen_ids=[],just_inserted = []):
        sequence_list_a = []
        sequence_list_b = []
        protDB = db_handling.ProteinDatabase()
        ncbi = NCBITaxa()
        chosen_id = ''
        insert_representative = False
        for entry_a, comparisons in compare_matrix.items():
            for entry_b, value in comparisons.items():
                if value:
                    if entry_a in lines_and_cols_sum.keys():
                        lines_and_cols_sum[entry_a] += value
                    if entry_b in lines_and_cols_sum.keys():
                        lines_and_cols_sum[entry_b] += value
                    else:
                        print('couldnt fit alignment data to descendant data')
                        print(entry_a, entry_b)
            sorted_dict = dict(sorted(lines_and_cols_sum.items(),
                    key=operator.itemgetter(1),
                    reverse=True))
            print(sorted_dict)
            max_match = 0
            chosen_id = ''
            for entry_id, matches in sorted_dict.items():
                if matches > max_match:
                    max_match = matches
                    chosen_id = entry_id
            
            if chosen_id != '':
                if chosen_id in compare_matrix.keys() and compare_matrix[chosen_id] != None:
                    for entry_id, value in compare_matrix[chosen_id].items():
                        entry = protDB.get_protein_entry(entry_id)
                        if value and entry[4] and entry[4] != chosen_id and entry[4] not in previous_chosen_ids and (chosen_id, entry_id) not in just_inserted:
                            print(entry_id, 'represented by', chosen_id)
                            sequence_list_a.append({'represented_by':chosen_id})
                            sequence_list_b.append(entry_id)
                            insert_representative = True
                            just_inserted.append((chosen_id , entry_id))
                        elif value and not entry[4] and (chosen_id, entry_id) not in previous_chosen_ids:
                            print(entry_id, 'represented by', chosen_id)
                            sequence_list_a.append({'represented_by':chosen_id})
                            sequence_list_b.append(entry_id)
                            insert_representative = True
                            just_inserted.append((chosen_id , entry_id))
                        else:
                            print('REPETITIVE REDUNDANT CYCLING')
                    for entry_id, comparisons_entry in compare_matrix.items():
                        if chosen_id in comparisons_entry.keys() and comparisons_entry[chosen_id] != None:
                            entry = protDB.get_protein_entry(entry_id)
                            if entry[4] and entry[4] != chosen_id and (chosen_id, entry_id) not in just_inserted:
                                print(entry_id, 'represented by', chosen_id)
                                sequence_list_a.append({'represented_by':chosen_id})
                                sequence_list_b.append(entry_id)
                                insert_representative = True
                                just_inserted.append((chosen_id , entry_id))
                            elif not entry[4] and entry_id not in previous_chosen_ids:
                                print(entry_id, 'represented by', chosen_id)
                                sequence_list_a.append({'represented_by':chosen_id})
                                sequence_list_b.append(entry_id)
                                insert_representative = True
                                just_inserted.append((chosen_id , entry_id))
                            else:
                                print('REPETITIVE REDUNDANT CYCLING')                    
                else:
                    for entry_id, comparisons_entry in compare_matrix.items():
                        if chosen_id in comparisons_entry.keys() and comparisons_entry[chosen_id] != None:
                            entry = protDB.get_protein_entry(entry_id)
                            if entry[4] and entry[4] != chosen_id and (chosen_id, entry_id) not in just_inserted:
                                print(entry_id, 'represented by', chosen_id)
                                sequence_list_a.append({'represented_by':chosen_id})
                                sequence_list_b.append(entry_id)
                                insert_representative = True
                                just_inserted.append((chosen_id , entry_id))
                            elif not entry[4] and entry_id not in previous_chosen_ids:
                                print(entry_id, 'represented by', chosen_id)
                                sequence_list_a.append({'represented_by':chosen_id})
                                sequence_list_b.append(entry_id)
                                insert_representative = True
                                just_inserted.append((chosen_id , entry_id))
                            else:
                                print('REPETITIVE REDUNDANT CYCLING') 

            if chosen_id != '' and insert_representative:
                print(chosen_id, 'representative_of_taxon', taxon_id)
                sequence_list_a.append({'representative_of_taxon':taxon_id})
                sequence_list_b.append(chosen_id)
                lineage = ncbi.get_lineage(taxon_id)
                lineage_ranks = ncbi.get_rank(lineage)
                lineage_translation = ncbi.get_taxid_translator(lineage)
                print(chosen_id, 'taxon_name_representative', lineage_translation[taxon_id])
                sequence_list_a.append({'taxon_name_representative':lineage_translation[taxon_id]})
                sequence_list_b.append(chosen_id)
                print(chosen_id, 'representative_taxon_rank', lineage_ranks[taxon_id])
                sequence_list_a.append({'representative_taxon_rank':lineage_ranks[taxon_id]})
                sequence_list_b.append(chosen_id)
        return sequence_list_a, sequence_list_b, chosen_id, just_inserted

    
    
    def ellect_representative_entries_of_taxon_level(self, entry_list, leaf=True,rank=None):
        protDB = db_handling.ProteinDatabase()
        params_a , params_b = self.populate_parameters_represented_by(entry_list, leaf=leaf,rank=rank)
        Parallelization.parallelize_7(protDB.update_protein_entry, [params_a, params_b])

        return params_b

    def ellect_representative_entries_of_taxon(self):
        protDB = db_handling.ProteinDatabase()
        align_list = protDB.get_comparisons_same_taxon()
        params_a, params_b = self.populate_parameters_represented_by(align_list)
        while params_b:
            self.ellect_representative_entries_of_taxon_level(align_list)

            protDB.renew_conn()

            align_list = protDB.get_comparisons_same_taxon()

            params_a, params_b = self.populate_parameters_represented_by(align_list)


    @staticmethod
    def sort_collection_by_taxon_rank(collection, key, rank='species', rank_id=None):
        ncbi = NCBITaxa()

        new_collection = collection

        i=0
        for item in collection:
            lineage = ncbi.get_lineage(item[key])
            lineage_ranks = ncbi.get_rank(lineage)
            if (rank in lineage_ranks.values() and item[key] in lineage_ranks.keys()
            and lineage_ranks[item[key]] != rank):
                    for taxon, taxon_rank in lineage_ranks.items():
                            if rank == taxon_rank:
                                new_rep = []
                                for k in range(len(item)):
                                    if k == key:
                                        new_rep.append(taxon)
                                    else:
                                        new_rep.append(item[k])
                                new_rep = tuple(new_rep)
                                if rank_id == None:
                                    new_collection[i] = new_rep
                                elif taxon == rank_id:
                                    new_collection[i] = new_rep

            i+=1
        sorted_collection = tuple(sorted(new_collection,
                        key=operator.itemgetter(key)))

        return sorted_collection
        
    def ellect_representative_entries_by_rank(self, align_list, rank='species'):
        protDB = db_handling.ProteinDatabase()

        sorted_align_list = self.sort_collection_by_taxon_rank(align_list,key=0,rank=rank)
        params_b = True
        while params_b:
            params_b = self.ellect_representative_entries_of_taxon_level(sorted_align_list, leaf=False, rank=rank)

            protDB.renew_conn()

            align_list = protDB.get_comparisons_same_taxon(leaf=False)

            sorted_align_list = self.sort_collection_by_taxon_rank(align_list,key=0,rank=rank)

    @staticmethod
    def correct_erroneous_repres_of_taxon_instances():
        protDB = db_handling.ProteinDatabase()
        ncbi = NCBITaxa()

        erroneous_instances = protDB.get_erroneous_repres_of_taxon_instances()

        for instance in erroneous_instances:
            lineage = ncbi.get_lineage(instance[1])
            lineage_ranks = ncbi.get_rank(lineage)
            representative_id = protDB.get_protein_entry(instance[2])[5]
            lineage_translation = ncbi.get_taxid_translator(lineage)
            if protDB.get_protein_entry(instance[2]):
                if representative_id != None:
                    representative_id = protDB.get_protein_entry(instance[2])[5]
                    print('representative_id', representative_id)
                    lineage_representative = ncbi.get_lineage(representative_id)
                    print('lineage_representative', lineage_representative)
                    lineage_representative_ranks = ncbi.get_rank(lineage_representative)
                    print('lineage_representative_ranks', lineage_representative_ranks)
                    if lineage_ranks[instance[1]] == lineage_representative_ranks[representative_id]:
                        count_instance = protDB.get_count_of_children_of_repres(instance[1])
                        count_representative = protDB.get_count_of_children_of_repres(representative_id)
                        if count_representative >= count_instance:
                            protDB.update_protein_entry({'representative_of_taxon':None},instance[1])
                        else:
                            protDB.update_protein_entry_by_repres_by({'represented_by':instance[1]},representative_id)
                            protDB.update_protein_entry({'represented_by':instance[1]},representative_id)
                            protDB.update_protein_entry({'representative_of_taxon':None},representative_id)
            else:
                protDB.update_protein_entry({'representative_of_taxon':None},instance[0])
                protDB.update_protein_entry({'representative_of_taxon':instance[1]},instance[2])
                protDB.update_protein_entry({'taxon_name_representative':lineage_translation[instance[1]]},representative_id) 
                protDB.update_protein_entry({'representative_taxon_rank':lineage_ranks[instance[1]]},representative_id)
        
        print('---------------------------------------------------------------------------------')
    

    def compare_sequences_same_taxon_level(self, taxon_rank, current_state):
        protDB = db_handling.ProteinDatabase()
        paralell = Parallelization()
        ncbi = ncbi_taxonomy.NCBITaxa()
        
        latest_taxon_id = current_state[3]

        latest_seq_a = current_state[7]

        latest_seq_b = current_state[8]

        
        taxon_rank_sequences = protDB.get_sequences_same_taxon_rank(taxon_rank)

        max_len = len(taxon_rank_sequences)

        start = 0

        if latest_taxon_id:
            for seq in taxon_rank_sequences:
                if seq[2] == latest_taxon_id:
                    start = taxon_rank_sequences.index(seq)
                    break
            taxon_rank_sequences = taxon_rank_sequences[start:]
        
        print('ctrl + c will stop the process at appropriate timing')
        sys.stdout.write("\n")

        bar = Bar(taxon_rank+' comparisons:', max=max_len)

        bar.index = start

        first_line = True

        params_a , params_b, sorted_repres_list, taxon_id = self.populate_parameters_same_taxon_comparisons(taxon_rank_sequences, latest_seq_a, latest_seq_b, partial=True, taxon_rank=taxon_rank)
        while params_a:
            taxon_name = ncbi.get_taxid_translator([taxon_id])[taxon_id]

            sys.stdout.write("\n")
            sys.stdout.write("\033[K")
            bar_taxon = Bar(taxon_rank+' '+taxon_name, max=len(params_a))
            paralell.parallelize_7(self.perform_needle_comparison, [params_a, params_b], log_param=taxon_rank, bar=bar_taxon, handle_signal=True)
            bar_taxon.finish()

            sys.stdout.write("\033[F")
            sys.stdout.write("\033[F")
            sys.stdout.write("\033[K")

            bar.next()

            protDB.update_state({'comparisons_latest_taxon_id': taxon_id}, taxon_rank)

            params_a , params_b, sorted_repres_list, taxon_id = self.populate_parameters_same_taxon_comparisons(sorted_repres_list, partial=True, taxon_rank=taxon_rank)

        bar.finish()
        
        protDB.assign_identity_paramters_comparisons()
        protDB.update_state({'comparisons_done':1},taxon_rank)


    def ellect_taxon_rank_representatives_step(self, taxon_rank_id, taxon_rank):
        protDB = db_handling.ProteinDatabase()
        ncbi = NCBITaxa()
        paralell = Parallelization()

                        
        taxon_name = ncbi.get_taxid_translator([taxon_rank_id])[taxon_rank_id]

        representative_seqs = protDB.get_comparisons_same_taxon_id(taxon_rank_id, taxon_rank)


        identity_matrix = {}

        for seq_tuple in representative_seqs:
            if seq_tuple[0] not in identity_matrix.keys():
                identity_matrix[seq_tuple[0]] = [seq_tuple[1]]
            elif seq_tuple[0] in identity_matrix.keys():
                identity_matrix[seq_tuple[0]].append(seq_tuple[1])
            
            if seq_tuple[1] not in identity_matrix.keys():
                identity_matrix[seq_tuple[1]] = [seq_tuple[0]]
            elif seq_tuple[1] in identity_matrix.keys():
                identity_matrix[seq_tuple[1]].append(seq_tuple[0])

        if identity_matrix:
            params_a, params_b = [], []

            chosen_item = max(identity_matrix, key=lambda k: len(identity_matrix[k]))
            # print(chosen_item, identity_matrix[chosen_item])


            for protein_id in identity_matrix[chosen_item]:
                # protDB.update_protein_entry({'represented_by':chosen_item},protein_id)
                params_a.append({'represented_by': chosen_item})
                params_b.append(protein_id)
            # print('paralellizing updates...')
            
            bar = Bar(taxon_rank+' '+taxon_name, max=len(params_a))
            
            paralell.parallelize_7(protDB.update_protein_entry, [params_a, params_b], bar=bar)
            
            bar.finish()
            
 

            # print('collapsed '+str(len(identity_matrix[chosen_item]))+' seqs of '+taxon_rank+' '+taxon_name+' into protein entry '+str(chosen_item))
            protDB.update_protein_entry({'representative_of_taxon':taxon_rank_id, 'representative_taxon_rank':taxon_rank, 'taxon_name_representative':taxon_name},chosen_item)

            return True
        else:
            sys.stdout.write("\033[K")
            print(taxon_rank+' '+taxon_name)

            return False

    stop_process = False

    def signal_handler(self, signal, frame):
        self.stop_process = True

    def iterate_ellection_taxon_rank_representatives(self, current_state, taxon_rank='species'):
        signal.signal(signal.SIGINT, self.signal_handler)

        protDB = db_handling.ProteinDatabase()
 
        latest_taxon_id = current_state[5]

        taxon_id_list = protDB.get_taxon_ids_by_rank(taxon_rank)

        

        taxon_rank_id_list = []
        for taxon_id_tuple in taxon_id_list:
            taxon_rank_id_list.append(taxon_id_tuple[0])
        
        start = 0

        if latest_taxon_id:
            start = taxon_rank_id_list.index(latest_taxon_id)
            taxon_rank_id_list = taxon_rank_id_list[start:]

        print('ctrl + c will stop the process at appropriate timing')
        sys.stdout.write("\n")    


        bar = Bar(taxon_rank+' classifications:', max=len(taxon_id_list))
        bar.index = start

        for taxon_rank_id in taxon_rank_id_list:
            bar.next()
            
            not_depleted = True
            while not_depleted:
                sys.stdout.write("\n")
                not_depleted = self.ellect_taxon_rank_representatives_step(taxon_rank_id, taxon_rank)
                sys.stdout.write("\033[F")
                sys.stdout.write("\033[F")
                if self.stop_process:
                    sys.exit('\nprocess terminated correctly')

                
            protDB.update_representatives_taxon_id_by_abstence(taxon_rank_id, taxon_rank)
            protDB.update_state({'classifications_latest_taxon_id': taxon_rank_id},taxon_rank),
          
        bar.finish()




    @staticmethod
    def previous_taxon_rank(taxon_rank):
        if taxon_rank == 'species':
            return None
        elif taxon_rank == 'genus':
            return 'species'
        elif taxon_rank == 'order':
            return 'genus'
        elif taxon_rank == 'family':
            return 'order'
        elif taxon_rank == 'phylum':
            return 'family'

    @staticmethod
    def bigger_than_rank_taxon(taxon_rank, rank):
        if taxon_rank == 'species':
            if rank == 'species':
                return False
            elif rank == 'genus':
                return False
            elif rank == 'order':
                return False
            elif rank == 'family':
                return False
            elif rank == 'phylum':
                return False
            else:
                return False
        elif taxon_rank == 'genus':
            if rank == 'species':
                return True
            elif rank == 'genus':
                return False
            elif rank == 'order':
                return False
            elif rank == 'family':
                return False
            elif rank == 'phylum':
                return False
            else:
                return False
        elif taxon_rank == 'order':
            if rank == 'species':
                return True
            elif rank == 'genus':
                return True
            elif rank == 'order':
                return False
            elif rank == 'family':
                return False
            elif rank == 'phylum':
                return False
            else:
                return False
        elif taxon_rank == 'family':
            if rank == 'species':
                return True
            elif rank == 'genus':
                return True
            elif rank == 'order':
                return True
            elif rank == 'family':
                return False
            elif rank == 'phylum':
                return False
            else:
                return False
        elif taxon_rank == 'phylum':
            if rank == 'species':
                return True
            elif rank == 'genus':
                return True
            elif rank == 'order':
                return True
            elif rank == 'family':
                return True
            elif rank == 'phylum':
                return False
            else:
                return False
        elif taxon_rank == 'no rank':
            return False


    def assign_gen_ranks(self, affected_entries):
        protDB = db_handling.ProteinDatabase()
        ncbi = NCBITaxa()
        paralell = Parallelization()
        params_a, params_b = [], []

        
        for entry in affected_entries:
            for rank in self.comparing_ranks:
                taxon_id = entry[1]
                entry_id = entry[0]
                lineage = ncbi.get_lineage(taxon_id)
                lineage_ranks = ncbi.get_rank(lineage)
                if (rank in lineage_ranks.values()):
                    rank_id = list(lineage_ranks.keys())[list(lineage_ranks.values()).index(rank)]
                    params_a.append({rank+'_rank_id': rank_id})
                    params_b.append(entry_id)
        
        paralell.parallelize_7(protDB.update_protein_entry, [params_a, params_b])

        


    def assign_rank_representation(self, rank='species'):
        protDB = db_handling.ProteinDatabase()
        ncbi = NCBITaxa()
        entries_no_representative = protDB.get_entries_no_representative()

        for entry in entries_no_representative:
            taxon_id = entry[1]
            with warnings.catch_warnings(record=True) as w:
                warn_msg = None
                warnings.simplefilter("always")
                lineage = ncbi.get_lineage(taxon_id)
                for a in w:
                    warn_msg = a.message
                if warn_msg:
                    warn_data = str(warn_msg).split()
                    taxon_id = int(warn_data[-1])
                    protDB.update_protein_entry({'representative_of_taxon':taxon_id},entry[0])

            lineage_ranks = ncbi.get_rank(lineage)
            lineage_translation = ncbi.get_taxid_translator(lineage)
            insert = True
            ellected_rank_id = ''
            for rank_id, lineage_rank in lineage_ranks.items():
                if rank == lineage_rank:
                    ellected_rank_id = rank_id
            print(entry[0])
            if lineage_ranks[taxon_id] != rank:
                if not self.bigger_than_rank_taxon(lineage_ranks[taxon_id], rank) and ellected_rank_id != '':
                    protDB.update_protein_entry({'representative_of_taxon':ellected_rank_id,
                    'representative_taxon_rank':rank, 'taxon_name_representative':lineage_translation[ellected_rank_id]},entry[0])
                    insert = False
            if entry[2] == None and insert and ellected_rank_id != '':
                protDB.update_protein_entry({'representative_of_taxon':ellected_rank_id,
                    'representative_taxon_rank':rank, 'taxon_name_representative':lineage_translation[ellected_rank_id]},entry[0])


