import db_handling
from memory_tempfile import MemoryTempfile
from parallelization import Parallelization
import subprocess
import os

class Alignment:

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
        if protDB.verify_comparison(seq_a[0], seq_b[0]):
            f1 = tempfile.NamedTemporaryFile(delete = False)
            f1.write(seq_a[1].encode('utf-8'))
            f1.close()
            f2 = tempfile.NamedTemporaryFile(delete = False)
            f2.write(seq_b[1].encode('utf-8'))
            f2.close()
            result = subprocess.run(['needle --asequence '+f1.name+' --bsequence '+f2.name+' --gapopen 10.0 --gapextend 0.5 --stdout < enter_param'], stdout=subprocess.PIPE, shell=True)
            lines = result.stdout.decode('utf-8').splitlines()
            values = self.get_values_from_needle_align(lines)
            values.insert(0,seq_a[0])
            values.insert(0,seq_b[0])
            print(seq_a[0]+' vs '+seq_b[0])
            print(values)
            protDB.insert_needle_alignment(values)
            os.unlink(f1.name)
            os.unlink(f2.name)
            

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

    @staticmethod
    def populate_parameters_equivalent_sequences(seq_list):
        previous_entry_id = ''
        previous_seq_id = ''
        previous_equivalent_seq_id = ''

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
            previous_equivalent_seq_id = current_equivalent_seq_id

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