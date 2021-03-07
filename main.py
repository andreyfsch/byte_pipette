from accession import Accession
from xml_handling import UniprotXMLHandler
from ete3 import NCBITaxa
import db_handling
import alignment
import operator
import numpy as np
import os
import sys

#Accession.download_uniprot_xml()

# uniprotHandler = UniprotXMLHandler(iterate=True)

ncbi = NCBITaxa()

align = alignment.Alignment()

protDB = db_handling.ProteinDatabase()

filtering_states = protDB.get_filtering_states()

for state in filtering_states:
        current_state = protDB.get_current_state()
        rank = state[0]
        if rank == current_state[0]:
                if current_state[2]:
                        if not current_state[4]:
                                align.iterate_ellection_taxon_rank_representatives(taxon_rank=rank, current_state=current_state)

                                protDB.set_next_filtering_state()

                                protDB.update_state({'classifications_done': 1, 'current_task':0},rank)
                        else:
                                protDB.set_next_filtering_state()
                else:
                        align.compare_sequences_same_taxon_level(taxon_rank=rank, current_state=current_state)

                        align.iterate_ellection_taxon_rank_representatives(taxon_rank=rank, current_state=current_state)

                        protDB.set_next_filtering_state()

                        protDB.update_state({'classifications_done': 1, 'current_task':0},rank)



# rank_list = protDB.get_rank_list('genus')

# print(len(rank_list))

# align.assign_gen_ranks(representative_entries)

# seqs = protDB.get_sequences_same_taxon()


# align_list = protDB.get_comparisons_same_taxon(leaf=False)

# align.ellect_representative_entries_by_rank(align_list)


#protDB.assign_identity_paramters_comparisons()

# entries_same_taxon_repres = protDB.get_sequences_same_repres()

# align.compare_sequences_taxon_rank(entries_same_taxon_repres, rank='genus')

# align_list = protDB.get_comparisons_same_taxon(leaf=False)

# align.ellect_representative_entries_by_rank(align_list, rank='genus')

# align.assign_rank_representation(rank='genus')

# align.correct_erroneous_repres_of_taxon_instances()


# protDB.set_sequences_values()

# protDB.fix_sequences_multiple_repres()

# acc = Accession()

# protDB.set_sequences_values()

# protDB.fix_sequences_multiple_repres()


# align.compare_sequences_same_entry(repres_seqs)


# for entry, data in uniprotHandler.entries:
        # if data[i][0] == 'taxon_id':
        #     taxid2name = ncbi.get_taxid_translator([data[i][1]])
        #     lineage = ncbi.get_lineage(data[i][1])
        #     print(taxid2name)
        #     print(lineage)

# print('------------ ccds ------------')
# for entry_ccds_ids in entries_ccds_ids:
#     entry_accession = entry_ccds_ids[0]
#     print("entry "+str(entry_accession)+"\n--------")
#     for entry, ccds_ids in entry_ccds_ids:
#         for id in ccds_ids:
#             print(str(id))
#             ccds_seq = UniprotXMLHandler.get_ccds_seq(Accession.get_html_ccds(id))
#             #print(ccds_seq)
#             print(str(len(ccds_seq)))

# print('------------ embl ------------')
# for entry_embl_ids in entries_embl_ids:
#     entry_accession = entry_embl_ids[0]
#     print("entry "+str(entry_accession)+"\n--------")
#     for entry, embl_ids in entry_embl_ids:
#         for id in embl_ids:
#             print(str(id))
#             embl_fasta = Accession.get_embl_seq(id)
#             embl_fasta_split_lines = embl_fasta.split('\n')
#             seq_len = 0
#             for line in embl_fasta_split_lines:
#                 if line.startswith('>'):
#                     pass
#                 else:
#                     seq_len += len(line)
#             print(str(seq_len))
            
# print('------------ taxonomy ------------')

# for entry_accession, entry_taxon_id in entries_taxon_ids:
#     print("entry "+str(entry_accession)+"\n--------")
#     taxid2name = ncbi.get_taxid_translator([entry_taxon_id])
#     lineage = ncbi.get_lineage(entry_taxon_id)
#     print(taxid2name)
#     print(lineage)
