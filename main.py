from accession import Accession
from xml_handling import UniprotXMLHandler
from ete3 import NCBITaxa
import db_handling
import alignment


#Accession.download_uniprot_xml()

# uniprotHandler = UniprotXMLHandler(iterate=True)

# ncbi = NCBITaxa()

align = alignment.Alignment()

protDB = db_handling.ProteinDatabase()

# acc = Accession()

# protDB.set_sequences_values()

# protDB.fix_sequences_multiple_repres()

repres_seqs = protDB.get_sequences_same_taxon()
align.compare_sequences_same_entry(repres_seqs)


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
