import db_handling
import biological_sequence
from typing import Any


class ProteinEntry:

    def __init__(self, protein_id:int = None) -> None:
        
        if protein_id:
            prot_db = db_handling.ProteinDatabase()

            protein_entry = prot_db.get_protein_entry(protein_id)

            if protein_entry:

                for attr_name, attr_value in protein_entry.items():
                    if attr_name in ProteinEntry.__dict__:
                        setattr(self, attr_name, attr_value)

                self.cds_seq = prot_db.get_ellected_seq(self.id)
                self.aa_seq = prot_db.get_uniprot_seq(self.id)


    @property
    def id(self) -> int:
        return self._id


    @id.setter
    def id(self, id:int) -> None:
        self._id = id


    @property
    def protein_existance(self) -> str:
        return self._protein_existance


    @protein_existance.setter
    def protein_existance(self, existance:str) -> None:
        self._protein_existance = existance


    @property
    def ncbi_taxonomy(self) -> int:
        return self._ncbi_taxonomy


    @ncbi_taxonomy.setter
    def ncbi_taxonomy(self, taxonomy_id) -> None:
        self._ncbi_taxonomy = taxonomy_id


    @property
    def transl_table(self) -> int:
        return self._transl_table


    @transl_table.setter
    def transl_table(self, table_id:int) -> None:
        self._transl_table = table_id


    @property
    def represented_by(self) -> int:
        return self._represented_by


    @represented_by.setter
    def represented_by(self, id:int) -> None:
        self._represented_by = id


    @property
    def representative_of_taxon(self) -> int:
        return self._representative_of_taxon


    @representative_of_taxon.setter
    def representative_of_taxon(self, id:int) -> None:
        self._representative_of_taxon = id


    @property
    def cds_seq(self) -> biological_sequence.BiologicalSequence:
        return self._cds_seq


    @cds_seq.setter
    def cds_seq(self, seq_id:int) -> None:
        self._cds_seq = biological_sequence.BiologicalSequence(seq_id)


    @property
    def aa_seq(self) -> biological_sequence.BiologicalSequence:
        return self._aa_seq


    @aa_seq.setter
    def aa_seq(self, seq_id:int) -> None:
        self._aa_seq = biological_sequence.BiologicalSequence(seq_id)