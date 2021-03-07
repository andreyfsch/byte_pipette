import mysql.connector
import alignment
from ete3 import NCBITaxa
import sys

class ProteinDatabase:
    host = "localhost"
    user="root"
    password="root"
    database="protein"
    conn = None
    cursor = None

    def __init__(self):
        self.conn = mysql.connector.connect(host=self.host,
         user=self.user, password=self.password, database=self.database)

        self.cursor = self.conn.cursor(buffered=True)
    
    def renew_conn(self):
        self.cursor.close()
        self.conn.close()
        self.__init__()

    def  get_protein_entry(self, id):
        sql = "select * from protein_entry where id = %s"

        self.cursor.execute(sql,[id])

        seqs = self.cursor.fetchone()

        return seqs

    

    def get_rank_list(self, rank_level):
        previous_rank = alignment.Alignment.previous_taxon_rank(rank_level)

        sql = "select distinct "+rank_level+"_rank_id from protein_entry\
        where representative_of_taxon is not null\
        and representative_taxon_rank = %s\
        and "+rank_level+"_rank_id is not null"

        self.cursor.execute(sql, [previous_rank])

        # print(self.cursor.statement)

        rank_ids = self.cursor.fetchall()

        return rank_ids

    def get_representative_entries(self):
        sql = "select pe.id, pe.ncbi_taxonomy\
        from protein_entry pe\
        inner join biological_sequence bs on(bs.protein_id = pe.id)\
        where bs.ellected_seq = true\
        and (pe.class_rank_id is null or  pe.family_rank_id is null or pe.genus_rank_id is null\
            or pe.phylum_rank_id is null or pe.species_rank_id is null\
            or pe.order_rank_id is null)\
        and bs.db != 'UniProt'"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def get_entries_no_representative(self):
        sql = "select pe.id, pe.representative_of_taxon, pe.taxon_name_representative, pe.ncbi_taxonomy \
        from protein_entry pe where pe.representative_of_taxon is not null \
        and pe.represented_by is null"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def update_protein_entry_by_repres_by(self, values, repres_by_id):
        sql = "update biological_sequence set "
        count = 0
        update_vals = []
        for col, val in values.items():
            count += 1
            update_vals.append(val)
            sql+= col+"=%s"
            if count < len(values):
                sql+= ","
        sql+=" where represented_by=%s"
        update_vals.append(id)

        self.cursor.execute(sql,update_vals)

        # print(self.cursor.statement)

        self.conn.commit()

    def get_erroneous_repres_of_taxon_instances(self):
        sql = "select pe.id, pe.representative_of_taxon, pe.represented_by \
        from protein_entry pe where pe.representative_of_taxon is not null \
        and pe.represented_by is not null"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def get_count_of_children_of_repres(self, repres_id):
        sql = "select count(*) from protein_entry pe where pe.represented_by = %s"

        self.cursor.execute(sql,[repres_id])

        # print(self.cursor.statement)

        seqs = self.cursor.fetchone()

        return seqs[0]

    def get_comparisons_same_taxon(self, leaf=True):
        mid_sql = ",pe.id, pe2.id, na.identical_sequences,\
		na.identity_90, na.identity_85, na.identity_80 from biological_sequence bs \
        inner join protein_entry pe on(pe.id = bs.protein_id)\
        inner join needle_alignment na on(bs.id = na.id_a)\
        inner join biological_sequence bs2 on(bs2.id = na.id_b)\
        inner join protein_entry pe2 on(bs2.protein_id = pe2.id)\
        where bs.ellected_seq = 1 and bs2.ellected_seq = 1 \
        and na.identity_80 = 1 and (pe.represented_by is null or pe2.represented_by is null)"

        if leaf:
            sql = "select pe.ncbi_taxonomy"
            sql += mid_sql
            sql += " and pe.representative_of_taxon is null \
            order by pe.ncbi_taxonomy"
        else:
            sql = "select pe.representative_of_taxon"
            sql += mid_sql
            sql += " and pe.representative_of_taxon is not null \
            order by pe.representative_of_taxon"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs
    
    def get_entry_ids_representative_of_taxon(self, taxon_id):
        sql = "select pe.id, pe.representative_of_taxon from biological_sequence bs \
        inner join protein_entry pe on(bs.protein_id = pe.id) where \
        pe.representative_of_taxon = %s"

        self.cursor.execute(sql,[taxon_id])

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def get_entry_ids_from_taxon(self, taxon_id, leaf=True):
        sql = "select pe.id, pe.ncbi_taxonomy from biological_sequence bs \
        inner join protein_entry pe on(bs.protein_id = pe.id) where"
        
        if leaf:
            sql += " pe.ncbi_taxonomy = %s and"
        else:
            sql += " pe.ncbi_taxonomy in ("
        
            for i in range(len(taxon_id)):
                sql += "%s"
                if i != len(taxon_id) - 1:
                    sql += ","
            
            sql += ") and"

        sql += " bs.ellected_seq = 1"

        if leaf:
            self.cursor.execute(sql,[taxon_id])
        else:
            self.cursor.execute(sql, taxon_id)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def get_sequences_same_taxon(self, representative=False):
        sql = "select bs.id, bs.value, pe.ncbi_taxonomy from biological_sequence bs \
        inner join protein_entry pe on(pe.id = bs.protein_id) \
        where bs.ellected_seq = 1"

        if representative:
            sql += " and pe.representative_of_taxon is not null"

        sql += " order by(pe.ncbi_taxonomy)"

        self.cursor.execute(sql)

        seqs = self.cursor.fetchall()

        return seqs

    def get_sequences_same_taxon_rank(self, taxon_rank):


        sql = "select bs.id, bs.value, pe."+taxon_rank+"_rank_id "

        if taxon_rank != 'species':
            previous_rank = self.get_previous_state_taxon_rank()
            sql += ", pe."+previous_rank+"_rank_id "

        sql += "from biological_sequence bs inner join protein_entry pe on(pe.id = bs.protein_id)\
        where bs.ellected_seq = 1 and pe.represented_by is null and pe."+taxon_rank+"_rank_id is not null "

        if taxon_rank != 'species':
            sql += "and (pe.representative_taxon_rank = '"+previous_rank+"' or pe.representative_taxon_rank is null)"

        sql += " order by pe."+taxon_rank+"_rank_id, bs.id"

        self.cursor.execute(sql)

        seqs = self.cursor.fetchall()

        return seqs
    
    def get_sequences_same_repres(self):
        sql = "select bs.id, bs.value, pe.representative_of_taxon from biological_sequence bs \
        inner join protein_entry pe on(pe.id = bs.protein_id) \
        where bs.ellected_seq = 1 and pe.representative_of_taxon is not null \
        and pe.represented_by is null \
        order by(pe.representative_of_taxon)"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs

    def update_biological_sequence(self, values, id):
        sql = "update biological_sequence set "
        count = 0
        update_vals = []
        for col, val in values.items():
            count += 1
            update_vals.append(val)
            sql+= col+"=%s"
            if count < len(values):
                sql+= ","
        sql+=" where id=%s"
        update_vals.append(id)

        self.cursor.execute(sql,update_vals)

        # print(self.cursor.statement)

        self.conn.commit()

    def update_protein_entry(self, values, id):
        sql = "update protein_entry set "
        count = 0
        update_vals = []
        for col, val in values.items():
            count += 1
            update_vals.append(val)
            sql+= col+"=%s"
            if count < len(values):
                sql+= ","
        sql+=" where id=%s"
        update_vals.append(id)

        self.cursor.execute(sql,update_vals)

        # print(self.cursor.statement)

        self.conn.commit()

    def insert_needle_alignment(self, values):
        sql = "insert into needle_alignment (id_a, id_b,\
        align_length, align_identity, align_similarity, \
        align_gaps, algin_score) values (%s,%s,%s,%s,%s,%s,%s)"

        self.cursor.execute(sql,values)

        # print(self.cursor.statement)

        self.conn.commit()

        self.renew_conn()
    
    def verify_comparison(self, seq_a, seq_b):
        possibility_1 = None
        possibility_2 = None
        sql = "select * from needle_alignment where id_a = %s and id_b = %s"

        self.cursor.execute(sql,[seq_a,seq_b])

        for whatever in self.cursor:
            possibility_1 = whatever
            # print(self.cursor.statement)
        
        if possibility_1:
            return False

        self.renew_conn()

        self.cursor.execute(sql,[seq_b,seq_a])

        for whatever in self.cursor:
            possibility_2 = whatever
            # print(self.cursor.statement)

        return (not possibility_1) and (not possibility_2)

    def get_comparison(self, seq_a, seq_b, treshold=80):
        sql = "select identity_"+str(treshold)+" from needle_alignment where id_a = %s and id_b = %s"

        self.cursor.execute(sql,[seq_a,seq_b])

        for whatever in self.cursor:
            possibility_1 = whatever
            # print(self.cursor.statement)
        
        if possibility_1:
            return self.cursor.fetchone()

        self.renew_conn()

        self.cursor.execute(sql,[seq_b,seq_a])

        for whatever in self.cursor:
            possibility_2 = whatever
            # print(self.cursor.statement)

        if possibility_2:
            return self.cursor.fetchone()
        else:
            return False

    def get_taxon_ids_by_rank(self, taxon_rank):
        sql = "select distinct "+taxon_rank+"_rank_id from protein_entry\
            where "+taxon_rank+"_rank_id is not null order by "+taxon_rank+"_rank_id"

        self.cursor.execute(sql)

        # print(self.cursor.statement)

        return self.cursor.fetchall()

    def update_representatives_taxon_id_by_abstence(self, taxon_rank_id, taxon_rank):
        ncbi = NCBITaxa()
        
        taxon_name = ncbi.get_taxid_translator([taxon_rank_id])[taxon_rank_id]

        sql = "update protein_entry set representative_of_taxon = %s, taxon_name_representative = %s\
            ,representative_taxon_rank = %s\
            where represented_by is null and "+taxon_rank+"_rank_id = %s "

        self.cursor.execute(sql, [taxon_rank_id, taxon_name, taxon_rank, taxon_rank_id])

        # print(self.cursor.statement)

        # if self.cursor.rowcount != 0:
        #     print("updated {} remaining sequences left represented by noone".format(self.cursor.rowcount))
        # else:
        #     print('already done')

        self.conn.commit()

    def get_comparisons_same_taxon_id(self, taxon_rank_id, taxon_rank='species', treshold=80):
        sql = "select pe.id, pe2.id\
        from biological_sequence bs\
            inner join protein_entry pe on(pe.id = bs.protein_id)\
            inner join needle_alignment na on(bs.id = na.id_a)\
            inner join biological_sequence bs2 on(bs2.id = na.id_b)\
            inner join protein_entry pe2 on(bs2.protein_id = pe2.id)\
        where bs.ellected_seq = 1 and bs2.ellected_seq = 1\
            and na.identity_80 = 1 and pe.represented_by is null and pe2.represented_by is null\
            and pe."+taxon_rank+"_rank_id = %s and  pe2."+taxon_rank+"_rank_id = %s\
        order by pe.id"

        self.cursor.execute(sql, [taxon_rank_id,taxon_rank_id])

        # print(self.cursor.statement)

        return self.cursor.fetchall()


    def get_previous_state_taxon_rank(self):
        current_seq = self.get_current_state()[6]

        sql = "select taxon_rank from filtering_process_log where seq = %s"

        self.cursor.execute(sql, [current_seq-1])

        self.conn.commit()

        return self.cursor.fetchone()[0]


    def get_current_state(self):
        sql = "select * from filtering_process_log where current_task = 1"
        
        self.cursor.execute(sql)

        self.conn.commit()

        return self.cursor.fetchone()

    def get_filtering_states(self):
        sql = "select * from filtering_process_log order by seq"
        
        self.cursor.execute(sql)

        self.conn.commit()

        return self.cursor.fetchall()

    def update_state(self, aspects, rank):
        sql = "update filtering_process_log set "

        count = 0
        update_vals = []
        for col, val in aspects.items():
            count += 1
            update_vals.append(val)
            sql+= col+"=%s"
            if count < len(aspects):
                sql+= ","
        sql+=" where taxon_rank=%s"
        update_vals.append(rank)

        self.cursor.execute(sql,update_vals)

        self.conn.commit()

    def set_next_filtering_state(self):

        sql = "update filtering_process_log set current_task = 1 where seq = %s"

        current_seq = self.get_current_state()[6]

        self.cursor.execute(sql,[current_seq+1])

        # print(self.cursor.statement)

        self.conn.commit()
        

    def assign_identity_paramters_comparisons(self):

        sql = "delete from needle_alignment where id_a = id_b"

        self.cursor.execute(sql)

        self.conn.commit()

        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identical_sequences = 0\
        where na.align_identity != na.align_length or \
        na.align_length != na.align_similarity or na.align_gaps != 0\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identical_sequences = 1\
        where na.align_identity = na.align_length and \
        na.align_length = na.align_similarity and na.align_gaps = 0\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_90 = 1 \
        where round(na.align_identity/na.align_length, 2) >= 0.90\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_90 = 0\
        where round(na.align_identity/na.align_length, 2) < 0.90\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_85 = 1 \
        where round(na.align_identity/na.align_length, 2) >= 0.85\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_85 = 0\
        where round(na.align_identity/na.align_length, 2) < 0.85\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_80 = 1 \
        where round(na.align_identity/na.align_length, 2) >= 0.80\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()
        
        self.renew_conn()

        sql = "update needle_alignment na\
        inner join biological_sequence bs on na.id_a = bs.id\
        inner join biological_sequence bs2 on na.id_b = bs2.id\
        set na.identity_80 = 0\
        where round(na.align_identity/na.align_length, 2) < 0.80\
        and bs.ellected_seq = 1 and bs2.ellected_seq = 1"

        self.cursor.execute(sql)

        self.conn.commit()

    def get_rank_id_sequences(self, taxon_rank, taxon_rank_id):
        sql = "select pe.id\
        from protein_entry pe\
            inner join biological_sequence bs on (bs.protein_id = pe.id)\
        where pe."+taxon_rank+"_rank_id = %s\
            and bs.ellected_seq = 1\
        order by pe.id"

        self.cursor.execute(sql,[taxon_rank_id])

        # print(self.cursor.statement)

        seqs = self.cursor.fetchall()

        return seqs
