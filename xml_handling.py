import xml.etree.ElementTree as Et
import os
import db_handling
import accession

class UniprotXMLHandler:

	current_entry = None

	def __init__(self, xml_filename="uniprot.xml", iterate=False):
		if self.local_file_presence() and iterate:
			self.iterate_tree(xml_filename)
	
	@staticmethod
	def local_file_presence(filename="uniprot.xml"):
		return os.path.exists(filename)
	
	def iterate_tree(self, source_filename="uniprot.xml"):
		protDB = db_handling.ProteinDatabase()
		level = 0
		count = 0
		for event, element in Et.iterparse(source_filename, events=("start", "end")):
			if event == "start":
				level += 1
			elif event == "end":
				level -= 1
				if level == 1 and element.tag == "{http://uniprot.org/uniprot}entry":
					self.current_entry = element
					count += 1
					print("xml tree element "+str(count)+" added")
					protDB.insert_data(self.get_entry_data())
					element.clear()

	def get_entry_data(self):
		data = []

		data.append(('accession', self.get_entry_accession(self.current_entry)))

		data.append(('protein_existance',
			self.get_entry_protein_existance(self.current_entry)))

		data.append(('taxon_id', self.get_entry_taxon_id(self.current_entry)))

		data.append(('sequence', self.get_entry_seq(self.current_entry)))

		data.append(('isoforms', self.get_entry_isoforms(self.current_entry)))

		data.append(('refseq_ids', self.get_entry_refseq_ids(self.current_entry)))
		
		data.append(('embl_ids', self.get_entry_embl_ids(self.current_entry)))

		data.append(('ccds_ids', self.get_entry_ccds_ids(self.current_entry))) 

		data.append(('ensembl_id', self.get_entry_ensembl(self.current_entry)))

		return data
			
	@staticmethod
	def get_entry_accession(entry):
		entry_accessions = []

		for accession in entry.iter('{http://uniprot.org/uniprot}accession'):
			entry_accessions.append(accession.text)
			entry_accession = entry_accessions[0]
		
		return entry_accession
	
	@staticmethod
	def get_entry_protein_existance(entry):
		for protein_existance in entry.iter('{http://uniprot.org/uniprot}proteinExistence'):
			existance = protein_existance.attrib['type']
		
		return existance


	@staticmethod
	def get_entry_ensembl(entry):
		ensembl_id = ''
		
		for esembl_ref in entry.iter('{http://uniprot.org/uniprot}dbReference'):
			if esembl_ref.attrib['type'].startswith('Ensembl'):
				ensembl_id = esembl_ref.attrib['id']
		return ensembl_id

	@staticmethod
	def get_entry_embl_ids(entry):
		embl_ids = []
		
		for embl_ref in entry.iter('{http://uniprot.org/uniprot}dbReference'):
			if embl_ref.attrib['type'] == 'EMBL':
				value = ''
				molecule_type = ''
				status = ''
				if 'id' in embl_ref.attrib:
					id = embl_ref.attrib['id']
				else:
					id = ''
				for prop in embl_ref.iter('{http://uniprot.org/uniprot}property'):
					
					if prop.attrib['type'].find('sequence ID') != -1:
						value = prop.attrib['value']

					if prop.attrib['type'].find('molecule') != -1:
						molecule_type = prop.attrib['value']

					if prop.attrib['type'] == 'status':
						status = prop.attrib['value']

				embl_ids.append((id, value, molecule_type, status))

		return embl_ids
	
	@staticmethod
	def get_entry_refseq_ids(entry):
		refseq_ids = []
		
		for refseq_ref in entry.iter('{http://uniprot.org/uniprot}dbReference'):
			if refseq_ref.attrib['type'] == 'RefSeq':

				id = ''
				if 'id' in refseq_ref.attrib:
					id = refseq_ref.attrib['id']

				value = ''
				for prop in refseq_ref.iter('{http://uniprot.org/uniprot}property'):
					if prop.attrib['type'].find('sequence ID') != -1:
						value = prop.attrib['value']

				isoform = ''
				for molecule in refseq_ref.iter('{http://uniprot.org/uniprot}molecule'):
					if 'id' in molecule.attrib:
						isoform = molecule.attrib['id']

				refseq_ids.append((id, value, isoform))

		return refseq_ids

	@staticmethod
	def get_entry_isoforms(entry):
		isoforms = []
		for isoform in entry.iter('{http://uniprot.org/uniprot}isoform'):
			for isoform_data in isoform:
				isoform_attr = []
				if isoform_data.tag == '{http://uniprot.org/uniprot}id':
					isoform_attr.append(('id', isoform_data.text))
				if isoform_data.tag == '{http://uniprot.org/uniprot}name':
					if 'evidence' in isoform_data.attrib:
						isoform_attr.append(('name', isoform_data.text, isoform_data.attrib['evidence']))
					else:
						isoform_attr.append(('name', isoform_data.text))
				if isoform_data.tag == '{http://uniprot.org/uniprot}sequence':
					sequence = []
					for alias, prop in isoform_data.attrib.items():
						if alias == 'type':
							sequence.append(('type', prop))
						if alias == 'ref':
							sequence.append(('ref', prop))

					isoform_attr.append(('sequence', sequence))

				isoforms.append(isoform_attr)

		return isoforms
	
	@staticmethod
	def get_entry_ccds_ids(entry):
		ccds_ids = []
		
		for ccds_ref in entry.iter('{http://uniprot.org/uniprot}dbReference'):
			if ccds_ref.attrib['type'] == 'CCDS':
				molecule_ref = ''
				for molecule in ccds_ref:
					if 'id' in molecule.attrib:
						molecule_ref = molecule.attrib['id']
				ccds_ids.append((ccds_ref.attrib['id'], molecule_ref))


		return ccds_ids
	
	@staticmethod
	def get_entry_taxon_id(entry):
		
		for taxon_ref in entry.iter('{http://uniprot.org/uniprot}dbReference'):
			if taxon_ref.attrib['type'] == 'NCBI Taxonomy':
				taxon_id = taxon_ref.attrib['id']

		return taxon_id

	@staticmethod
	def get_ccds_seq(html_content):
		ccds_seq = ''

		soup = BeautifulSoup(html_content, 'html.parser')

		tt = soup.find("tt") 

		spans = tt.find_all("span")

		for span in spans:
			font = span.find("font")
			ccds_seq += font.get_text()

		return ccds_seq


	@staticmethod
	def get_entry_seq(entry):
		entry_seq = ''
		version = 0
		for sequence_element in entry.iter('{http://uniprot.org/uniprot}sequence'):
			if 'version' in sequence_element.attrib:
				if int(sequence_element.attrib['version']) > version:
					entry_seq = sequence_element.text
			else:
				entry_seq = sequence_element.text
		
		return entry_seq

	@staticmethod
	def get_entry_seq_caution(entry):
		pass
	
	@staticmethod
	def get_entry_evidences(entry):
		pass