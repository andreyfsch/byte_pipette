from urllib.request import Request, urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup
import xml.etree.ElementTree as Et
import gzip
from xml_handling import UniprotXMLHandler
import db_handling
import re


class Accession:
	
	@staticmethod
	def download_uniprot_xml(limit=False, force=False):
		if (not UniprotXMLHandler.local_file_presence("uniprot.xml")) or force:
			url = 'https://www.uniprot.org/uniprot/?query=family%3A%22peptidase+m24a+family%22&format=xml&compress=yes'
			if limit:
				url += '&limit='+str(limit)
			print('downloading from '+url)
			req = Request(url)
			req.add_header('Accept-Encoding', 'gzip')
			response = urlopen(req)
			content = gzip.decompress(response.read())
			decomp_req = content.splitlines()
			with open("uniprot.xml", "wb") as f:
				for line in decomp_req:
					f.write(line+'\n'.encode())
				f.close()

	@staticmethod
	def get_embl_seq(accession_id):
		url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"+str(accession_id)
		req = Request(url)
		try:
			response = urlopen(req)
		except HTTPError as err:
			return err
		fasta = response.read().decode('utf-8')
		lines = fasta.split('\n')
		lines.pop(0)
		return ''.join(lines)

	@staticmethod
	def get_ddbj_seq(accession_id):
		url = "http://getentry.ddbj.nig.ac.jp/getentry/dad/"+str(accession_id)
		req = Request(url)
		try:
			response = urlopen(req)
		except HTTPError as err:
			return err
		fasta = response.read().decode('utf-8')
		lines = fasta.split('\n')
		lines.pop(0)
		return ''.join(lines)

	def get_refseq_value_from_id(self, accession_id, protein=False):
		params = self.get_html_refseq(accession_id,protein=protein)

		if type(params[0]) is not HTTPError:
			seq_value = self.get_refseq_seq(params[0],params[1],protein=protein)
		else:
			seq_value = params[0]

		return seq_value

	
	@staticmethod
	def get_html_refseq(accession_id, protein=False):
		if not protein:
			url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id="+accession_id+"&db=nuccore&report=genbank&conwithfeat=on&withparts=on&hide-cdd=on&retmode=html"
		else:
			url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id="+accession_id+"&db=protein&report=genbank&conwithfeat=on&withparts=on&hide-cdd=on&retmode=html"
		print(url)
		req = Request(url)
		try:
			response = urlopen(req)
		except HTTPError as err:
			return err, err
		
		if not protein:
			req_main = Request("https://www.ncbi.nlm.nih.gov/nuccore/"+accession_id)
			print("https://www.ncbi.nlm.nih.gov/nuccore/"+accession_id)
		else:
			req_main = Request("https://www.ncbi.nlm.nih.gov/protein/"+accession_id)
			print("https://www.ncbi.nlm.nih.gov/protein/"+accession_id)
		try:
			response_main = urlopen(req_main)
		except HTTPError as err:
			return err, err

		return response.read(), response_main.read()
	
	@staticmethod
	def get_refseq_html(accession_id, params):
		url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id="+str(accession_id)+"&"+params
		req = Request(url)
		try:
			response = urlopen(req)
		except HTTPError as err:
			return err
		fasta = response.read().decode('utf-8')
		lines = fasta.split('\n')
		lines.pop(0)
		return ''.join(lines)

	@staticmethod
	def get_html_ccds(accession_id):
		url = "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA="+str(accession_id)
		req = Request(url)
		response = urlopen(req)
		return response.read()


	def get_transl_table(self, accession_id):
		html, main_html = self.get_html_refseq(accession_id,protein=True)

		soup = BeautifulSoup(html, 'html.parser')

		a_links = soup.find_all(href=re.compile('wprintgc.cgi'))
		
		transl_table = 1

		for a_link in a_links:
			if 'sfeat' in a_link:
				print(a_link['sfeat'])
			if a_link['href'].find('Taxonomy/Utils') != -1:
				transl_table = int(a_link.get_text())

		return transl_table


	def get_refseq_seq(self, html_content_features, html_content_main, protein=False):
		soup_features = BeautifulSoup(html_content_features, 'html.parser')
		soup_main = BeautifulSoup(html_content_main, 'html.parser')

		if not protein:
			warning = soup_main.findAll('li', {'class':'hi_warn icon'})
			if warning:
				for li in warning:
					span = li.find("span", {'class': 'icon'})
					if span and span.get_text().startswith('Record removed.'):
						detail = span.find('em').get_text()
						return 'sequence_removed'
			
			div_description = soup_main.find("div", {'id', 'rprtheader'})
			h1_description = div_description.find("h1")
			description = h1_description.get_text()
			if (description.find('cds') == -1 and description.find('complete') != -1) or description.find('WGS') != -1 or description.find('chromosome') != -1 or description.find('whole') != -1:
				return 'not_cds_seq'
			elif description.find('partial') != -1:
				return 'partial'
			else:
				a_links = soup_features.findAll("a")

				for a_link in a_links:
					if a_link.get_text().find('CDS') != -1:
						link = a_link
			
				href = link['href']

				href = link['href']		
				params = href.split('?')[-1]
				div = soup_main.find("div", {'id':'viewercontent1'})

				reference_id = div['val']

				refseq_seq = self.get_refseq_html(reference_id, params)
		else:
			a_links = soup_features.findAll("a")

			for a_link in a_links:
				if a_link.get_text().find('CDS') != -1:
					link = a_link

		
			href = link['href']		
			reference_id = href.split('/')[-1].split('?')[0]
			params = href.split('/')[-1].split('?')[-1]
			refseq_seq = self.get_refseq_html(reference_id, params)

		return refseq_seq

	def handle_404_error(self, reference_id):
		db_handler = db_handling.ProteinDatabase()
		embl_again = self.get_embl_seq(reference_id)
		if type(embl_again) is not HTTPError:
			print(reference_id, 'EMBL', embl_again)
			db_handler.update_biological_sequence({'value': embl_again, 'error_404':False}, reference_id)
		else:
			refseq_trial = self.get_refseq_value_from_id(reference_id, protein=True)
			if type(refseq_trial) is not HTTPError:
				print(reference_id, 'RefSeq', refseq_trial)
				db_handler.update_biological_sequence({'value': refseq_trial,'error_404':False}, reference_id)
			else:
				ddbj_trial = self.get_ddbj_seq(reference_id)
				if type(ddbj_trial) is not HTTPError:
					print(reference_id, 'DDBJ', ddbj_trial)
					db_handler.update_biological_sequence({'value': ddbj_trial,'error_404':False}, reference_id)
		