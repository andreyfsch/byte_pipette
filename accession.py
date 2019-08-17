from urllib.request import Request, urlopen
import xml.etree.ElementTree as Et
import gzip
import os

if not (os.path.exists("uniprot.xml") and False):
	req = Request("https://www.uniprot.org/uniprot/?query=methionine+aminopeptidase&format=xml&compress=yes&limit=10")
	req.add_header('Accept-Encoding', 'gzip')
	response = urlopen(req)
	content = gzip.decompress(response.read())
	decomp_req = content.splitlines()
	with open("uniprot.xml", "wb") as f:
		for line in decomp_req:
			f.write(line)
		f.close()
uniprot_tree = Et.parse('uniprot.xml')
uniprot_root = uniprot_tree.getroot()

entries = uniprot_root.findall("./{http://uniprot.org/uniprot}entry")

new_root = Et.Element('Life')
current_node = new_root
for entry in entries:
	for org_attr in entry.iter("{http://uniprot.org/uniprot}organism"):
		for taxon in org_attr.iter("{http://uniprot.org/uniprot}taxon"):
			child_node = current_node.find(taxon.text.replace(" ", "_").replace("/", "_"))
			if child_node is None:
				if new_root.find(".//"+current_node.tag.replace(" ", "_").replace("/", "_")) is not None:
					print(new_root.find(".//"+current_node.tag.replace(" ", "_").replace("/", "_")).tag)
				while new_root.find(".//"+current_node.tag.replace(" ", "_").replace("/", "_")) is not None:
					parent = new_root.find(".//"+current_node.tag.replace(" ", "_").replace("/", "_"))
					new_child_node = parent.find(taxon.text.replace(" ", "_").replace("/", "_"))
					if new_child_node is None:
						current_node = parent
					else:
						next_node = new_child_node
						break
				next_node = Et.SubElement(current_node, taxon.text.replace(" ", "_").replace("/", "_"))
			else:
				next_node = child_node
			current_node = next_node
	new_entry = Et.SubElement(current_node, entry.tag)
	#for child_tag in entry:
	#	new_entry.append(child_tag)
new_tree = Et.ElementTree(new_root)
new_tree.write("new_uniprot.xml")
