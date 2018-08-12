#!/usr/bin/env python3.6

## MODULE FOR PARSING XML TO TAB


## the functions below extract info from the "description" lines of the XML file...called from the BioPython parser by tags like "hsp.hit_description"

def get_gene(x):
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "gene=" in str(i):
			descrip = (str(i)[5:])
	descrip = str("none")
	return descrip

def get_location(x):
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "location=" in str(i):
			descrip = (str(i)[9:])
			return descrip
	descrip = str("none")
	return descrip

def get_GC(x):
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "GCpct=" in str(i):
			descrip = (str(i)[6:])
			return descrip
	descrip = str("none")
	return descrip

def get_locusid(x): #this will give you the full locus id of the object from the description line [locus_id=XXXX####]
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "locus_tag=" in str(i):
			locus_id = (str(i)[10:])
			return locus_id
	locus_id = "no_locus_id"
	return locus_id

def get_locus_tag(x): #this will give you the locus tag at the front of the locusid.  There are two general formats of locus tags.
	import re
	x = str(x.split("|")[0]) #sometimes we pass items like "Y75_RS14370|WP_000220066.1|0.84" and we only want the front part.
	if "_" in x:
	    locus_tag = str(x.split("_")[0]) #this works for locus tags formatted like STM14_#######
	    return locus_tag
	else:
	    locus = re.match("\D+", x) #this works for locus tags formatted like STM#### (no hyphen...letters before a string of #s)
	    locus_tag = locus.group(0)
	    return locus_tag

def get_protein(x):
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "protein=" in str(i):
			descrip = (str(i)[8:])
			return descrip
		if "product=" in str(i):
			descrip = (str(i)[8:])
			return descrip
		if "gbkey=tRNA" in str(i):
			return "tRNA"
	descrip = str("none")
	return descrip

def get_proteinid(x):
	descriptors = str(x[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for i in descriptor_list:
		if "protein_id=" in str(i):
			descrip = (str(i)[11:])
			return descrip
		if "pseudo=true" in str(i):
			descrip = "PSUEDOGENE"
			return descrip
		if "product=" in str(i):
			descrip = (str(i)[8:])
			return descrip
		if "gbkey=tRNA" in str(i):
			return "tRNA"
	descrip = str("none")
	return descrip


def pxml(infile, outfile, strain_dictionary, pctid_threshold):
	import numpy
	import re
	from Bio import SearchIO
	from Bio import SeqIO
	from Bio.SeqUtils import GC
	
	output_tab_file = open(outfile, "w")
	
	##Below we set the percent identity required for inclusion in our list.  The default value is set to 70%.
	try:
		pctid_threshold = float(pctid_threshold) / float('100')
	except TypeError:
		pctid_threshold = float("0.7")  #set the default to 70%


	### Here we read the "dictionary" file and make our list.  
	### The format should be: Locus tag | BiosampleID | species and strain

	strain_list = []
	strain_locus_tag_list = []

	with open(strain_dictionary, "r") as f:
		for line in f:
			strain_list.append(line.rstrip())
			dict_locus_tag = line.split(" | ")[0]
			strain_locus_tag_list.append(dict_locus_tag)
		
	strains = "\t".join(strain_list)
	### HERE WE IMPORT DATA FROM THE FILE WITH THE NAME GIVEN FROM THE COMMAND LINE
	qresults = SearchIO.parse(infile, 'blast-xml')

	output_tab_file.write("LOCUS ID\tACCESSION\tGENE NAME\tDNA ACCESSION\tLOCATION\tPROTEIN ID\tGENE SIZE\tPCT_GC\tNUM OF HITS\tSTRAINS WITH HITS\tAVG % ID\tDESCRIPTION\t" + strains +"\n")

	for qresult in qresults:
		query_accessions = str(qresult.id)
		query_DNA_accession = str(query_accessions.split("_cds_")[0])[4:] #retrieves only the nucleotide portion
		query_description_line = ()
		query_length = int()
		query_length = qresult.seq_len
		query_description_line = qresult.description
		query_locus_id = get_locusid(qresult.description)
		query_protein_id = str(get_proteinid(qresult.description))
		query_location = str(get_location(qresult.description))
		query_gene_name = str(get_gene(qresult.description))
		query_GC = str(get_GC(qresult.description))
		percent_id_list = []
		locuslist = set()
		descriptions = set()
		hit = str()
		hsp = str()

		for hit in qresult:
			for hsp in hit:
				param = float(qresult.seq_len)/float(hit.seq_len)
				pctid = float(hsp.ident_num)/float(hit.seq_len)
				coverage = str(float(hsp.aln_span)/float(qresult.seq_len))
				percent_id = str(pctid)[0:4]
				coverage = coverage[0:4]
				if param < 1.2 and param > 0.8 and float(pctid) > pctid_threshold and float(coverage) >0.8:
					hit_accession = str(hsp.hit_id)
					hit_DNA_accession = str(hit_accession.split("_cds_")[0])
					percent_id_list.append(pctid)
					locus_id = get_locusid(hsp.hit_description)
					protein = get_protein(hsp.hit_description)
					location = get_location(hsp.hit_description)
					proteinid = get_proteinid(hsp.hit_description)				#hitinfo = str(hit_DNA_accession) + " (" + str(locus_id) + ") |" + str(proteinid) + "|" + percent_id
					hitinfo = str(hit_accession) + "|" + str(locus_id) + "|" + str(proteinid) + "|" + percent_id
					locuslist.add(hitinfo)
					descriptions.add(protein)
		if len(percent_id_list) >= 1:
			average_percent_id = str(numpy.mean(percent_id_list))[0:5]
		else:
			average_percent_id = "N/A"

		if hit == "":
			descriptions.add("no hits that exceeded evalue")
		else:
			descriptions.add(get_protein(hsp.fragment.query.description))
	
		descriptionlist = list(descriptions)
		descriptionlist.sort()
		list_of_descriptions = ", ".join(descriptionlist)

		hit_locus_list = list(locuslist)
		hit_locus_list.sort()

		strain_hits = []
		for i in strain_locus_tag_list:
			matchlist = []
			for j in hit_locus_list:
				j_locus = str(j.split("|")[1])
				j_locus_tag = get_locus_tag(j_locus)
				if j_locus_tag == i:
					matchlist.append(j)
			if len(matchlist) >= 1 : matches = "["+str(len(matchlist))+"]|"+str(", ".join(matchlist) + "\t")
			else: matches = str("*\t")
			strain_hits.append(matches)
		number_of_strains_with_hits = len(strain_hits) - strain_hits.count("*\t")
		list_of_hit_loci = "".join(strain_hits)
		output_tab_file.write(str(query_locus_id) + "\t" + str(query_accessions)[4:] + "\t" + query_gene_name + "\t" + query_DNA_accession + "\t" + query_location + "\t" + query_protein_id + "\t" + str(qresult.seq_len) + "\t" + query_GC + "\t" + str(len(hit_locus_list)) + "\t" + str(number_of_strains_with_hits) + "\t" + average_percent_id + "\t" + list_of_descriptions + "\t" + list_of_hit_loci + "\n")
	return
