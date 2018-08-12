#!/usr/bin/env python3.6

def convert_gbk_to_fna(filename):
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	from Bio.SeqUtils import GC
	input_handle = open(filename, "r")
	output_handle = open(filename.split(".")[0] + ".fna", "w")

	index=0
	rna_index=0
	seq_accession_index = ""
	for seq_record in SeqIO.parse(input_handle, "genbank") :
		title = (seq_record.annotations["accessions"][0])
		print("Dealing with GenBank record %s" % seq_record.id)
		seq_accession_prefix = seq_record.id.split("_")[0]
		if seq_accession_prefix == "NC":
			index = 0
			rnaindex = 0
		if seq_accession_prefix == "NZ":
			new_seq_accession_index = seq_record.id[:7]
			if new_seq_accession_index != seq_accession_index: #must restart index numbering for any new cluster of sequences
				index = 0
				rnaindex = 0
		locus_number = 5 #setting a base locus number in case the gbk file does not automatically have one for its CDS features.
		for seq_feature in seq_record.features :
			if seq_feature.type=="CDS":
				index += 1
				seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
				GC_content = str(GC(seq_feature.extract(seq_record).seq))
				#location = location_format(str(seq_feature.location))
				location = str(seq_feature.location)
				if "locus_tag" in seq_feature.qualifiers:
					locus_tag = seq_feature.qualifiers["locus_tag"][0]    
				if "locus_tag" not in seq_feature.qualifiers:
					title = (seq_record.annotations["accessions"][0])
					locus_tag = (title + "_" + str(locus_number).zfill(5))
					seq_feature.qualifiers["locus_tag"] = [locus_tag]
					locus_number += 5
				if "protein_id" in seq_feature.qualifiers:
					protein_id = seq_feature.qualifiers["protein_id"][0]
					protein_accession = "_"+protein_id
				if "translation" in seq_feature.qualifiers:
					translation = (seq_feature.qualifiers["translation"][0])
				if "product" in seq_feature.qualifiers:
					product = seq_feature.qualifiers["product"][0]
				if "pseudo" in seq_feature.qualifiers:
					X = (seq_feature.extract(seq_record).seq)
					translation = (X.translate()[:-1])
					protein_id = "PSEUDOGENE"
					protein_accession=""
				output_handle.write(">lcl|%s_cds%s%s [locus_tag=%s] [protein=%s] [protein_id=%s] [location=%s] [GCpct=%s] [gbkey=CDS]\n%s\n" % (
					seq_record.id,
					protein_accession,
					str("_"+str(index)),
					locus_tag,
					product,
					protein_id,
					location,
					'{:.4}'.format(GC_content),
					seq_feature_sequence))
					
			if seq_feature.type=="tRNA" or seq_feature.type=="rRNA" or seq_feature.type=="ncRNA" or seq_feature.type=="tmRNA":
				rna_index += 1
				seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
				if "locus_tag" in seq_feature.qualifiers:
					locus_tag = seq_feature.qualifiers["locus_tag"][0]
				if "product" in seq_feature.qualifiers:
					product = seq_feature.qualifiers["product"][0]                    
				#location = location_format(str(seq_feature.location))
				location = str(seq_feature.location)
				GC_content = str(GC(seq_feature.extract(seq_record).seq))
				protein_id = product
				protein_accession = ""
				output_handle.write(">lcl|%s_%s%s [locus_tag=%s] [product=%s] [location=%s] [GCpct=%s] [gbkey=%s]\n%s\n" % (
					seq_record.id,
					seq_feature.type.lower(),
					str("_"+str(rna_index)),
					locus_tag,
					product,
					location,
					'{:.4}'.format(GC_content),
					seq_feature.type,
					seq_feature_sequence))

	output_handle.close()
	return

def convert_fna_to_faa(filename):
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	output_handle = open("./faa/" + filename.split(".")[0] + ".faa", "w")
	for seq_record in SeqIO.parse(filename, "fasta"):
		if get_gbkey(seq_record.description) == "CDS":
			output_handle.write(">" + seq_record.description[4:])
			output_handle.write(line_format(str(seq_record.seq.translate()[:-1]))+"\n") #the [:-1] removes the * for the stop codon from the sequence record
	output_handle.close()
	print("Done")

def get_gbkey(description_line):
	descriptors = str(description_line[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for item in descriptor_list:
		if "gbkey=" in str(item):
			print(item[6:])
			descrip = (str(item)[6:])
			return descrip
	descrip = str("none")
	return descrip

def get_locusid(seq_record): 
	for seq_feature in seq_record.features:
		if seq_feature.type=="CDS":
			if "locus_tag" in seq_feature.qualifiers:
				locus_id = seq_feature.qualifiers["locus_tag"][0]
				return locus_id

def get_locus_tag(x): #this will give you the locus tag at the front of the locusid.  There are two general formats of locus tags.
	import re
	if x is None:
		return
	if "_" in x:
	    locus_tag = str(x.split("_")[0]) #this works for locus tags formatted like STM14_#######
	    return locus_tag
	else:
	    locus = re.match("\D+", x) #this works for locus tags formatted like STM#### (no hyphen...letters before a string of #s)
	    locus_tag = locus.group(0)
	    return locus_tag

def line_format(sequence): #this function returns a FASTA seqeunce with 80 characters per line
        newline_list = []
        while sequence:
                newline = sequence[0:80]
                newline_list.append(newline)
                sequence = sequence[80:]
        output = "\n".join(newline_list)
        return output

def dictmaker(filename):
	from Bio import SeqIO
	input_handle = open(filename, "r")
	locus_tag_set = set()
	locus_descrip_list = []
	for seq_record in SeqIO.parse(input_handle, "genbank"):
		description = seq_record.description
		accession = (" | ".join(seq_record.dbxrefs))
		locusid = get_locusid(seq_record)
		locustag = get_locus_tag(locusid)
		if locustag is not None:
			if locustag not in locus_tag_set:
				locus_tag_set.add(locustag)
				dict_entry = str(locustag + " | " + accession + " | " + description)
				locus_descrip_list.append(dict_entry)
	return dict_entry
