#!/usr/bin/env python3.6
import os
import nav
import parser

project_name = input("What is the name of the dict and blast db for this project? ")
pct_id = input("What pct ID cutoff for the tab files? ")

genbank_file_list =[]
directory_list = []
for entry in os.scandir('.'):
    if entry.is_file():        	
        if str(entry.name)[0] != "." and str(entry.name.split(".")[-1]) == "gbff":
            with open(entry) as fp:
                line = str(fp.readline())
                if "LOCUS" in line:
                    #print(str(i)+". "+entry.name + "\t{}".format(line.strip()))
                    genbank_file_list.append(entry.name)
    if entry.is_dir():
    	directory_list.append(entry.name)
    	
#if not already in the current directory make the directories gbk and fna
if "gbk" not in directory_list:
	os.mkdir("gbk")
if "fna" not in directory_list:
	os.mkdir("fna")
if "xml" not in directory_list:
	os.mkdir("xml")
if "tab" not in directory_list:
	os.mkdir("tab")


#now that we have a list of gbk like files we will make a dictionary out of them
# and rename them to something better than what NCBI does

project_dictionary = []

for gbfffile in genbank_file_list:
	dictentry = nav.dictmaker(gbfffile)
	newname = dictentry.split(" | ")[0]
	os.rename(gbfffile, str("gbk/" + newname + ".gbk"))
	project_dictionary.append(dictentry)

dictionaryfile = open(project_name + ".dict", "w")
dictionaryfile.write("\n".join(project_dictionary))
dictionaryfile.close()

#part 2 - now create fna files and put in the fna directory

for entry in os.scandir('gbk'):
	if entry.is_file():
		if str(entry.name)[0] != "." and str(entry.name.split(".")[-1]) == "gbk":
			passfile = "gbk/" + str(entry.name)
			nav.convert_gbk_to_fna(passfile)
			newname = str(entry.name.split(".")[0]) + ".fna"
			os.rename("gbk/" + newname, "fna/" + newname)
			
			
#make a blast db	
cmd = "cat fna/*.fna > fna/temp.fna"
os.system(cmd)
cmd2 = "makeblastdb -parse_seqids -in fna/temp.fna -dbtype nucl -out " + project_name
os.system(cmd2) 
cmd3 = "rm fna/temp.fna" #remove the concatenated fna file
os.system(cmd3) 

#now go through each fna file and blast against the database
fna_list = []
for entry in os.scandir('fna'):
    if entry.is_file():        	
        if str(entry.name)[0] != "." and str(entry.name.split(".")[-1]) == "fna":
        	fna_list.append(entry.name)
for fna_file in fna_list:
	xml_name = "xml/" + fna_file.split(".")[0] + ".xml"
	tab_name = "tab/" + fna_file.split(".")[0] + ".tab"
	dict_name = project_name + ".dict"
	blast_command = "blastn -outfmt 5 -query fna/" + fna_file + " -db " + project_name + " -out " + xml_name
	os.system(blast_command)
	parser.pxml(xml_name, tab_name, dict_name, pct_id)
	os.system("rm " + xml_name)
	
import pandas as pd
	