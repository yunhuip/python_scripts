# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import pandas as pd 
import requests
import re
import csv

## Retrive histone sequence by uniprotID and perform MSA using clustal-omega
def uniprotID_to_fa_msa(uniprotID_file, fa_file):
    uniprot_ID = pd.read_table(uniprotID_file,header = None)
    with open(fa_file, "w") as fh:
        for i in range(0, len(uniprot_ID)):    
            response = requests.post("https://www.uniprot.org/uniprot/" + uniprot_ID.iloc[i][0] + ".fasta")
            fh.write(response.text) 
    os.system("./clustal-omega-1.2.3-macosx -i " + fa_file + " -o " + fa_file[0:-3] + "_msa.fa --auto -v --force")

def joint_list_to_string(input_list,sting_seq):       
    # Converting integer list to string list 
    s = [str(i) for i in input_list]       
    # Join list items using join() 
    res = str(sting_seq.join(s)) 
    return(res)

## read aligned histone sequences from msa file
def read_prot_sequences(prot_seq_file):
    prot_seq = {}     # dictionary to store all histone sequences from  MSA
    with open(prot_seq_file, 'r') as reader:  
        for lines in (reader.read().split(">")[1:]):
            uniprotID = lines.split("\n")[0].split("|")[1].lower()
            aligned_seq = joint_list_to_string(lines.split("\n")[1:],"")
            prot_seq[uniprotID] = aligned_seq
    return prot_seq 
        
### mapping binding site in the aligned histone sequences and count binding partners per sites   
def map_binding_sites_to_msa(binding_sites_file, prot_seq_file, output_resi_map_per_his,output_resi_map_consensus): 
#    binding_sites_file = "../raw_binding_mutation/PDB_H2A_binding_sites.txt"
#    prot_seq_file = "combined_H2A_msa.fa"
    
    ## read data files contain the binding sites sites informations 
    df = pd.read_table(binding_sites_file, skiprows = 1, header=None,     \
                          names = ["histone_uniprot_ID", "histone_type", "binding_site_index", \
                                   "partner_uniprot_ID", "sources"]) 
    #read aligned histone sequences from msa file
    all_msa_seq = read_prot_sequences(prot_seq_file) 
    
    ## map old residue index to new index
    new_binding_sites_msa = pd.DataFrame() # dataframe to save new binding sites indexes on aligned sequences
    for i in range(0,len(df)):
        his_uniprot_ID = df["histone_uniprot_ID"][i]
        old_bindsite_index = df["binding_site_index"][i]
        BP_uniprot_ID = df["partner_uniprot_ID"][i]
        data_source = df["sources"][i]
        
        ## get the aligned seq and original seq for a histone
        seq_aligned = all_msa_seq[his_uniprot_ID]
        seq_original = seq_aligned.replace('-', '')
        
        ## find new index in aligned seq
        for resdi in range(old_bindsite_index, len(seq_aligned) +1 ):
            if(seq_original[0:old_bindsite_index] == seq_aligned[0:resdi].replace('-', '')):
                new_bindsite_index = resdi
                break
         ## save new residue index to dataframe
        if data_source == "PDB" :
            new_binding_sites_msa = new_binding_sites_msa.append({"his_uniprot_ID": his_uniprot_ID, \
            "new_bindsite_index": new_bindsite_index,"partner_uniprot_ID": BP_uniprot_ID}, ignore_index=True)
                     
        elif data_source == "XLMS":
            new_binding_sites_msa = new_binding_sites_msa.append({"his_uniprot_ID": his_uniprot_ID, \
            "new_bindsite_index": new_bindsite_index -1,"partner_uniprot_ID": BP_uniprot_ID}, ignore_index=True)
            new_binding_sites_msa = new_binding_sites_msa.append({"his_uniprot_ID": his_uniprot_ID, \
            "new_bindsite_index": new_bindsite_index,"partner_uniprot_ID": BP_uniprot_ID}, ignore_index=True)
            new_binding_sites_msa = new_binding_sites_msa.append({"his_uniprot_ID": his_uniprot_ID, \
            "new_bindsite_index": new_bindsite_index+1,"partner_uniprot_ID": BP_uniprot_ID}, ignore_index=True)

    ## count number of binding partners per residue per histone variant
    with open(output_resi_map_per_his, 'w') as myfile:
        all_his_uniprotIDs= np.unique(new_binding_sites_msa["his_uniprot_ID"])
        
        for his_uniprotID in all_his_uniprotIDs:
            his_seq_binding_site_map = [0]*len(seq_aligned)
            his_seq_binding_site_map.insert(0,his_uniprotID)
            
            
            all_bindind_site = new_binding_sites_msa[new_binding_sites_msa["his_uniprot_ID"] == his_uniprotID]["new_bindsite_index"]
                   
            bind_sites_index = np.unique(all_bindind_site,return_counts=True)[0]
            bind_sites_BP_count = np.unique(all_bindind_site,return_counts=True)[1]
            
            for res_idex, bp_count  in zip(bind_sites_index, bind_sites_BP_count):
                his_seq_binding_site_map[int(res_idex)] = int(bp_count)
                
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(his_seq_binding_site_map)
            
   ## count number of binding partners per residue on  consensus sequence    
    with open(output_resi_map_consensus, 'w') as myfile:

        binding_sites_consensus = new_binding_sites_msa.drop("his_uniprot_ID",axis=1).drop_duplicates().reset_index()
        
        for resi in range(1,len(seq_aligned)+1):
            bp_count = 0
            for i in range(0,len(binding_sites_consensus)):
                if (resi == binding_sites_consensus["new_bindsite_index"][i]):
                    bp_count +=1
            myfile.write(str(resi) + "," + str(bp_count) + '\n' )  
        
if __name__=="__main__":
    uniprotID_to_fa_msa("../raw_binding_mutation/combined_uniprot_H1.txt", "combined_H1.fa")
    uniprotID_to_fa_msa("../raw_binding_mutation/combined_uniprot_H2A.txt", "combined_H2A.fa")
    uniprotID_to_fa_msa("../raw_binding_mutation/combined_uniprot_H2B.txt", "combined_H2B.fa")
    uniprotID_to_fa_msa("../raw_binding_mutation/combined_uniprot_H3.txt", "combined_H3.fa")
    uniprotID_to_fa_msa("../raw_binding_mutation/combined_uniprot_H4.txt", "combined_H4.fa")
    
    map_binding_sites_to_msa("../raw_binding_mutation/PDB_H1_binding_sites.txt","combined_H1_msa.fa", './test/PDB_H1_mapping_binding_sites.txt',"./test/H1_binding_sites_consens_counts_PDB_unique.txt")   
    map_binding_sites_to_msa("../raw_binding_mutation/PDB_H2A_binding_sites.txt","combined_H2A_msa.fa", './test/PDB_H2A_mapping_binding_sites.txt',"./test/H2A_binding_sites_consens_counts_PDB_unique.txt")   
    map_binding_sites_to_msa("../raw_binding_mutation/PDB_H2B_binding_sites.txt","combined_H2B_msa.fa", './test/PDB_H2B_mapping_binding_sites.txt',"./test/H2B_binding_sites_consens_counts_PDB_unique.txt")   
    map_binding_sites_to_msa("../raw_binding_mutation/PDB_H3_binding_sites.txt","combined_H3_msa.fa", './test/PDB_H3_mapping_binding_sites.txt',"./test/H3_binding_sites_consens_counts_PDB_unique.txt")   
    map_binding_sites_to_msa("../raw_binding_mutation/PDB_H4_binding_sites.txt","combined_H4_msa.fa", './test/PDB_H4_mapping_binding_sites.txt',"./test/H4_binding_sites_consens_counts_PDB_unique.txt")   
    
