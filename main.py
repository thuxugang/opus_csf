# -*- coding: utf-8 -*-

"""

Created on Wed May 11 08:43:50 2016

@author: Xu Gang

"""

from __future__ import division
import pdb_reader
import residue_info
import analyze
import mongoDB


#########################################################
#  Parameters  #
#########################################################

scoring_function = "csf" #(or ”ssf" or "dasf")

db_name = "csf_db"
collection_name_5 = "csf_5"
collection_name_7 = "csf_7"
collection_name_9 = "csf_9"
collection_name_11 = "csf_11"

list_path = "./list_path.txt"

#########################################################

if __name__ == "__main__":
    
    print ("Using Scoring Function: ", scoring_function)
    
    assert scoring_function in ["csf", "dasf", "ssf"]
    
    db5 = mongoDB.getConnect(db_name,collection_name_5)
    db7 = mongoDB.getConnect(db_name,collection_name_7)
    db9 = mongoDB.getConnect(db_name,collection_name_9)
    db11 = mongoDB.getConnect(db_name,collection_name_11)
    
    db_cash = {}
    files_path = []

    with open(list_path, 'r') as f:
        contents = f.readlines()
        for content in contents:
            files_path.append(content.strip())
    
    for file_path in files_path:
        atomsDatas = pdb_reader.readPDB(file_path) 
        scores = 0
        for atomsData in atomsDatas:
            residuesData = residue_info.getResidueData(atomsData)      
            score5,db_cash = analyze.getConservedInfo(atomsData,residuesData,db5,db_cash,scoring_function,5)
            score7,db_cash = analyze.getConservedInfo(atomsData,residuesData,db7,db_cash,scoring_function,7)
            score9,db_cash = analyze.getConservedInfo(atomsData,residuesData,db9,db_cash,scoring_function,9)
            score11,db_cash = analyze.getConservedInfo(atomsData,residuesData,db11,db_cash,scoring_function,11)
            scores = scores + score5 + score7 + score9 + score11
    
        print (file_path, scores)

                
                

                    
                                       