# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:03:01 2016

@author: XuGang
"""

from __future__ import division
import numpy as np

def getConservedInfo(atomsData, residuesData, db, db_cash, scoring_function, window_len):

    scores = 0 
    
    m = int((window_len-1)/2)
    for i in range(m,len(residuesData)-m):
        group = []
        
        k = m
        while(k > 0):             
            group.append(residuesData[i-k])                
            k = k - 1
            
        group.append(residuesData[i])  

        j = 1
        while(j < int((window_len+1)/2)):           
            group.append(residuesData[i+j])         
            j = j + 1  
            
        assert len(group) == window_len
        
        resids = []
        for g in group:
            resids.append(g.resid)
        
        if np.mean(resids) == group[m].resid:
            if scoring_function == "csf":
                score, db_cash = getCSFScore(group, window_len, db, db_cash)
            elif scoring_function == "ssf":
                score, db_cash = getSSFScore(group, window_len, db, db_cash)
            else:
                score, db_cash = getDASFScore(group, window_len, db, db_cash)
            scores = scores + score*window_len
        
    return scores,db_cash
        
def getCSFScore(group, window_len, db, db_cash):
    
    seq = []
    atoms = []
    features = []
    
    score = 0
    INCOMPLETE = False
    
    for indiv in group:
        seq.append(indiv.resname) 
        atoms.append(indiv.atoms) 

    if(window_len%4 == 1): #5,9
        i = 0
    else:  #7,11
        i = 1
    ref_id = int((window_len-1)/2)
    while(i < window_len):
        if(i == ref_id):
            i = i + 2
            continue

        try:
            features.append(transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,atoms[i]["C"].position))
        except:
            INCOMPLETE = True
            break
        
        i = i + 2

    if(INCOMPLETE):
        return 0, db_cash
    
    seq_str = "".join(seq)
    
    if seq_str in db_cash:
        score = calZscore1(features,seq_str,db_cash[seq_str])
    elif db.find({"key":seq_str}) is not None:
        stander_str = db.find({"key":seq_str})['value']
        standers = []
        for i in stander_str.split(';'):
            stander = []
            indiv = i.split('#')
            stander.append(np.array([float(j) for j in indiv[0].split('_')]))
            stander.append(np.array([float(j) for j in indiv[1].split('_')]))
            stander.append(float(indiv[2]))
            standers.append(stander) 
        db_cash[seq_str] = standers
        score = calZscore1(features,seq_str,db_cash[seq_str])
    else:
        score = 0
    
    return score, db_cash
    
def getDASFScore(group, window_len, db, db_cash):
    
    seq = []
    atoms = []
    features = []
    
    score = 0
    INCOMPLETE = False
    
    for indiv in group:
        seq.append(indiv.resname) 
        atoms.append(indiv.atoms) 

    if(window_len%4 == 1): #5,9
        i = 0
    else:  #7,11
        i = 1
    
    while(i < window_len):
        
        f1 = f2 = f3 = f4 = None
        try:
            if(seq[i] == 'G'):
                pass
            
            elif(seq[i] == 'A'):
                pass
            
            elif(seq[i] == 'V'):
                f1 = atoms[i]["CG1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])
            elif(seq[i] == 'I'):
                f1 = atoms[i]["CG1"].position
                try:
                    f2 = atoms[i]["CD"].position
                except:
                    f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                             
            elif(seq[i] == 'L'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                         
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                               
            elif(seq[i] == 'S'):
                f1 = atoms[i]["OG"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
            elif(seq[i] == 'T'):
                f1 = atoms[i]["OG1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                               
            elif(seq[i] == 'D'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["OD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                             
            elif(seq[i] == 'N'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["OD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                             
            elif(seq[i] == 'E'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD"].position
                f3 = atoms[i]["OE1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f3),window_len])   
            elif(seq[i] == 'Q'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD"].position
                f3 = atoms[i]["OE1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f3),window_len])    
            elif(seq[i] == 'K'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD"].position
                f3 = atoms[i]["CE"].position
                f4 = atoms[i]["NZ"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f3),window_len])   
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f4),window_len])  
            elif(seq[i] == 'R'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD"].position
                f3 = atoms[i]["NE"].position
                f4 = atoms[i]["CZ"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f3),window_len])   
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f4),window_len])                            
            elif(seq[i] == 'C'):
                f1 = atoms[i]["SG"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
            elif(seq[i] == 'M'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["SD"].position
                f3 = atoms[i]["CE"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f3),window_len])   
            elif(seq[i] == 'F'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])       
            elif(seq[i] == 'Y'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                              
            elif(seq[i] == 'W'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])       
            elif(seq[i] == 'H'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["ND1"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])       
            elif(seq[i] == 'P'):
                f1 = atoms[i]["CG"].position
                f2 = atoms[i]["CD"].position
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f1),window_len])                       
                features.append([transCoordinate(atoms[i]["CA"].position,atoms[i]["C"].position,atoms[i]["O"].position,\
                               f2),window_len])                        
            else:
                print ("Residue Name Not Found: ", seq[i])
        except:
            INCOMPLETE = True
            break
        
        i = i + 2

    if(INCOMPLETE):
        return 0, db_cash
    
    seq_str = "".join(seq)
    
    if seq_str in db_cash:
        score = calZscore2(features,seq_str,db_cash[seq_str])
    elif db.find({"key":seq_str}) is not None:
        stander_str = db.find({"key":seq_str})['value']
        standers = []
        for i in stander_str.split(';'):
            stander = []
            indiv = i.split('#')
            stander.append(np.array([float(j) for j in indiv[0].split('_')]))
            stander.append(np.array([float(j) for j in indiv[1].split('_')]))
            stander.append(float(indiv[2]))
            standers.append(stander) 
        db_cash[seq_str] = standers
        score = calZscore2(features,seq_str,db_cash[seq_str])
    else:
        score = 0
    
    return score, db_cash
    
def getSSFScore(group, window_len, db, db_cash):
    
    seq = []
    atoms = []
    features = []
    
    score = 0
    INCOMPLETE = False
    
    for indiv in group:
        seq.append(indiv.resname) 
        atoms.append(indiv.atoms) 

    if(window_len%4 == 1): #5,9
        i = 0
    else:  #7,11
        i = 1
    ref_id = int((window_len-1)/2)
    while(i < window_len):
        if(i == ref_id):
            i = i + 2
            continue
        
        try:
            if(seq[i] == 'G'):
                f1 = atoms[i]["C"].position
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])
            elif(seq[i] == 'A'):
                f1 = atoms[i]["CB"].position
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                                f1),1])
            elif(seq[i] == 'V'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG1"].position) + np.array(atoms[i]["CG2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])
            elif(seq[i] == 'I'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG1"].position) + np.array(atoms[i]["CG2"].position))/3
                try:
                    f2 = atoms[i]["CD"].position
                except:
                    f2 = atoms[i]["CD1"].position
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])                             
            elif(seq[i] == 'L'):
                f1 = atoms[i]["CB"].position
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CD1"].position) + np.array(atoms[i]["CD2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])                               
            elif(seq[i] == 'S'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["OG"].position))/2
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])                       
            elif(seq[i] == 'T'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["OG1"].position) + np.array(atoms[i]["CG2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])                               
            elif(seq[i] == 'D'):
                f1 = atoms[i]["CB"].position
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["OD1"].position) + np.array(atoms[i]["OD2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])                             
            elif(seq[i] == 'N'):
                f1 = atoms[i]["CB"].position
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["OD1"].position) + np.array(atoms[i]["ND2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])                             
            elif(seq[i] == 'E'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG"].position))/2
                f2 = (np.array(atoms[i]["CD"].position) + np.array(atoms[i]["OE1"].position) + np.array(atoms[i]["OE2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])                              
            elif(seq[i] == 'Q'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG"].position))/2
                f2 = (np.array(atoms[i]["CD"].position) + np.array(atoms[i]["OE1"].position) + np.array(atoms[i]["NE2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2])  
            elif(seq[i] == 'K'):
                f1 = atoms[i]["CB"].position
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CD"].position))/2
                f3 = (np.array(atoms[i]["CE"].position) + np.array(atoms[i]["NZ"].position))/2
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),3])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),3]) 
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f3),3]) 
            elif(seq[i] == 'R'):
                f1 = atoms[i]["CB"].position
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CD"].position))/2
                f3 = (np.array(atoms[i]["NE"].position) + np.array(atoms[i]["CZ"].position) + np.array(atoms[i]["NH1"].position) + np.array(atoms[i]["NH2"].position))/4
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),3])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),3]) 
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f3),3])                            
            elif(seq[i] == 'C'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["SG"].position))/2
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])                       
            elif(seq[i] == 'M'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG"].position))/2
                f2 = (np.array(atoms[i]["SD"].position) + np.array(atoms[i]["CE"].position))/2
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2]) 
            elif(seq[i] == 'F'):
                f1 = np.array(atoms[i]["CB"].position)
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CE1"].position) + np.array(atoms[i]["CE2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2]) 
            elif(seq[i] == 'Y'):
                f1 = np.array(atoms[i]["CB"].position)
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CE1"].position) + np.array(atoms[i]["CE2"].position))/3
                f3 = np.array(atoms[i]["OH"].position)
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),3])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),3]) 
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f3),3])                         
            elif(seq[i] == 'W'):
                f1 = np.array(atoms[i]["CB"].position)
                f2 = (np.array(atoms[i]["NE1"].position) + np.array(atoms[i]["CZ2"].position) + np.array(atoms[i]["CZ3"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2]) 
            elif(seq[i] == 'H'):
                f1 = np.array(atoms[i]["CB"].position)
                f2 = (np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CE1"].position) + np.array(atoms[i]["NE2"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),2])                       
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f2),2]) 
            elif(seq[i] == 'P'):
                f1 = (np.array(atoms[i]["CB"].position) + np.array(atoms[i]["CG"].position) + np.array(atoms[i]["CD"].position))/3
                features.append([transCoordinate(atoms[ref_id]["CA"].position,atoms[ref_id]["C"].position,atoms[ref_id]["O"].position,\
                               f1),1])                    
            else:
                print ("Residue Name Not Found: ", seq[i])
        except:
            INCOMPLETE = True
            break
        
        i = i + 2

    if(INCOMPLETE):
        return 0, db_cash
    
    seq_str = "".join(seq)
    
    if seq_str in db_cash:
        score = calZscore2(features,seq_str,db_cash[seq_str])
    elif db.find({"key":seq_str}) is not None:
        stander_str = db.find({"key":seq_str})['value']
        standers = []
        for i in stander_str.split(';'):
            stander = []
            indiv = i.split('#')
            stander.append(np.array([float(j) for j in indiv[0].split('_')]))
            stander.append(np.array([float(j) for j in indiv[1].split('_')]))
            stander.append(float(indiv[2]))
            standers.append(stander) 
        db_cash[seq_str] = standers
        score = calZscore2(features,seq_str,db_cash[seq_str])
    else:
        score = 0
    
    return score, db_cash
    
def calZscore1(features,seq, standers):    
    
    if standers[0][2] < 3:
        return 0  

    feature_num = len(features)
    results = 0
    
    i = 0
    while(i < feature_num):
        
        central_stander = standers[i][0]
        sd_stander = standers[i][1]
        
        r = np.abs((features[i]-central_stander)/sd_stander)
        r[np.where(r>5)] = 5
        result = np.sum(r)
     
        results = results + result
        i = i + 1  
        
    return results

def calZscore2(features,seq, standers):    
    
    if standers[0][2] < 3:
        return 0  

    feature_num = len(features)
    results = 0
    
    i = 0
    while(i < feature_num):
        
        central_stander = standers[i][0]
        sd_stander = standers[i][1]


        result = np.sum(np.abs((features[i][0]-central_stander)/sd_stander))
        if result > 15:
            result = 15
            
        result = result/features[i][1]
        
        results = results + result
        i = i + 1  
        
    return results
    
def transCoordinate(atom_ca_ref, atom_c_ref, atom_o_ref, atom_c):
    
    ref = atom_ca_ref
    c_ref_new = atom_c_ref - ref
    o_ref_new = atom_o_ref - ref
    c_new = atom_c - ref
  
    #c-ca
    x_axis = c_ref_new/np.linalg.norm(c_ref_new)
    
    c_o = o_ref_new - c_ref_new
    
    #o-c perpendicular to x_axis
    y_axis = c_o - (x_axis.dot(c_o)/x_axis.dot(x_axis) * x_axis)
    y_axis = y_axis/np.linalg.norm(y_axis)

    z_axis = np.cross(x_axis,y_axis)
    
    rotation_matrix = np.array([x_axis[0],y_axis[0],z_axis[0],0,x_axis[1],y_axis[1],z_axis[1],0,x_axis[2],y_axis[2],z_axis[2],0,0,0,0,1]).reshape(4,4)

    new = np.array([c_new[0],c_new[1],c_new[2],1]).dot(rotation_matrix)

    return np.array([new[0],new[1],new[2]])
    










