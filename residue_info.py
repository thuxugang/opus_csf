# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:47:45 2016

@author: XuGang
"""

class Residue:
    def __init__(self, resid, resname):
        #res_id 178A..
        if resid[-1].isalpha():
            self.resid = -100
        else:
            self.resid = int(resid)      
        self.resname = resname
        self.atoms = {}
        
def getResidueData(atomsData):
    
    residuesData = []
    last_resid = None
    
    for atom in atomsData:
        #first
        if(atom == atomsData[0]):           
            residue = Residue(atom.resid,atom.resname) 
            residue.atoms[atom.name1] = atom
            
        #last
        elif(atom == atomsData[-1]):
            
            if(last_resid == atom.resid):                
                residue.atoms[atom.name1] = atom  
                residuesData.append(residue)
            else:                
                residue = Residue(atom.resid,atom.resname) 
                residue.atoms[atom.name1] = atom
                residuesData.append(residue)
                
        else:
            if(last_resid == atom.resid):                
                residue.atoms[atom.name1] = atom          
            else:
                residuesData.append(residue)
                residue = Residue(atom.resid,atom.resname) 
                residue.atoms[atom.name1] = atom                  
        
        last_resid = atom.resid
    
    return residuesData

def singleResname(AA):
    if(len(AA) == 1):
        return AA
    else:
        if(AA in ['GLY','AGLY']):
            return "G"
        elif(AA in ['ALA','AALA']):
            return "A"
        elif(AA in ['SER','ASER']):
            return "S"
        elif(AA in ['CYS','ACYS']):
            return "C"
        elif(AA in ['VAL','AVAL']):
            return "V"
        elif(AA in ['ILE','AILE']):
            return "I"
        elif(AA in ['LEU','ALEU','LEF','ALEF']):
            return "L"
        elif(AA in ['THR','ATHR']):
            return "T"
        elif(AA in ['ARG','AARG']):
            return "R"
        elif(AA in ['LYS','ALYS']):
            return "K"
        elif(AA in ['ASP','AASP']):
            return "D"
        elif(AA in ['GLU','AGLU']):
            return "E"
        elif(AA in ['ASN','AASN']):
            return "N"
        elif(AA in ['GLN','AGLN']):
            return "Q"
        elif(AA in ['MET','AMET']):
            return "M"
        elif(AA in ['HIS','AHIS']):
            return "H"
        elif(AA in ['PRO','APRO']):
            return "P"
        elif(AA in ['PHE','APHE']):
            return "F"
        elif(AA in ['TYR','ATYR', "PTR"]):
            return "Y"
        elif(AA in ['TRP','ATRP']):
            return "W"
        else:
            pass
            #print "Residue.singleResname() false" + AA

def triResname(AA):
    if(len(AA) == 3):
        return AA
    else:
        if(AA == "G"):
            return "GLY"
        elif(AA == "A"):
            return "ALA"
        elif(AA == "S"):
            return "SER"
        elif(AA == "C"):
            return "CYS"
        elif(AA == "V"):
            return "VAL"
        elif(AA == "I"):
            return "ILE"
        elif(AA == "L"):
            return "LEU"
        elif(AA == "T"):
            return "THR"
        elif(AA == "R"):
            return "ARG"
        elif(AA == "K"):
            return "LYS"
        elif(AA == "D"):
            return "ASP"
        elif(AA == "E"):
            return "GLU"
        elif(AA == "N"):
            return "ASN"
        elif(AA == "Q"):
            return "GLN"
        elif(AA == "M"):
            return "MET"
        elif(AA == "H"):
            return "HIS"
        elif(AA == "P"):
            return "PRO"
        elif(AA == "F"):
            return "PHE"
        elif(AA == "Y"):
            return "TYR"
        elif(AA == "W"):
            return "TRP"
        else:
            print ("residue_info.triResname() false" + AA)