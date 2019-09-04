# -*- coding: utf-8 -*-
"""
Created on Fri May 29 18:32:13 2015

@author: XuGang
"""

'''
ATOM   1   N     SER     A    9    9.932  76.723  40.808  1.00 47.18  N
pass   id  name1 resname pass pass x      y       z       pass pass   name2
self.rigidBodyBlocks = belongs to which rigidBodyBlocks
self.name = id in rigidBodyBlocks r1~r6
only uss A conformation
'''

import atom_info
import numpy as np    

def readPDB(filename):
    f = open(filename,'r')
    atomsDatas = []
    atomsData = []
    count = 0
    while True:
        count = count + 1
        line = f.readline()   
        if (line == "" or line == "\n"):
            atomsDatas.append(atomsData)
            atomsData = []
            break
        elif(line.split()[0] == "TER"):
            atomsDatas.append(atomsData)
            atomsData = []
        else:
            if (line.split()[0] == 'ATOM' ):
                temp = list(line)
                id = ""
                for i in [6,7,8,9,10]:
                    if(str(temp[i]) != " "):
                        id = id + temp[i]
        
                name1 = ""
                for i in [11,12,13,14,15]:
                    if(str(temp[i]) != " "):
                        name1 = name1 + temp[i]
    
                resname = ""
                for i in [16,17,18,19]:
                    if(str(temp[i]) != " "):
                        resname = resname + temp[i]
                        
                #B confomation        
                if(len(resname) == 4 and resname[0] != "A"):
                    continue
                
                resid = ""
                for i in [22,23,24,25,26]:
                    if(str(temp[i]) != " "):
                        resid = resid + temp[i] 
                        
                x = ""
                for i in [30,31,32,33,34,35,36,37]:
                    if(str(temp[i]) != " "):
                        x = x + temp[i]
        
                y = ""
                for i in [38,39,40,41,42,43,44,45]:
                    if(str(temp[i]) != " "):
                        y = y + temp[i]
        
                z = ""
                for i in [46,47,48,49,50,51,52,53]:
                    if(str(temp[i]) != " "):
                        z = z + temp[i]                        
                if(name1[0] in ["N","O","C","S"]):
                    position = np.array([float(x), float(y), float(z)])
                    atom = atom_info.Atom(id, name1, resname, resid, position)
                    atomsData.append(atom)
    #atomsData = sorted(atomsData,key = lambda x:int(x.id))
#    if(lost):print "lost atom data in", filename
    f.close()
    return atomsDatas