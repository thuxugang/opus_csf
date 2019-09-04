# -*- coding: utf-8 -*-
"""
Created on Sat May 30 07:14:18 2015

@author: XuGang
"""

import numpy as np
import residue_info

class Atom:
    def __init__(self, id, name1, resname, resid, position):
        self.id = id
        self.name1 = name1
        self.resname = residue_info.singleResname(resname)
        self.resid = resid
        self.position = np.array([position[0],position[1],position[2]])
 
