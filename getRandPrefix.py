# -*- coding: utf-8 -*-
"""
Created on Wed May 16 23:39:51 2012

@author: Maxim
"""

def getRandPrefix(fileExt = "", addSymbol = ""):
    from random import randrange
    from time import gmtime, strftime
    
    Time = int(strftime("%H%M%S", gmtime()))
    NamePrefix = str(Time+randrange(0,1e6,1)) + addSymbol
    if fileExt != "":
        NamePrefix  = NamePrefix + "." + fileExt
    
    return(NamePrefix)
    
    
h = getRandPrefix()
print(h)