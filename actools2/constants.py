# -*- coding: utf-8 -*-
"""
Class for loading global constants from the excel constants database file
Provides functions for loading cosntant values.
"""

import paths
import dbTools

g = 9.81

def load(sheetName,overrideColumn=0):
    const = Constants()
    const.load(sheetName,overrideColumn)
    return const


class Constants():
    def __init__(self):
        pth=paths.Database()
        self.db=dbTools.loadDB(pth.constants)
        self.data=list()
        self.dict=dict()
        self._headers=1
    def load(self,sheetName,overrideColumn=0):
        data=self.db.readSheet(sheetName)
        self.data=data
        for i in range(1,len(self.data),1):
            row=data[i]
            try:
                value=float(row[overrideColumn+self._headers])
            except:
                value=float(row[self._headers])
            key=str(row[0])
            self.dict.update([(key,value)])
    def listConstants(self):
        for key in self.dict.keys():print key
    def getValue(self,name):
        value=None
        try:
            value=self.dict.get(name)
        except:
            print('error in constants(): entry not found')
        return value

def run_test1():
    const = Constants()
    const.load('mass',0)
    const.listConstants()
    print const.getValue('fuseCGratio')

if __name__=='__main__':
    run_test1()
