# -*- coding: utf-8 -*-
"""
Created on Tue Jan 07 16:21:10 2014

@author: Maxim
"""

import os
import xlrd
import xlutils.copy
import xlwt
import numpy as np
from paths import MyPaths

def read_txt_table(filePath,hasHeader=True):
    pass

def write_txt_table(filePath,hasHeader=True):
    pass

def load_xls_read(xlsPath,sheetName):
    pass

def load_xls_write(xlsPath,sheetName):
    pass


class SectionMap:
    r"""
    Creates map of given xls sheet by searching locations of cells containing 
    keyword. Useful when working with "infinite" list of data when number of 
    lines is not specified and requred data is bounded by cells containing 
    keyword.
    
    Parameters
    ----------
    dbSheet : sheet object
        xls sheet to be scanned to detect sections
    keyword : string
        keyword for map building
    
    getItem method returns *range* iterable object that contains row indexes 
    of required section.
    """
    def __init__(self,dbSheet,keyword='SECTION: '):
        self.sheet = dbSheet
        self._build_map(keyword)

    def __getitem__(self,name):
        idx = self.sectionName.index(name)
        return range(self.rowIdx[idx],self.rowIdx[idx+1]-1)

    def _build_map(self,keyword='SECTION: '):
        n           = len(keyword)
        rowIdx      = list()
        sectionName = list()
        for row in range(self.sheet.nrows):
            cellValue = str(self.sheet.row(row)[0].value)
            if cellValue[:n]==keyword:
                rowIdx.append(row+1)
                sectionName.append(cellValue[n:])
        rowIdx.append(self.sheet.nrows+1)
        self.rowIdx = rowIdx
        self.sectionName = sectionName


class LoadDatabase(object):
    def __init__(self,xlsPath,mode='r'):
        self.mode = mode
        self.filePath = xlsPath
        self.workbook = xlrd.open_workbook(self.filePath)
        self.sheetList,self.nameList = self._load_data()
        if mode=='w':
            self.workbook = xlutils.copy.copy(self.workbook)

    def _load_data(self):
        """ loads names of all datasheets and sheet objects
        """
        sheetList=list()
        nameList =list()
        for sheetName in self.workbook.sheet_names():
            worksheet = self.workbook.sheet_by_name(sheetName)
            sheetList.append(worksheet)
            nameList.append(str(sheetName))
        return sheetList,nameList

    def select_by_name(self,sheetName):
        """
        Returns xls sheet object with specified name. If sheet is not found 
        then None value is returned
        
        Parameters
        ----------

        sheetName : string
            name of xls sheet to open
        """
        try:
            sheetNum=self.nameList.index(sheetName)
            sheet=self.sheetList[sheetNum]
        except:
            print "Error: specified sheet doesn't exist"
            sheet=None
        return sheet
    
    def _check_name(self):
        """
        Optional function used to check if name of sheet to save already exists. 
        If name exists then name will be updated in format name_copyN where 
        N - number.
        Returns name if sheet with given name doesn't exist and name_copyN 
        if sheet exist.
        
        Parameters
        ----------
        name : string
        """
        def name_exists(name):
            for existName in self.nameList:
                if existName==name:
                    return True
            else:
                return False
        if name_exists(name):
            newName = name
            idx = 0
            while name_exists(newName):
                idx+=1
                newName = '%s_copy%d'%(name,idx)
        else:
            newName = name
        return newName

    def add_sheet(self, newSheetName):
        """
        Adds new sheet with specified name to workbook. To save new sheet 
        program runs check_name method first, so final name can be different 
        from newSheetName
        
        Parameters
        ----------
        
        newSheetName : string
            name of new sheet that will be added
        """
        newSheetName = self._check_name(newSheetName)
        newSheet = self.workbook.add_sheet(newSheetName)
        return newSheet
    
    def delete_sheet(self, sheetName):
        pass
    
    def save(self):
        """Saves open database. It is neccessary to save db after writing 
        data to it.
        """
        if self.mode=='w':
            self.workbook.save(self.filePath)


class ReadDatabase():
    def __init__(self,filePath,sheetName):
        pass

class WriteDatabase():
    pass

# --- debug section ---
def run_test1():
    pth = MyPaths()
    sh = ReadDatabase(pth.db.aircraft,'V0510')

if __name__=="__main__":
    run_test1()