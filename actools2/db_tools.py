# -*- coding: utf-8 -*-
"""
Created on Tue Jan 07 16:21:10 2014

Provides easy interface to XLS and TXT format data.

@author: Maxim
"""

import xlrd
import xlutils.copy
import numpy as np
from paths import MyPaths

# --- text files ---

def read_txt_table(filePath,hasHeader=True):
    pass

def write_txt_table(filePath,hasHeader=True):
    pass

# --- xls operations ---
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
    
    def _check_name(self, name):
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
    """
    This class provides easy interface for xls database reading. 
    
    Parameters
    ----------
    
    xlsPath : string
        path to xls file
    sheetName : string
        name of workbook sheet to be read
    secKeyword : string
        name of the keyword to create section map. If not specified then 
        no section map will be created.
    """
    def __init__(self,xlsPath,sheetName,secKeyword=None):
        db = LoadDatabase(xlsPath,'r')
        self._inputSheet = db.select_by_name(sheetName)
        self._irowPrev = -1
        self.nrows = self._inputSheet.nrows
        self.ncols = self._inputSheet.ncols
        if not secKeyword==None:
            self.sectionMap = SectionMap(self._inputSheet,secKeyword)
    
    def find_header(self,header):
        """
        Finds specific header in a first row of the open sheet and 
        returns row number of this header if exists otherwise returns -1
        """
        found = False
        out = list()
        for irow in range(self._inputSheet.nrows):
            if str(self._inputSheet.row(irow)[0].value)==header:
                out.append(irow)
                found = True
        if not found:
            print "Error: specified header not found"
        if len(out)==1:
            return out[0]
        else:
            return out

    def read_row(self,rowIdx=-1,colIdx=0,iterable=False):
        """
        Reads row of data specified by row index and first column.
        
        Parameters
        ----------
        
        rowIdx : integer
            row index. If rowIdx=-1 then next row from previous function call 
            will be used.
        colIdx : integer
            first column index. To skip the first column with parameters 
            description colIdx should be 1.
        iterable : bool
            if row contains single value then single value will be returned if 
            iterable=False. If iterable=True then array with single value will 
            be returned
        """
        row = self._get_row_by_index(rowIdx)
        rowValues = self._filter_cells(row[colIdx:])
        output = np.array(rowValues)
        if len(output)==1 and not iterable:
            output = output[0]
        return output
    
    def read_column(self,colIdx,rowIdx=0,numRow=1,iterable=False):
        """
        Reads column of data specified by column index, starting row index and 
        number of rows to be read
        
        Parameters
        ----------
        
        colIdx : integer
            column index
        rowIdx : integer
            starting row index
        numRow : integer
            number of rows to be read
        """
        col = self._inputSheet.col_slice(colIdx,rowIdx,rowIdx+numRow)
        colValues = self._filter_cells(col)
        output = np.array(colValues)
        if len(output)==1 and not iterable:
            output = output[0]
        return output
    
    def read_row_range(self,startRow,startCol,numRow):
        """
        Reads range of rows. Useful for tabulated data: for example airfoil 
        C81 tables. 
        
        Parameters
        ----------
        
        startRow : integer
            starting row index
        startCol : integer
            starting column index
        numRow : integer
            number of rows to read
        
        Note
        ----
        
        All values of the table to be read should have same data type 
        (float, string etc.)
        """
        line1 = self.read_row(startRow,startCol,True)
        output = np.zeros([numRow,len(line1)])
        output[0] = line1
        for i in range(numRow-1):
            output[i+1] = self.read_row(startRow+i+1,startCol,True)
        return output
    
    def read_section(self,sectionName,startCol=1):
        """\
        Creates map of given input sheet using SectionMap and reads all data 
        in given section starting from startCol column
        
        Paramters
        ---------
        
        sectionName : string
            name of the section to be read
        startCol : integer
            index of the first column
        """
        data = list()
        for irow in self.sectionMap[sectionName]:
            data.append(self.read_row(irow,startCol))
        return data
    
    def read_cell(self,rowIdx,colIdx):
        return self._inputSheet.cell(rowIdx,colIdx).value

    def _filter_cells(self,cells):
        values = list()
        for cell in cells:
            if not cell.value=='':
                values.append(cell.value)
        return values
    
    def _get_row_by_index(self,rowIdx):
        if rowIdx==-1:
            self._irowPrev += 1
            rowIdx = self._irowPrev
        else:
            self._irowPrev = rowIdx
        return self._inputSheet.row(rowIdx)


class WriteDatabase():
    def __init__(self,filePath,sheetName):
        self.db = LoadDatabase(filePath,'w')
        self._irowPrev = -1
        self._outputSheet = self.db.add_sheet(sheetName)
    
    def save(self):
        self.db.save()
    
    def write_row(self,label='',data=None,rowIdx=-1,colIdx=0):
        """
        Writes label to the first cell of the row and data starting from second.
        If rowIdx = -1 then data will be written to the next row from 
        previous function call
        
        Parameters
        ----------
        
        label : string
            value of the first cell in row
        data : array
            data that will be written into the cells next by label cell
        rowIdx : integer
            index of the row the data will be written to
        colIdx : integer
            index of the first column the data will be written to
        """
        rowIdx = self._get_next_row_index(rowIdx)
        icol = colIdx
        if not label=='':
            self._outputSheet.write(rowIdx,icol,label)
            icol += 1
        if not data==None:
            if hasattr(data,'__iter__'):
                for i,cellValue in enumerate(data):
                    self._outputSheet.write(rowIdx,icol+i,cellValue)
            else:
                self._outputSheet.write(rowIdx,icol,data)
    
    def write_column(self,data,startRow=-1,colIdx=0):
        """
        Writes data to the column specified
        """
        startRow = self._get_next_row_index(startRow)
        if hasattr(data,'__iter__'):
            for i,cellValue in enumerate(data):
                self._outputSheet.write(i+startRow,colIdx,cellValue)
                self._irowPrev += 1
        else:
            self._outputSheet.write(startRow,colIdx,data)
    
    def write_row_range(self,data,startRow=-1,startCol=0):
        """
        Writes data to the range of rows
        """
        for line in data:
            self.write_row('',line,-1,1)
#        startRow = self._get_next_row_index(startRow)
#        if hasattr(data,'__iter__'):
#            for irow,line in enumerate(data):
#                self._irowPrev += 1
#                for icol,cellValue in enumerate(line):
#                    self._outputSheet.write(startRow+irow,startCol+icol,cellValue)
#        else:
#            self._outputSheet.write(startRow,startCol,data)
    
    def _get_next_row_index(self,rowIdx):
        if rowIdx==-1:
            self._irowPrev += 1
            rowIdx = self._irowPrev
        else:
            self._irowPrev = rowIdx
        return rowIdx

# --- debug section ---
def run_test1():
    pth = MyPaths()
    sh1 = ReadDatabase(pth.db.airfoil,'Clark-Y')
    data = sh1.read_row_range(10,1,35)
    print data

if __name__=="__main__":
    run_test1()