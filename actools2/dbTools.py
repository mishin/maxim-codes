# -*- coding: utf-8 -*-
"""
Set of functions to communicate with xls database.
This set provides easy way to read and write different types of data as:
    
    1) single cell
    2) row
    3) column
    4) range
    5) section specified by Keyword
"""
import os
import xlrd
import xlutils.copy
import xlwt
import numpy
import paths

def read_txt_data(filePath):
    from numpy import vstack,zeros
    fid = open(filePath,'rt')
    header = fid.readline()
    del header
    line = fid.readline()
    segLine = line.split()
    n = len(segLine)
    data = zeros(n)
    for i in range(n):
        data[i] = float(segLine[i])
    for line in fid:
        if line.strip()!='':
            segLine = line.split()
            newData = zeros(n)
            for i in range(n):
                newData[i] = float(segLine[i])
            data = vstack([data,newData])
    fid.close()
    return data

def write_text_data(filePath,data,header=[]):
    fid = open(filePath,'wt')
    if header!=[]:
        fid.write(header+'\n')
    for line in data:
        for entry in line:
            fid.write('%.6f\t'%entry)
        fid.write('\n')
    fid.close()

def read_xls_doe(filePath,sheetName):
    workbook = loadDB(filePath)
    sheet = workbook.selectByName(sheetName)
    nrow = sheet.nrows
    db = readDB(sheet)
    data = db.readRange(0,0,nrow)
    return data

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
    
    getItem method returns *range* iterable object that contains row indexes 
    of required section.
    """
    def __init__(self,dbSheet):
        self.sheet = dbSheet
        self._build_map()
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

def read(dbPath,sheetName):
    workbook = loadDB(dbPath)
    sheet = workbook.selectByName(sheetName)
    db = readDB(sheet)
    return db

class loadDB:
    """
    opens database file (workbook) and creates connection to write/read data.
    
    Parameters
    ----------
    fullPath : string
        path to *.xls* database
    mode : string
        Database file open mode. *mode='r'* - open file for reading, 
        *mode='w'* - open file for writing (saving).
    """
    def __init__(self,fullPath,mode='r'):
        self.mode = mode
        self.pathName=fullPath #=os.getcwd()+relativePath
        try:
            self.workbook =xlrd.open_workbook(self.pathName)
            self.sheetList,self.nameList=self.loadData()
            if mode=='w':
                self.workbook = xlutils.copy.copy(self.workbook)
        except:
            print "Error: database file not found "
    def loadData(self):
        sheetList=list()
        nameList =list()
        for sheetName in self.workbook.sheet_names():
            worksheet = self.workbook.sheet_by_name(sheetName)
            sheetList.append(worksheet)
            nameList.append(str(sheetName))
        return sheetList,nameList
    def selectByName(self,sheetName):
        """
        Opens specified xls sheet
        
        Parameters
        ----------
        sheetName : string
            name of xls sheet to open
        
        Returns
        -------
        sheet : xls sheet object
            xls sheet with specified name
        """
        try:
            sheetNum=self.nameList.index(sheetName)
            sheet=self.sheetList[sheetNum]
        except:
            print "Error: specified sheet doesn't exist"
            sheet=list()
        return sheet
        
    def readSheet(self,sheetName):
        sheet=self.selectByName(sheetName)
        rd=readDB(sheet)
        rowList=list()
        for i,row in enumerate(range(sheet.nrows)):            
            row=rd.readRowAsList(i,1)
            rowList.append(row)
        return rowList
        
    def selectFromList(self):
        """
        Opens xls sheet specified by keyboard input
        
        Returns:
            xls sheet object
        """
        i=0
        for name in self.nameList:
            i=i+1
            print "[",i,"]  ",name
        selection=int(raw_input("Enter selection >"))
        try:        
            sheet=self.sheetList[selection-1]
        except:
            print "Invalid selection"
            sheet=list()
        return sheet
    def check_name(self,name):
        """
        Optional function used to check if name of sheet to save already exists. 
        If name exists then name will be updated in format name_copyN where 
        N - number
        
        Parameters
        ----------
        name : string
        
        Returns:
            name if if doesn't exist and name_copyN if exist.
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
    def add_sheet(self,newSheetName):
        """
        Adds new sheet with specified name to workbook. To save new sheet 
        program runs check_name method first, so final name can be different 
        from newSheetName
        
        :param string newSheetName:
            name of new sheet that will be added
        """
        newSheetName = self.check_name(newSheetName)
        newSheet = self.workbook.add_sheet(newSheetName)
        return newSheet

    def delete_sheet(self,sheetName):
        self.__init__(self.pathName,'r')
        names = self.workbook.sheet_names()
        del names[names.index(sheetName)]
        newWb = xlwt.Workbook()
        
        
        self.workbook
    def save_db(self):
        """
        Saves open database. It is neccessary to save db after writing 
        data to it.
        """
        if self.mode=='w':
            self.workbook.save(self.pathName)

class readDB:
    """
    This class is used to read in data from xls sheets.
    
    :param sheet inputSheet:
        xls sheet object that can be received from loadDB.selectByName method
    """
    def __init__(self,inputSheet):
        self._inputSheet=inputSheet
        self.sectionMap = []
    def findHeader(self,header):
        """
        Finds specific header in a first row of the open sheet and 
        returns row number of this header if exists otherwise returns -1
        
        :param string header:
            string to be found in first column of the sheet
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
    def readRow(self,rowNum,headerNum,iterable=False):
        """
        Reads specified row starting from headerNum to the end of the row. 
        Returns float or string list depending on content.

        :param integer rowNum:
            row number to read
        :param integer headerNum:
            number of first column
        :param bool iterable:
            if row contains single cell of data will return iterable 
            object (array or list)
       """
        row=self._inputSheet.row(rowNum)
        values=list()
        for cell in row[headerNum:]:
            cellValue = cell.value
            if cellValue<>str(""):
                values.append(cellValue)
        index = len(values)
        if type(values[0])==type(float()) and index>1:
            output=numpy.array(values[0:index])
        elif type(values[0])==type(float()) and index==1:
            if iterable:
                output=numpy.array([float(values[0])])
            else:
                output=float(values[0])
        elif type(values[0])==type(unicode()) and index>1:
            output=values[0:index]
        elif type(values[0])==type(unicode()) and index==1:
            if iterable:
                output=[values[0]]
            else:
                output=values[0]
        else:
            output=list()
            print "Error: input row not recognized on row ",rowNum
        return output

    def readRowAsArray(self,rowNum,headerNum):
        row=self._inputSheet.row(rowNum)
        values=list()
        for cell in row[headerNum:]:
            cellValue = cell.value
            if cellValue<>str(""):
                values.append(cellValue)
        return numpy.array(values)
    def readRowAsList(self,rowNum,headerNum):
        row=self._inputSheet.row(rowNum)
        values=list()
        for cell in row[headerNum:]:
            cellValue = cell.value
            if cellValue<>str(""):
                values.append(cellValue)
        return values
    def readRange(self,startRow,startCol,numRow):
        """
        Reads range of cells and returns array of float numbers. 
        Useful for reading tabulated data for example airfoil aerodynamic 
        tables.
        
        :param integer startRow:
            first row of the range.
        :param integer startCol:
            first column of the range
        :param integer numRow:
            number of rows to read
        """
        values = self.readRow(startRow,startCol)
        for i in range(startRow+1,startRow+numRow):
            values = numpy.vstack([values,self.readRow(i,startCol)])
        return numpy.array(values)
    def readCol(self,startRow,colNum,numRow):
        """
        Reads column of data and returns array of float numbers.
        
        Parameters
        ----------
        startRow : integer
            start row
        colNum : integer
            column number
        numRow : integer
            number of rows to read
        """
        column = self._inputSheet.col_slice(colNum,startRow,startRow+numRow)
        values = list()
        for cell in column:
            cellValue = cell.value
            if cellValue<>str(""):
                values.append(cellValue)
        if type(values[0])==type(float()):
            output = numpy.array(values)
        else:
            output = values
        return output

    def read_section(self,sectionName,startCol=1):
        r"""
        creates map using SectionMap and reads all data in given section 
        starting from startCol column.
        
        :param integer startCol:
            first column to read
        
        Returns:
            list of data read from db
            
        Note:
            data type check should be perfromed to be sure that data type is 
            right
        """
        if self.sectionMap==[]:
            self.sectionMap = SectionMap(self._inputSheet)
        data = list()
        for irow in self.sectionMap[sectionName]:
            data.append(self.readRow(irow,startCol))
        return data
     
class writeDB:
    """
    This class is used to write data to xls sheets. Input sheet can be created 
    using loadDB.add_sheet or use existing sheet. NOTE in case of using 
    existing sheet attempt to overwrite existing data will cause an error.
    
    :param sheet inputSheet:
        xls sheet object to write into
    """
    def __init__(self,inputSheet):
        self._inputSheet=inputSheet
        self._prevRow = -1
    def _row_num(self,rowNum):
        if rowNum<0:
            self._prevRow += abs(rowNum)
            return self._prevRow
        else:
            self._prevRow = rowNum
            return rowNum
    def writeRow(self,label="",data=[],rowNum=-1,startCol=0):
        """
        Writes label to the first cell of the row and data starting from 
        second cell in row. 
        
        If rowNum is negative will write the data into row 
        following the last written row with offset of abs(ronNum). 
        For example if rowNum=-1 will write data into next row after previous 
        write operation.

        :param string label:
            string to write into first cell of the row
        :param array data: 
            array of data to be written. Data may be string or float 
            (NOTE: numpy.float64 type data is not supported)
        :param rowNum:
            number of row to write data. If rowNum=-1 writes data to the 
            row next after previous use of this function
        :param integer startCol:
            first column to write data
        """
        j = startCol
        rowNum = self._row_num(rowNum)
        if label<>"":
            self._inputSheet.write(rowNum,j,label)
            j += 1
        if data<>[]:
            if type(data)==type(str()):
                self._inputSheet.write(rowNum,j,data)
            elif hasattr(data,'__iter__'):
                for i, cellVal in enumerate(data):
                    self._inputSheet.write(rowNum,i+j,cellVal)
            else:
                self._inputSheet.write(rowNum,j,data)
    def writeCol(self,data,startRow=-1,colNum=0):
        """
        Writes data to specified column starting from startRow.

        :param array data: array of data to be written in a column
        :param integer startRow: row number of first cell to write
        :param integer colNum: column number to write the data
        """
        startRow = self._row_num(startRow)
        if hasattr(data,'__iter__'):
            for i,cellVal in enumerate(data):
                self._inputSheet.write(i+startRow,colNum,cellVal)
            self._prevRow +=i
        else:
            self._inputSheet.write(startRow,colNum,data)
    def writeRange(self,data,startRow=-1,startCol=0):
        """
        Writes 2D array into xls sheet.

        :param array data: array of data to write
        :param integer startRow: first row to write data
        :param integer startCol: first column to write data
        """
        startRow = self._row_num(startRow)
        for i,row in enumerate(data):
            for j,cellVal in enumerate(row):
                self._inputSheet.write(i+startRow,j+startCol,cellVal)
        self._prevRow += i

def run_test():
    dbPath = '/database/airfoil3.xls'
    data = numpy.ones([4,5]) * 7.5
    db = loadDB(dbPath,mode='w')
    sheet = db.add_sheet('new airfoil')
    sh = writeDB(sheet)
    sh.writeRow('testlabel',data)
    db.save_db()

def run_test2():
    dbPath = paths.Database().aircraft
    db = loadDB(dbPath)
    sh = db.selectByName('sampleInput2')
    sh = readDB(sh)
    data = sh.read_section('PROPULSION',1)
    for line in data:
        print line
    print sh.readRowAsArray(1,1)

def run_test3():
    header = 'hello'
    filePath = r'D:\light aircraft\V05\propeller analysis\new 2.txt'
    filePath1 = r'D:\light aircraft\V05\propeller analysis\new 1.txt'
    data = read_txt_data(filePath)
    write_text_data(filePath1,data,header)

def debug1():
    filePath = r'D:\codes\actools\pyAC\actools\database\prop_airfoil2.xls'
    db = loadDB(os.path.abspath(filePath),'w')
    db.delete_sheet('neuform1')

def debug2():
    filePath = 'flapDoE.xls'
    data = read_xls_doe(filePath,'Sheet1')
    print data
    
if __name__=="__main__":
    debug2()
