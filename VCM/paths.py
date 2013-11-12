# -*- coding: utf-8 -*-
"""
Created on Mon Jun 04 10:25:42 2012

@author: Maxim
"""
import os
class MyPaths:
    def __init__(self):
        try:
            with open('javapath.txt','rt') as javapath: pass
        except IOError:
            print 'Error: The root folder \n['+os.getcwd()+']\nmust contain a file called javapath.txt\nOn Windows machines, it must contain the full path to java.exe\nOn linux machines, if Java is properly installed, the path is just [java].'
        javapath=open('javapath.txt','rt')
        jpath = javapath.readline()
        javapath.close()
        if jpath[:4]=='java':
            self.java = jpath
        else:
            self.java = os.path.abspath(jpath)
        self.NamePrefix = []
        self.wdir      = os.getcwd()
        self.javafoil  = os.path.abspath(self.wdir + '/Jfoil/javafoil.jar')
        self.mhclasses = os.path.abspath(self.wdir + '/Jfoil/mhclasses.jar')
        self.tmpdir    = os.path.abspath(self.wdir + '/temp')        
        self.platform  =determinePlatform()
        if self.platform=='lin32':
            self.avl       = self.wdir + '/Avl/avl_lin32'
            self.Xfoil     = self.wdir + '/Xfoil/xfoil_lin32'
        elif self.platform=='lin64':
            self.avl       = self.wdir + '/Avl/avl_lin64'
            self.Xfoil     = self.wdir + '/Xfoil/xfoil_lin64'
        elif self.platform=='nt':
            self.avl       = self.wdir + '\\Avl\\avl_win32.exe'
            self.Xfoil     = self.wdir + '\\Xfoil\\xfoil_win32.exe'
        else:
            print 'Error: unrecognized platform' 
        self.setRandPrefix()
    def setRandPrefix(self,N=7):
        from random import randrange
        prefix = ''
        for ii in range(N):
            r = randrange(0,3)
            if r == 0:
                prefix += str(randrange(0,10))
            elif r == 1:
                prefix += chr(randrange(65,91))
            else:
                prefix += chr(randrange(97,123))
        self.NamePrefix = prefix
    def getTmpFile(self,fileExt = "", addSymbol = ""):
        fileName = self.tmpdir + '/' + self.NamePrefix + addSymbol
        if fileExt != "" : fileName  = fileName + "." + fileExt
        return os.path.abspath(fileName)

class Database:
    def __init__(self):
        self.wdir     = os.path.abspath(os.getcwd() + '/database')
        self.aircraft = os.path.abspath(self.wdir + '/aircraft.xls')
        self.airfoil  = os.path.abspath(self.wdir + '/airfoil.xls')
        self.engine   = os.path.abspath(self.wdir + '/piston_engine.xls')
        self.engineExp= os.path.abspath(self.wdir + '/piston_engine_experimental.xls')
        self.prop     = os.path.abspath(self.wdir + '/prop.xls')
        self.propAirfoil = os.path.abspath(self.wdir + '/prop_airfoil.xls')
        self.drag     = os.path.abspath(self.wdir + '/drag_components.xls')
        self.constants= os.path.abspath(self.wdir + '/constants.xls')

def determinePlatform():
    import platform
    out=''
    if os.name=='nt':
        out='nt'
    elif os.name=='posix':
        if platform.architecture()[0]=='32bit':
            out='lin32'
        elif platform.architecture()[0]=='64bit':
            out='lin64'
        else:
            print 'Error: OS not recognized'
    return out
def fixPaths(pth):
    if os.name=='nt':
        pth=pth.replace('/','\\')
    elif os.name=='posix':
        pth=pth.replace('\\','/')
    else:
        print 'Error: OS not recognized'
    return pth
def testFcn():
    a = MyPaths()
    a.setRandPrefix()
    print a.NamePrefix
    print a.avl
    print fixPaths('C:/test/test2\\test3')
    print os.path.abspath('C:/test/test2\\test3')
    print a.getTmpFile()

class CFD_paths:
    def __init__(self, name_prefix = 'tmp_file'):
        self.wdir = os.getcwd()
        self.name_prefix = name_prefix
        self.fluent =  r'C:\Program Files\ANSYS Inc\v130\fluent\ntbin\win64\fluent.exe'
        self.pointwise = r'C:\Program Files (x86)\Pointwise\PointwiseV16.03R2\win64\bin\Pointwise.exe'
        self.tmpdir = self.wdir + '\\temp'
        self.templates = self.wdir + '\\templates'
        self.template_pw = self.templates + '\\meshJournalTemplate.glf'
        self.template_pw2 = self.templates + '\\CmeshJournalTemplate.glf'
        self.template_pw3 = self.templates + '\\OmeshJournalTemplate.glf'
        self.template_fl = self.templates + '\\fluentJournalTemplate.jou'
        self.set_name()
        
    def _work_file(self, file_ext = "", add_symbol = ""):
        fileName = self.tmpdir + '\\' + self.name_prefix + add_symbol
        if file_ext != "" : fileName  = fileName + "." + file_ext
        return fileName

    def set_name(self,name_prefix = ""):
        if name_prefix != "": self.name_prefix = name_prefix
        self.file_igs = self._work_file('igs')
        self.file_glf = self._work_file('glf')
        self.file_jou = self._work_file('jou')
        self.file_cas = self._work_file('cas')
        self.file_pol = self._work_file('pol')
        self.file_airfoil = self._work_file('dat','coord')

    def set_name_alpha(self,alpha):
        alpha = '_%.0fdeg'%float(alpha)
        self.file_cl_hist = self._work_file('cl', alpha)
        self.file_cd_hist = self._work_file('cd', alpha)
        self.file_cm_hist = self._work_file('cm', alpha)
        self.file_cp_dist = self._work_file('cp', alpha)
        self.list_hist_files = [self.file_cl_hist, self.file_cd_hist, self.file_cm_hist]
    
    def clean(self):
        os.remove(self.file_cas)
        #os.remove(self.file_igs)
        os.remove(self.file_glf)
        os.remove(self.file_jou)
        os.remove(self.file_cd_hist)
        os.remove(self.file_cl_hist)
        os.remove(self.file_cm_hist)

if __name__=="__main__":
    testFcn()