# -*- coding: utf-8 -*-
"""
Created on Tue Jan 07 14:37:02 2014

Module for manipulation with standard files and paths used in software

TODO:
extend with CFD, CAD and other paths if necessary

@author: Maxim
"""
import os
from random import randrange
import platform
from warnings import warn

class dbPaths:
    """ structure with database paths"""
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


class MyPaths:
    def __init__(self):
        self._init_javapath()
        self.db   = dbPaths()
        self.wdir = os.getcwd()
        self.javafoil  = os.path.abspath(self.wdir + '/Jfoil/javafoil.jar')
        self.mhclasses = os.path.abspath(self.wdir + '/Jfoil/mhclasses.jar')
        self.tmpDir    = os.path.abspath(self.wdir + '/temp')
        self.platform  = determine_platform()
    
    def set_file_prefix(self,prefix):
        """
        sets temporal filename prefix that will be used with get_tmp_file
        """
        self.namePrefix = prefix

    def set_file_prefix_random(self,n=4):
        """
        sets random filename prefix that will be used with get_tmp_file
        """
        prefix = ''
        for i in range(n):
            r = randrange(0,3)
            if r == 0:
                prefix += str(randrange(0,10))
            elif r == 1:
                prefix += chr(randrange(65,91))
            else:
                prefix += chr(randrange(97,123))
        self.namePrefix = prefix
    
    def get_tmp_file(self,fileExt=None,addSymbol=''):
        """
        Returns temporal filename with given file extension and file prefix. 
        Absolute path is returned.
        """
        fileName = self.tmpDir + '/' + self.namePrefix + addSymbol
        if not fileExt==None:
            fileName += '.'+fileExt
        return os.path.abspath(fileName)

    def _init_javapath(self):
        try:
            with open('javapath.txt','rt') as javapath: pass
        except IOError:
            print 'Error: The root folder \n['+os.getcwd()+']'
            print 'must contain a file called javapath.txt'
            print 'On Windows, it must contain the full path to java.exe'
            print 'On linux, if Java is properly installed, the path is just [java].'
        javapath=open('javapath.txt','rt')
        jpath = javapath.readline()
        javapath.close()
        if jpath[:4]=='java':
            self.java = jpath
        else:
            self.java = os.path.abspath(jpath)


def determine_platform():
    if os.name=='nt':
        out='nt'
    elif os.name=='posix':
        if platform.architecture()[0]=='32bit':
            out='lin32'
        elif platform.architecture()[0]=='64bit':
            out='lin64'
    else:
        warn('OS not recognized. Using default value: WinNT')
        out = 'nt'
    return out


# --- debug section ---
def run_test1():
    pth = MyPaths()
    pth.set_file_prefix_random()
    print pth.namePrefix
    print pth.get_tmp_file('txt')
    pth.set_file_prefix_random()
    print pth.get_tmp_file('txt','pol')
    pth.set_file_prefix('my_file')
    print pth.get_tmp_file('psd')
    print pth.db.airfoil
    print pth.platform

if __name__=="__main__":
    run_test1()