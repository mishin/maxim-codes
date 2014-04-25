# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 12:45:34 2014

@author: Maxim
"""

class DragItem:
    """
    Drag item is similar with mass component. Describes aerodynamic drag of 
    single drag component. It can be obtained from analysis (wing, tail) or 
    table lookup (antenna, cowling)
    """
    n = 0
    def __init__(self):
        DragItem.n += 1
        self.name        = 'drag_component_%d'%DragItem.n
        self.quantity    = 0
        self.frontArea   = 0.0
        self.CD          = 0.0
        self.CDtimesArea = 0.0
        self.CDtotal     = 0.0
    def display(self,header=True):
        """
        Displays short information about current drag item.

        Parameters
        ----------
        header : bool
            if header is False then name and drag coefficient are displayed 
            only. False option is useful for displaying multiple components.
        """
        if header:
            print '{0:26}|{1:1}Drag coef.{1:1}|'.format('Name',' ')
            print '-'*40
        print '{0:27} {1:6.4e}'.format(self.name,self.CDtotal)
    def display_full(self):
        """
        Displays full information about current drag item including its name, 
        drag coefficient, frontal area and quantity
        """
        print '-'*15
        print 'name \t\t: %s'%self.name
        print 'drag coefficient\t: %.2e'%self.CD
        if self.quantity>1:
            print 'quantity\t\t: %d'%self.quantity
        if self.frontArea != 0.0:
            print 'frontal area\t: %.4f'%self.frontArea
            print 'Cd times area\t: %.2e'%self.CDtimesArea

class DragList:
    n = 0
    """
    List of drag components (DragItem)
    
    Parameters
    ----------
    
    name : string
        name of new drag list. Used only for report generations
    """
    def __init__(self,name=''):
        DragList.n += 1
        if name=='':
            name = 'Drag list #%d'%DragList.n
        self.name = name
        self.items = list()
        self.numItems = 0
    def __getitem__(self,i):
        return self.items[i]
    def add_from_db(self,dbEntry):
        """
        Append drag item to the list by parameters

        Parameters
        ----------
        
        dbEntry : list
            List extracted from aircraft db in format [name,quantity,
            frontal area]
      
        Examples
        --------
        
        >>> dbEntry = ['landing gear', 1, 0.085] # line from aircraft.xls
        >>> dragList = DragList()
        >>> dragList.add_from_db(dbEntry)
        """
        newItem           = DragItem()
        newItem.name      = dbEntry[0]
        newItem.quantity  = int(dbEntry[1])
        newItem.frontArea = float(dbEntry[2])
        self.add_item(newItem)
        
    def add_item(self,dragItem):
        """
        Append drag item to the list

        Parameters
        ----------
        
        dragItem : DragItem
        """
        self.items.append(dragItem)
        self.numItems += 1
    def add_list(self,dragList):
        """
        Append all drag items from given list into current one.
        
        Parameters
        ----------
        
        dragList : DragList
            same object type as current class
        """
        for item in dragList:
            self.add_item(item)
    def remove_item(self,name):
        """
        Removes item specified by name from drag list. If item name does not 
        exist, displays message that item was not found without ending 
        execution.
        
        Parameters
        ----------
        
        name : string
            name of the drag item to delete
        """
        try:
            for i,item in enumerate(self.items):
                if item.name==name:
                    del self.items[i]
                    break
            print '%s\tremoved'%name
        except:
            print '%s is not found'%name
    def get_total_drag(self):
        """
        Returns total drag of all items in the list.

        .. math:
            C_{D_{total}} = \sum{C_D*quantity}
        """
        CDtotal = 0.0
        for item in self.items:
            CDtotal += item.CD*item.quantity
        return CDtotal
    def display(self):
        """
        Displays report of the drag components breakdown in tabulated form.
        """
        limiter1 = '='*40
        limiter2 = '-'*40
        print limiter1
        print '%s drag components breakdown'%self.name
        print limiter1
        print '{0:26}|{1:1}Drag coef.{1:1}|'.format('Name',' ')
        print limiter2
        for item in self.items:
            item.display(False)
        print limiter2
        print '{0:27} {1:8.4e}'.format('TOTAL',self.get_total_drag())
        print limiter2
    def read_txt(self,path):
        # was used before creating xls db to read components from text file
        fid = open(path,'rt')
        for line in fid:
            if line!='':
                segLine = line.split()
                item = DragItem()
                item.name = segLine[0]
                item.quantity = float(segLine[1])
                item.frontArea = float(segLine[2])
                self.add_item(item)
        fid.close()


class AircraftDrag:
    """ Data structure with parasite drag breakdown and easy access to aircraft 
    drag components. This object is used in aircraft.drag
    """
    def __init__(self):
        self.friction   = DragList('Friction')
#        self.components = DragList('Components')
        self.total      = DragList('Total')
        self.totalDrag  = 0.0
    def update_total(self):
        self.total = DragList('Total')
        self.total.add_list(self.friction)
        self.total.add_list(self.components)
        self.totalDrag = self.total.get_total_drag()
    def __call__(self):
        return self.totalDrag
    def display(self):
        self.update_total()
        self.total.display()