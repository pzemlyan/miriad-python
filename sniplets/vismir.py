# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 18:15:48 2014

@author: pete
"""

from scipy.optimize import leastsq
import numpy as np
import miriad
import mirtask
from mirexec import *

def get_name(task, name=None):
    if ((name == None) and (task._name!=None)):
        name = task._name+".def"
    if name == None:
        name = task.__class__.__name__[4:].lower()+".def"
    return name
    

def parse_file(task,name=None):
    if name!=None:
        name = get_name(task)        
    dic = {}
    for i in task._keywords:
        dic[i]=''

    print(name+' in '+os.getcwd())
    try :
        data = open(name).read()
    except (TypeError, FileNotFoundError):
        print('not found')
        data = ''
    data = data.split('\n')
    for s in data:
        line = s.split()
        if (np.size(line) < 3): continue
        if (line[0][0]in'0123456789')or(line[0][0:2]=='in'):line[0]=line[0]+'_'
        dic[line[0]]=line[2]
    task.set(**dic)
    return (task,dic)
    
def toArr(dic):
    arr = []
    arr2 = []
    for key,val in dic.items():
        arr.append(key)
        arr2.append(val)
    return np.array([arr,arr2]).T
import sys

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtQuick import QQuickView
from PyQt5.QtQml import qmlRegisterType, QQmlComponent, QQmlEngine
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5 import QtOpenGL
from ui import Ui_MainWindow
from listitem import Ui_ListItem
import sys

        
        
        
class Controller(QObject):
    @pyqtSlot(QObject)
    def thingSelected(self, wrapper):
        print ('User clicked on:', wrapper)
            
    @pyqtSlot(QObject)
    def edit(self, modelData):
        print ('Edit', modelData)
            
            

class EditBoxDelegate(QItemDelegate):
    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        print ('QEditor')
        return editor

    def setEditorData(self, spinBox, index):
        value = index.model().data(index, Qt.EditRole).value()
        spinBox.setText(value)

    def setModelData(self, editBox, model, index):
        value = editBox.text()
        print ('Save')
        model.setData(index, value, Qt.EditRole)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)

class MyTableModel(QAbstractTableModel):
    def __init__(self, parent=None, data=[],*args):
        QAbstractTableModel.__init__(self, parent, *args)
        self._data = data        
        
    def rowCount(self, parent):
        #кол-во строк
        return self._data.shape[0]
    def columnCount(self, parent):
        #кол-во колонок
        return self._data.shape[1]
        
    def flags(self, index):
        flags = super(MyTableModel, self).flags(index)
        if index.column() in (1,2):
            flags |= Qt.ItemIsEditable | Qt.ItemIsEnabled
        return flags
    def data(self, index, role):
        #print(index.row(),index.column(),role)
        ##тут фунция вытягивания данных
        if not index.isValid():
            return QVariant()
        elif role == Qt.EditRole:
            print (str(u"эдит ")+str(index.row())+"-"+str(index.column()))
            return QVariant(str(self._data[index.row(),index.column()]))
        elif role != Qt.DisplayRole:
            return QVariant()
            
        ## для отладки-можно видеть обращения к функции
        #print QString(u"клетка ")+str(index.row())+"-"+str(index.column())
#        print (str(u"клетка ")+str(index.row())+"-"+str(index.column()))
        return QVariant(str(self._data[index.row(),index.column()]))
        
    def setData(self, index, value, role=Qt.EditRole):
        if not index.isValid():
            return
        elif role == Qt.EditRole:
            print (str(u"эдит c")+str(index.row())+"-"+str(index.column()))
            self._data[index.row(),index.column()] = value
            self.dataChanged.emit(index, index)
            return 
        elif role != Qt.DisplayRole:
            return QVariant()
        
    def headerData(self, col, orientation, role):
        ## тут задаются заголовки
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(str(u"Колонка ")+str(col))
        if orientation == Qt.Horizontal and role == Qt.EditRole:
            return QVariant(str(u"Колонка ")+str(col))
        return QVariant()
        
        
    def getItem(self, index):
        if index.isValid():
            item = index.internalPointer()
            print(itme)
            if item:
                return item

        return self.rootItem


        
    def setArr(self,data):
        topLeft = self.createIndex(0,0)
        bottomRight = self.createIndex(data.shape[0],2)
        print(topLeft.row(),bottomRight.row())
        self._data=data
        self.resetInternalData()


class CalculatorForm(QMainWindow):

    def on_item_changed(self, item):
        # If the changed item is not checked, don't bother checking others
        if not item.checkState():
            return
     
        # Loop through the items until you get None, which
        # means you've passed the end of the list
        i = 0
        while self.model.item(i):
            if not self.model.item(i).checkState():
                return
            i += 1
     
        app.quit()

    def __init__(self, arr, task, parent=None):
        QMainWindow.__init__(self, parent)
        self._task = task
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        list = self.ui.list
        self.model = MyTableModel(list,data=arr)
 
        list.setItemDelegate(EditBoxDelegate())
        list.setModel(self.model)
        
        list.clicked.connect(self.table_pressed) 
        
        combo = self.ui.task
        
        self.tasks = TaskBase.__subclasses__()
        for i in self.tasks:
            if (i._name != None):
                combo.addItem(i._name)
            else :
                combo.addItem(i.__name__)
        combo.currentIndexChanged.connect(self.task_changed)
        self.ui.name.setText(task._name)        
        self.ui.dir.clicked.connect(self.dirdial)
        self.ui.save.clicked.connect(self.save)
        
        self.ui.bexit.clicked.connect(lambda y: app.quit())
        self.ui.run.clicked.connect(self.launch)


    def save(self):
        data = self.model._data
        f = open(get_name(self._task),'w')
        for i in range(data.shape[0]):
            f.write(data[i,0]+"\t\t=\t"+data[i,1]+'\n')
        f.close()
        print(get_name(self._task)+" saved!")

    def dirdial(self):
        string = QFileDialog.getExistingDirectory(self, "Open Directory","",                                                
                                                QFileDialog.ShowDirsOnly)
        print('chdir:'+string)
        if string:
            os.chdir(string)
        
    @pyqtSlot(int)
    def task_changed(self, x):
        self._task = self.tasks[x]()
        if (self._task._name != None):
            self.ui.name.setText(self._task._name)
        else :
            self.ui.name.setText(self._task.__class__.__name__)
            
        self._task, dic = parse_file(self._task)
        print(dic)
        self.model.setArr(toArr(dic))

    @pyqtSlot(QModelIndex)
    def table_pressed(self, x):
        print(x.column(),x.row())
        
    def launch(self,x):
        self.save()
        arr = self.model._data
        print(arr)
        dic = {}
        for i in arr:
            dic[i[0]]=i[1]
        self._task.set(**dic)
        self._task.run()

    @pyqtSlot(int)
    def on_exit_pressed(self, value):
        self.ui.list.setText(str(value + self.ui.inputSpinBox2.value()))

    @pyqtSlot(int)
    def on_inputSpinBox2_valueChanged(self, value):
        self.ui.outputWidget.setText(str(value + self.ui.inputSpinBox1.value()))


class ListItem(QStandardItem):

    def on_item_changed(self, item):
        # If the changed item is not checked, don't bother checking others
        if not item.checkState():
            return
     
        # Loop through the items until you get None, which
        # means you've passed the end of the list
        i = 0
        while self.model.item(i):
            if not self.model.item(i).checkState():
                return
            i += 1
     
        app.quit()

    def __init__(self, parent=None, data=('','')):
        QStandardItem.__init__(self, parent)
        (name,text) = data
        self.ui = Ui_ListItem()
        
        self.ui.setupUi(QWidget())
        self.ui.name_box.setText(name)
        self.ui.edit_box.setText(text)
        self.ui.edit_box.textChanged.connect(lambda x: print (x))
        
        
    
        


import os
if __name__ == "__main__":
    app = QApplication(sys.argv)
    os.chdir('../Documents/astro/sio/n_reg/')
    task, dic = parse_file(TaskCgDisp(),name='cgdisp.def')
    calculator = CalculatorForm(task=task,arr=toArr(dic))
    calculator.show()
    app.exec_()

