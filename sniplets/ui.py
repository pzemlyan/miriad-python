# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'untitled1/mainwindow.ui'
#
# Created: Thu Mar 20 17:22:58 2014
#      by: PyQt5 UI code generator 5.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(705, 472)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.centralWidget.sizePolicy().hasHeightForWidth())
        self.centralWidget.setSizePolicy(sizePolicy)
        self.centralWidget.setMaximumSize(QtCore.QSize(2048, 2048))
        self.centralWidget.setObjectName("centralWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralWidget)
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.name = QtWidgets.QLabel(self.centralWidget)
        self.name.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.name.setAlignment(QtCore.Qt.AlignCenter)
        self.name.setObjectName("name")
        self.horizontalLayout_2.addWidget(self.name)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.task = QtWidgets.QComboBox(self.centralWidget)
        self.task.setMinimumSize(QtCore.QSize(200, 0))
        self.task.setObjectName("task")
        self.horizontalLayout_2.addWidget(self.task)
        self.dir = QtWidgets.QPushButton(self.centralWidget)
        self.dir.setAutoRepeatDelay(100)
        self.dir.setObjectName("dir")
        self.horizontalLayout_2.addWidget(self.dir)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.list = QtWidgets.QTableView(self.centralWidget)
        self.list.setGridStyle(QtCore.Qt.DashLine)
        self.list.setSortingEnabled(True)
        self.list.setObjectName("list")
        self.list.horizontalHeader().setStretchLastSection(True)
        self.list.verticalHeader().setStretchLastSection(True)
        self.verticalLayout.addWidget(self.list)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.bexit = QtWidgets.QPushButton(self.centralWidget)
        self.bexit.setObjectName("bexit")
        self.horizontalLayout.addWidget(self.bexit)
        self.save = QtWidgets.QPushButton(self.centralWidget)
        self.save.setObjectName("save")
        self.horizontalLayout.addWidget(self.save)
        self.run = QtWidgets.QPushButton(self.centralWidget)
        self.run.setObjectName("run")
        self.horizontalLayout.addWidget(self.run)
        self.horizontalLayout.setStretch(0, 1)
        self.horizontalLayout.setStretch(1, 1)
        self.horizontalLayout.setStretch(2, 1)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.setStretch(1, 1)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 705, 29))
        self.menuBar.setObjectName("menuBar")
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        self.mainToolBar.setObjectName("mainToolBar")
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.name.setText(_translate("MainWindow", "TaskName"))
        self.dir.setText(_translate("MainWindow", "..."))
        self.bexit.setText(_translate("MainWindow", "Exit"))
        self.save.setText(_translate("MainWindow", "Save"))
        self.run.setText(_translate("MainWindow", "Run"))

