# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'untitled1/listitem.ui'
#
# Created: Tue Mar 18 15:06:33 2014
#      by: PyQt5 UI code generator 5.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ListItem(object):
    def setupUi(self, ListItem):
        ListItem.setObjectName("ListItem")
        ListItem.resize(400, 69)
        self.horizontalLayoutWidget = QtWidgets.QWidget(ListItem)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 401, 71))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.name_box = QtWidgets.QLabel(self.horizontalLayoutWidget)
        self.name_box.setContentsMargins(10, 10, 10, 10)
        self.name_box.setObjectName("name_box")
        self.horizontalLayout.addWidget(self.name_box)
        self.edit_box = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.edit_box.setObjectName("edit_box")
        self.horizontalLayout.addWidget(self.edit_box)
        self.path_btn = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.path_btn.sizePolicy().hasHeightForWidth())
        self.path_btn.setSizePolicy(sizePolicy)
        self.path_btn.setObjectName("path_btn")
        self.horizontalLayout.addWidget(self.path_btn)

        self.retranslateUi(ListItem)
        QtCore.QMetaObject.connectSlotsByName(ListItem)

    def retranslateUi(self, ListItem):
        _translate = QtCore.QCoreApplication.translate
        ListItem.setWindowTitle(_translate("ListItem", "Form"))
        self.name_box.setText(_translate("ListItem", "Name"))
        self.path_btn.setText(_translate("ListItem", "Path"))

