# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pyqtcontinuum.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(780, 520)
        self.button_fit = QtWidgets.QPushButton(Dialog)
        self.button_fit.setGeometry(QtCore.QRect(590, 20, 75, 21))
        self.button_fit.setObjectName("button_fit")
        self.button_draw = QtWidgets.QPushButton(Dialog)
        self.button_draw.setGeometry(QtCore.QRect(675, 20, 75, 21))
        self.button_draw.setObjectName("button_draw")
        self.label_function = QtWidgets.QLabel(Dialog)
        self.label_function.setGeometry(QtCore.QRect(590, 60, 121, 16))
        self.label_function.setObjectName("label_function")
        self.edit_function = QtWidgets.QLineEdit(Dialog)
        self.edit_function.setGeometry(QtCore.QRect(590, 80, 160, 23))
        self.edit_function.setObjectName("edit_function")
        self.edit_knots = QtWidgets.QLineEdit(Dialog)
        self.edit_knots.setGeometry(QtCore.QRect(590, 130, 75, 23))
        self.edit_knots.setObjectName("edit_knots")
        self.label_knots = QtWidgets.QLabel(Dialog)
        self.label_knots.setGeometry(QtCore.QRect(590, 110, 121, 16))
        self.label_knots.setObjectName("label_knots")
        self.edit_niter = QtWidgets.QLineEdit(Dialog)
        self.edit_niter.setGeometry(QtCore.QRect(590, 230, 75, 23))
        self.edit_niter.setObjectName("edit_niter")
        self.label_niter = QtWidgets.QLabel(Dialog)
        self.label_niter.setGeometry(QtCore.QRect(590, 210, 121, 16))
        self.label_niter.setObjectName("label_niter")
        self.edit_lowrej = QtWidgets.QLineEdit(Dialog)
        self.edit_lowrej.setGeometry(QtCore.QRect(590, 180, 75, 23))
        self.edit_lowrej.setObjectName("edit_lowrej")
        self.label_lowrej = QtWidgets.QLabel(Dialog)
        self.label_lowrej.setGeometry(QtCore.QRect(590, 160, 60, 16))
        self.label_lowrej.setObjectName("label_lowrej")
        self.edit_highrej = QtWidgets.QLineEdit(Dialog)
        self.edit_highrej.setGeometry(QtCore.QRect(675, 180, 75, 23))
        self.edit_highrej.setObjectName("edit_highrej")
        self.label_highrej = QtWidgets.QLabel(Dialog)
        self.label_highrej.setGeometry(QtCore.QRect(675, 160, 60, 16))
        self.label_highrej.setObjectName("label_highrej")
        self.edit_nave = QtWidgets.QLineEdit(Dialog)
        self.edit_nave.setGeometry(QtCore.QRect(675, 230, 75, 23))
        self.edit_nave.setObjectName("edit_nave")
        self.label_nave = QtWidgets.QLabel(Dialog)
        self.label_nave.setGeometry(QtCore.QRect(675, 210, 80, 16))
        self.label_nave.setObjectName("label_nave")
        self.label_samples = QtWidgets.QLabel(Dialog)
        self.label_samples.setGeometry(QtCore.QRect(590, 360, 121, 16))
        self.label_samples.setObjectName("label_samples")
        self.label_grow = QtWidgets.QLabel(Dialog)
        self.label_grow.setGeometry(QtCore.QRect(590, 260, 70, 16))
        self.label_grow.setObjectName("label_grow")
        self.edit_grow = QtWidgets.QLineEdit(Dialog)
        self.edit_grow.setGeometry(QtCore.QRect(590, 280, 75, 23))
        self.edit_grow.setMouseTracking(True)
        self.edit_grow.setObjectName("edit_grow")
        self.edit_samples = QtWidgets.QPlainTextEdit(Dialog)
        self.edit_samples.setGeometry(QtCore.QRect(590, 380, 160, 120))
        self.edit_samples.setObjectName("edit_samples")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.button_fit.setText(_translate("Dialog", "Fit"))
        self.button_draw.setText(_translate("Dialog", "Redraw"))
        self.label_function.setText(_translate("Dialog", "function"))
        self.label_knots.setText(_translate("Dialog", "dwvl_knots"))
        self.label_niter.setText(_translate("Dialog", "n_iterate"))
        self.label_lowrej.setText(_translate("Dialog", "low_rej"))
        self.label_highrej.setText(_translate("Dialog", "high_rej"))
        self.label_nave.setText(_translate("Dialog", "n_average"))
        self.label_samples.setText(_translate("Dialog", "samples"))
        self.label_grow.setText(_translate("Dialog", "grow"))

