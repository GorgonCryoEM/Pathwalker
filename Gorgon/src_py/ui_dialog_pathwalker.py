# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_dialog_pathwalker.ui'
#
# Created by: PyQt4 UI code generator 4.12.dev1605051544
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_DialogPathwalker(object):
    def setupUi(self, DialogModelVisualization):

        DialogModelVisualization.setObjectName(_fromUtf8("DialogModelVisualization"))
        DialogModelVisualization.resize(304, 338)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DialogModelVisualization.sizePolicy().hasHeightForWidth())
        DialogModelVisualization.setSizePolicy(sizePolicy)
        self.gridlayout = QtGui.QGridLayout(DialogModelVisualization)
        self.gridlayout.setObjectName(_fromUtf8("gridlayout"))
        self.tabWidget = QtGui.QTabWidget(DialogModelVisualization)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.gridLayout = QtGui.QGridLayout(self.tab)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))

        spacerItem = QtGui.QSpacerItem(20, 16, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 5, 1, 1, 1)
        self.label_4 = QtGui.QLabel(self.tab)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)

        self.preprocessPushButton = QtGui.QPushButton(self.tab)
        self.preprocessPushButton.setObjectName(_fromUtf8("preprocessPushButton"))
        self.gridLayout.addWidget(self.preprocessPushButton, 5, 0, 1, 1)

        self.label_3 = QtGui.QLabel(self.tab)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)


        self.lineEditAtomScore = QtGui.QLineEdit(self.tab)
        self.lineEditAtomScore.setObjectName(_fromUtf8("lineEditAtomScore"))
        self.gridLayout.addWidget(self.lineEditAtomScore, 0, 0, 1, 1)
        self.pushButtonBrowseAtomScore = QtGui.QPushButton(self.tab)
        self.pushButtonBrowseAtomScore.setObjectName(_fromUtf8("pushButtonBrowseAtomScore"))
        self.gridLayout.addWidget(self.pushButtonBrowseAtomScore, 0, 1, 1, 1)
        

        self.normalizeLine = QtGui.QLineEdit(self.tab)
        self.normalizeLine.setObjectName(_fromUtf8("normalizeLine"))
        self.gridLayout.addWidget(self.normalizeLine, 4, 0, 1, 2)
        self.zeroThresholdLine = QtGui.QLineEdit(self.tab)
        self.zeroThresholdLine.setObjectName(_fromUtf8("zeroThresholdLine"))
        self.gridLayout.addWidget(self.zeroThresholdLine, 2, 0, 1, 2)

        self.tabWidget.addTab(self.tab, _fromUtf8(""))

        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.gridlayout1 = QtGui.QGridLayout(self.tab_2)
        self.gridlayout1.setMargin(0)
        self.gridlayout1.setObjectName(_fromUtf8("gridlayout1"))
        self.label_5 = QtGui.QLabel(self.tab_2)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridlayout1.addWidget(self.label_5, 3, 0, 1, 1)
        self.label_9 = QtGui.QLabel(self.tab_2)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridlayout1.addWidget(self.label_9, 8, 0, 1, 1)
        self.label_7 = QtGui.QLabel(self.tab_2)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridlayout1.addWidget(self.label_7, 5, 0, 1, 1)
        self.label_8 = QtGui.QLabel(self.tab_2)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridlayout1.addWidget(self.label_8, 6, 0, 1, 1)
        self.lineEdit_7 = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_7.setObjectName(_fromUtf8("lineEdit_7"))
        self.gridlayout1.addWidget(self.lineEdit_7, 6, 1, 1, 1)
        self.label_6 = QtGui.QLabel(self.tab_2)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridlayout1.addWidget(self.label_6, 4, 0, 1, 1)
        self.lineEdit_8 = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_8.setObjectName(_fromUtf8("lineEdit_8"))
        self.gridlayout1.addWidget(self.lineEdit_8, 8, 1, 1, 1)
        self.lineEdit_3 = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.gridlayout1.addWidget(self.lineEdit_3, 3, 1, 1, 1)
        self.lineEdit_6 = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_6.setObjectName(_fromUtf8("lineEdit_6"))
        self.gridlayout1.addWidget(self.lineEdit_6, 5, 1, 1, 1)
        self.lineEdit_5 = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.gridlayout1.addWidget(self.lineEdit_5, 4, 1, 1, 1)
        self.pushButton_2 = QtGui.QPushButton(self.tab_2)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.gridlayout1.addWidget(self.pushButton_2, 9, 0, 1, 2)
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.gridlayout2 = QtGui.QGridLayout(self.tab_3)
        self.gridlayout2.setMargin(0)
        self.gridlayout2.setObjectName(_fromUtf8("gridlayout2"))
        self.lineEdit_9 = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_9.setObjectName(_fromUtf8("lineEdit_9"))
        self.gridlayout2.addWidget(self.lineEdit_9, 1, 0, 1, 1)
        self.label_10 = QtGui.QLabel(self.tab_3)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridlayout2.addWidget(self.label_10, 4, 0, 1, 1)
        self.label_2 = QtGui.QLabel(self.tab_3)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridlayout2.addWidget(self.label_2, 2, 0, 1, 1)
        self.lineEdit_12 = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_12.setObjectName(_fromUtf8("lineEdit_12"))
        self.gridlayout2.addWidget(self.lineEdit_12, 7, 0, 1, 1)
        self.label_11 = QtGui.QLabel(self.tab_3)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridlayout2.addWidget(self.label_11, 6, 0, 1, 1)
        self.lineEdit_11 = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_11.setObjectName(_fromUtf8("lineEdit_11"))
        self.gridlayout2.addWidget(self.lineEdit_11, 5, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridlayout2.addItem(spacerItem1, 9, 0, 1, 1)
        self.lineEdit_10 = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_10.setObjectName(_fromUtf8("lineEdit_10"))
        self.gridlayout2.addWidget(self.lineEdit_10, 3, 0, 1, 1)
        self.label = QtGui.QLabel(self.tab_3)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridlayout2.addWidget(self.label, 0, 0, 1, 1)
        self.pushButton_4 = QtGui.QPushButton(self.tab_3)
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.gridlayout2.addWidget(self.pushButton_4, 8, 0, 1, 1)
        self.tabWidget.addTab(self.tab_3, _fromUtf8(""))
        self.gridlayout.addWidget(self.tabWidget, 0, 0, 1, 1)

        self.retranslateUi(DialogModelVisualization)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(DialogModelVisualization)

    def retranslateUi(self, DialogModelVisualization):
        DialogModelVisualization.setWindowTitle(QtGui.QApplication.translate("DialogModelVisualization", "Preprocess", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonBrowseAtomScore.setText(QtGui.QApplication.translate("DialogModelVisualization", "Volume...", None, QtGui.QApplication.UnicodeUTF8))
        self.preprocessPushButton.setText(QtGui.QApplication.translate("DialogModelVisualization", "Preprocess", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QtGui.QApplication.translate("DialogModelVisualization", "Preprocess", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(_translate("DialogModelVisualization", "Processor", None))
        self.label_9.setText(_translate("DialogModelVisualization", "Segment Separation", None))
        self.label_7.setText(_translate("DialogModelVisualization", "Amplitude Weight", None))
        self.label_8.setText(_translate("DialogModelVisualization", "Threshold", None))
        self.label_6.setText(_translate("DialogModelVisualization", "Segments", None))
        self.pushButton_2.setText(_translate("DialogModelVisualization", "Generate PseudoAtoms", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("DialogModelVisualization", "Pseudoatom Generation", None))
        self.lineEdit_7.setText(_translate("DialogModelVisualization", "10", None))
        self.lineEdit_8.setText(_translate("DialogModelVisualization", "1", None))
        self.lineEdit_3.setText(_translate("DialogModelVisualization", "segment.kmeans", None))
        self.lineEdit_6.setText(_translate("DialogModelVisualization", "1", None))
        self.lineEdit_5.setText(_translate("DialogModelVisualization", "1022", None))
        self.lineEdit_9.setText(_translate("DialogModelVisualization", "1", None))
        self.label_10.setText(_translate("DialogModelVisualization", "Map Threshold", None))
        self.label_2.setText(_translate("DialogModelVisualization", "Maximum Distance", None))
        self.lineEdit_12.setText(_translate("DialogModelVisualization", "200", None))
        self.label_11.setText(_translate("DialogModelVisualization", "Map Weight", None))
        self.lineEdit_11.setText(_translate("DialogModelVisualization", "12", None))
        self.lineEdit_10.setText(_translate("DialogModelVisualization", "10", None))
        self.label.setText(_translate("DialogModelVisualization", "Minimum Distance", None))
        self.pushButton_4.setText(_translate("DialogModelVisualization", "Pathwalk", None))
        self.normalizeLine.setText(_translate("DialogModelVisualization", "normalize.edgemean", None))
        self.zeroThresholdLine.setText(_translate("DialogModelVisualization", "threshold.belowtozero", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), QtGui.QApplication.translate("DialogModelVisualization", "Pathwalk", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(_translate("DialogModelVisualization", "Processor 2", None))
        self.label_3.setText(_translate("DialogModelVisualization", "Processor 1", None))


from colored_push_button import ColoredPushButton

