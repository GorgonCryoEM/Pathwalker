import sys
import os
import subprocess
import csv
import ctypes
from PyQt4 import QtCore, QtGui
from ui_dialog_pathwalker import Ui_DialogPathwalker
from seq_model.Chain import Chain
from base_dock_widget import BaseDockWidget
from base_viewer import *
from libpyGORGON import VolumeRenderer
from libpyGORGON import PDBBond
from volume_viewer import *
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from seq_model.Chain import Chain

from string import split, upper

#todo: if pathwalker mode, allow for calpha_viewer to delete or add using ctrl + click
class Pathwalker(BaseDockWidget):

  def __init__(self, main, parent=None):
    BaseDockWidget.__init__(self, main,"&Pathwalker","Find path with pathwalker","perform_pathwalker","actions-calpha-pathwalker","actions-calpha", QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea | QtCore.Qt.BottomDockWidgetArea, QtCore.Qt.RightDockWidgetArea,parent)
    self.app = main
    self.renderer = VolumeRenderer()
    self.volumeName = "volume name"
    self.nobonds = ""
    self.createUi()
    self.pathWalkerMode = True
    self.connect(self.ui.preprocessPushButton, QtCore.SIGNAL("clicked()"), self.preprocessButtonPress)
    self.connect(self.ui.pushButtonBrowseAtomScore, QtCore.SIGNAL("clicked (bool)"), self.loadVolume)
    self.connect(self.ui.pushButton_2, QtCore.SIGNAL("clicked()"), self.generateAtomsButtonPress)
    self.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.pathwalkButtonPress)
    self.connect(self.ui.pushButton_5, QtCore.SIGNAL("clicked()"), self.addDeletedBond)
    self.connect(self.ui.pushButton_6, QtCore.SIGNAL("clicked()"), self.removeBondFromList)
    self.connect(self.ui.pushButton_7, QtCore.SIGNAL("clicked()"), self.clearDeleteListRows)

    self.connect(self.ui.pushButton_8, QtCore.SIGNAL("clicked()"), self.addNewBond)
    self.connect(self.ui.pushButton_9, QtCore.SIGNAL("clicked()"), self.removeNewBondFromList)
    self.connect(self.ui.pushButton_10, QtCore.SIGNAL("clicked()"), self.clearAddListRows)
    self.connect(self.ui.pushButton_11, QtCore.SIGNAL("clicked()"), self.setCTermini)
    self.connect(self.ui.pushButton_12, QtCore.SIGNAL("clicked()"), self.setNTermini)

    self.connect(self.ui.checkBox, QtCore.SIGNAL("clicked()"), self.viewDeletedBonds)
    self.connect(self.ui.checkBox_2, QtCore.SIGNAL("clicked()"), self.viewAddedBonds)    
    #self.connect(self.ui.pushButton_5, QtCore.SIGNAL("clicked()"), self.deleteBonds)      
    #self.connect(self.ui.pushButton_6, QtCore.SIGNAL("clicked()"), self.createBonds)
    #self.connect(self.ui.pushButton_7, QtCore.SIGNAL("clicked()"), self.setCTermini)
    #self.connect(self.ui.pushButton_8, QtCore.SIGNAL("clicked()"), self.setNTermini)   

  def preprocessButtonPress(self):
      print 'Preprocessing..'
      self.preprocess()
      print 'Finished Preprocessing'

  def generateAtomsButtonPress(self):
      print 'Generating Pseudo-Atoms...'
      self.findPseudoAtoms()
      print 'Finished Generating Pseudo-Atoms'

  def clearDeleteListRows(self):
      self.ui.tableWidget.setRowCount(0)

  def clearAddListRows(self):
      self.ui.tableWidget_2.setRowCount(0)

  def findPseudoAtoms(self):
    paramStrings = "--process="+str(self.ui.lineEdit_3.text())+":ampweight="+str(self.ui.lineEdit_6.text())+":nseg="+str(self.ui.lineEdit_5.text())+":verbose=1:minsegsep="+str(self.ui.lineEdit_8.text())+":pseudoatom=1:thr="+str(self.ui.lineEdit_7.text())
    subprocess.call(['python','EMAN2/bin/e2segment3d.py','EMAN2/bin/map.mrc','--pdbout=pseudoatoms.pdb',paramStrings])
    self.generateAtoms("pseudoatoms.pdb")

  def pathwalkButtonPress(self):
      print 'Pathwalking...'
      self.pathWalk()
      print 'Finished Pathwalking'

  def pathWalk(self):
      dmin = "--dmin="+str(self.ui.lineEdit_9.text())
      dmax = "--dmax="+str(self.ui.lineEdit_10.text())
      threshold = "--mapthresh="+str(self.ui.lineEdit_11.text())
      mapweight = "--mapweight="+str(self.ui.lineEdit_12.text())

      bondsToPrevent = "--nobonds="
      allRows = self.ui.tableWidget.rowCount()
      for row in range(allRows):
        currentRow = str(self.ui.tableWidget.item(row, 0).text()) + " " + str(self.ui.tableWidget.item(row, 1).text()) + " "
        bondsToPrevent += currentRow

      newBonds = "--newbonds="
      allRowsNewBonds = self.ui.tableWidget_2.rowCount()
      for row in range(allRowsNewBonds):
        currentRow = str(self.ui.tableWidget_2.item(row, 0).text()) + " " + str(self.ui.tableWidget_2.item(row, 1).text()) + " "
        newBonds += currentRow

      cTerminus = "--start="+str(self.ui.lineEdit_13.text())
      nTerminus = "--end="+str(self.ui.lineEdit_14.text())
      if cTerminus == "--start=":
        if nTerminus == "--end=":
          subprocess.call(['python','EMAN2/bin/e2pathwalker.py','pseudoatoms.pdb', '--mapfile=EMAN2/bin/map.mrc','--output=path0.pdb','--solver=lkh','--overwrite',dmin,dmax,threshold,mapweight,bondsToPrevent,newBonds])
        else:
          subprocess.call(['python','EMAN2/bin/e2pathwalker.py','pseudoatoms.pdb', '--mapfile=EMAN2/bin/map.mrc','--output=path0.pdb','--solver=lkh','--overwrite',dmin,dmax,threshold,mapweight,nTerminus,bondsToPrevent,newBonds])
      else:
        if nTerminus == "--end=":
          subprocess.call(['python','EMAN2/bin/e2pathwalker.py','pseudoatoms.pdb', '--mapfile=EMAN2/bin/map.mrc','--output=path0.pdb','--solver=lkh','--overwrite',dmin,dmax,threshold,mapweight,cTerminus,bondsToPrevent,newBonds])
        else:
          subprocess.call(['python','EMAN2/bin/e2pathwalker.py','pseudoatoms.pdb', '--mapfile=EMAN2/bin/map.mrc','--output=path0.pdb','--solver=lkh','--overwrite',dmin,dmax,threshold,mapweight,nTerminus,cTerminus,bondsToPrevent,newBonds])
      self.app.viewers['calpha'].unloadData()
      self.generatePathwalkedAtoms("path0.pdb")
      self.app.viewers['calpha'].emitModelChanged()

  def generateAtoms(self, fileName):
      def setupChain(mychain):            
            self.app.viewers['calpha'].main_chain = mychain
            self.app.viewers['calpha'].loadedChains.append(mychain)
            mychain.setViewer(self.app.viewers['calpha'])
            #mychain.addCalphaBondsPathwalker()
            #mychain.addSideChainBonds()
            renderer = self.app.viewers['calpha'].renderer
            for i in mychain.unsortedResidueRange():
                for atomName in mychain[i].getAtomNames():
                    atom = mychain[i].getAtom(atomName)
                    if atom:
                        atom = renderer.addAtom(atom)
                        mychain[i].addAtomObject(atom)
      
      self.atomFileName = fileName
      fileNameTemp = self.atomFileName
      self.app.viewers['calpha'].whichChainID = None
      filename = unicode(self.atomFileName)
      if filename.split('.')[-1].lower() == 'pdb':
          dlg = CAlphaChooseChainToLoadForm(unicode(self.atomFileName))
          if dlg.exec_():
              self.app.viewers['calpha'].whichChainID = dlg.whichChainID
              #if not self.atomFileName.isEmpty():
              if(self.app.viewers['calpha'].loaded):
                  self.app.viewers['calpha'].unloadData()
              
              self.atomFileName = fileNameTemp
                    
              if self.app.viewers['calpha'].whichChainID == 'ALL':
                  mychainKeys = Chain.loadAllChainsPathwalker(str(self.atomFileName), qparent=self.app)
                  for chainKey in mychainKeys:
                    setupChain(Chain.getChain(chainKey))
              else:
                  mychain = Chain.__loadFromPDBPathwalker(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  #mychain = Chain.load(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  setupChain(mychain)
        
              if not self.app.viewers['calpha'].loaded:
                  #self.app.viewers['calpha'].setDisplayStyle(6)
                  #self.app.viewers['calpha'].displayStyle = 6
                  #self.app.viewers['calpha'].renderer.setDisplayStyle(6)
                  #self.app.viewers['calpha'].setAtomColorsAndVisibility(6)
                  #self.app.viewers['calpha'].modelChangedPathwalker()
                  self.app.viewers['calpha'].dirty = False
                  self.app.viewers['calpha'].loaded = True
                  self.app.viewers['calpha'].setAtomColorsAndVisibility(self.app.viewers['calpha'].displayStyle)                        
                  self.app.viewers['calpha'].emitModelLoadedPreDraw()
                  #self.app.viewers['calpha'].emitModelPathwalker()
                  self.app.viewers['calpha'].emitModelLoaded()
                  self.app.viewers['calpha'].emitViewerSetCenter()

  def printDeleted(self):
       self.app.viewers['calpha'].main_chain.printDeletedBonds()

  def generatePathwalkedAtoms(self, fileName):
      def setupChain(mychain):            
            self.app.viewers['calpha'].main_chain = mychain
            self.app.viewers['calpha'].loadedChains.append(mychain)
            mychain.setViewer(self.app.viewers['calpha'])
            mychain.addCalphaBondsPathwalker()
            #mychain.addSideChainBonds()
            renderer = self.app.viewers['calpha'].renderer
            for i in mychain.unsortedResidueRange():
                for atomName in mychain[i].getAtomNames():
                    atom = mychain[i].getAtom(atomName)
                    if atom:
                        atom = renderer.addAtom(atom)
                        mychain[i].addAtomObject(atom)
      
      self.atomFileName = fileName
      fileNameTemp = self.atomFileName
      self.app.viewers['calpha'].whichChainID = None
      filename = unicode(self.atomFileName)
      if filename.split('.')[-1].lower() == 'pdb':
          dlg = CAlphaChooseChainToLoadForm(unicode(self.atomFileName))
          if dlg.exec_():
              self.app.viewers['calpha'].whichChainID = dlg.whichChainID
              #if not self.atomFileName.isEmpty():
              if(self.app.viewers['calpha'].loaded):
                  self.app.viewers['calpha'].unloadData()
              
              self.atomFileName = fileNameTemp
                    
              if self.app.viewers['calpha'].whichChainID == 'ALL':
                  mychainKeys = Chain.loadAllChainsPathwalker(str(self.atomFileName), qparent=self.app)
                  for chainKey in mychainKeys:
                    setupChain(Chain.getChain(chainKey))
              else:
                  mychain = Chain.__loadFromPDBPathwalker(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  setupChain(mychain)
        
              if not self.app.viewers['calpha'].loaded:
                  self.app.viewers['calpha'].dirty = False
                  self.app.viewers['calpha'].loaded = True
                  self.app.viewers['calpha'].setAtomColorsAndVisibility(self.app.viewers['calpha'].displayStyle)                        
                  self.app.viewers['calpha'].emitModelLoadedPreDraw()
                  self.app.viewers['calpha'].emitModelLoaded()
                  self.app.viewers['calpha'].emitViewerSetCenter()
           

  def viewDeletedBonds(self):
      allRows = self.ui.tableWidget.rowCount()
      
      if self.ui.checkBox.checkState() == QtCore.Qt.Checked:
        for row in range(allRows):
          self.app.viewers['calpha'].renderer.drawDeletedBond(int(self.ui.tableWidget.item(row, 0).text()), int(self.ui.tableWidget.item(row, 1).text()))
        self.app.viewers['calpha'].emitModelChanged()
      else:
        for row in range(allRows):
          self.app.viewers['calpha'].renderer.undrawDeletedBond(int(self.ui.tableWidget.item(row, 0).text()), int(self.ui.tableWidget.item(row, 1).text()))
        self.app.viewers['calpha'].emitModelChanged()

  def viewAddedBonds(self):
      allRows = self.ui.tableWidget_2.rowCount()
      
      if self.ui.checkBox_2.checkState() == QtCore.Qt.Checked:
        for row in range(allRows):
          self.app.viewers['calpha'].renderer.drawAddedBond(int(self.ui.tableWidget_2.item(row, 0).text()), int(self.ui.tableWidget_2.item(row, 1).text()))
        self.app.viewers['calpha'].emitModelChanged()
      else:
        for row in range(allRows):
          self.app.viewers['calpha'].renderer.undrawAddedBond(int(self.ui.tableWidget_2.item(row, 0).text()), int(self.ui.tableWidget_2.item(row, 1).text()))
        self.app.viewers['calpha'].emitModelChanged()

  def removeBondFromList(self):
      self.ui.tableWidget.removeRow(self.ui.tableWidget.currentRow())

  def removeNewBondFromList(self):
      self.ui.tableWidget_2.removeRow(self.ui.tableWidget_2.currentRow())

  def setCTermini(self):
      selectedAtoms = self.app.viewers['calpha'].main_chain.getSelectionSerial()
      selectedAtoms = str(selectedAtoms)
      selectedAtoms = selectedAtoms.translate(None, '![@#]$,')
      selectedAtoms = selectedAtoms.split()
      if len(selectedAtoms) != 1:
        return
      self.ui.lineEdit_13.setText(str(selectedAtoms[0]))

  def setNTermini(self):
      selectedAtoms = self.app.viewers['calpha'].main_chain.getSelectionSerial()
      selectedAtoms = str(selectedAtoms)
      selectedAtoms = selectedAtoms.translate(None, '![@#]$,')
      selectedAtoms = selectedAtoms.split()
      if len(selectedAtoms) != 1:
        return
      self.ui.lineEdit_14.setText(str(selectedAtoms[0]))

  def addDeletedBond(self):
      selectedAtoms = self.app.viewers['calpha'].main_chain.getSelectionSerial()
      selectedAtoms = str(selectedAtoms)
      selectedAtoms = selectedAtoms.translate(None, '![@#]$,')
      selectedAtoms = selectedAtoms.split()
      if len(selectedAtoms) != 2:
        return
      rowPosition = self.ui.tableWidget.rowCount()
      self.ui.tableWidget.insertRow(rowPosition)
      self.ui.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem(selectedAtoms[0]))
      self.ui.tableWidget.setItem(rowPosition , 1, QtGui.QTableWidgetItem(selectedAtoms[1]))
      currentDistance = self.app.viewers['calpha'].renderer.findDistance(int(selectedAtoms[0]), int(selectedAtoms[1]))
      self.ui.tableWidget.setItem(rowPosition , 2, QtGui.QTableWidgetItem(currentDistance))

  def addNewBond(self):
      selectedAtoms = self.app.viewers['calpha'].main_chain.getSelectionSerial()
      selectedAtoms = str(selectedAtoms)
      selectedAtoms = selectedAtoms.translate(None, '![@#]$,')
      selectedAtoms = selectedAtoms.split()
      if len(selectedAtoms) != 2:
        return
      rowPosition = self.ui.tableWidget_2.rowCount()
      self.ui.tableWidget_2.insertRow(rowPosition)
      self.ui.tableWidget_2.setItem(rowPosition , 0, QtGui.QTableWidgetItem(selectedAtoms[0]))
      self.ui.tableWidget_2.setItem(rowPosition , 1, QtGui.QTableWidgetItem(selectedAtoms[1]))
      currentDistance = self.app.viewers['calpha'].renderer.findDistance(int(selectedAtoms[0]), int(selectedAtoms[1]))
      self.ui.tableWidget_2.setItem(rowPosition , 2, QtGui.QTableWidgetItem(currentDistance))

  def createBonds(self):
      selectedCreated = self.app.viewers['calpha'].main_chain.getSelectionSerial()
      selectedCreated = str(selectedCreated)
      selectedCreated = selectedCreated.translate(None, '![@#]$,')
      #self.ui.lineEdit_14.setText(str(selectedCreated))
      self.newBonds = selectedCreated
      if self.newBonds:
        self.app.viewers['calpha'].renderer.addSelectedBonds(self.newBonds)
        self.app.viewers['calpha'].emitModelChanged() 

  def createUi(self):
      self.ui = Ui_DialogPathwalker()
      self.ui.setupUi(self)      
      self.ui.preprocessPushButton.setEnabled(True)


  def loadVolume(self, temp):
    self.volumeName = str(QtGui.QFileDialog.getOpenFileName(self.app.viewers["volume"], self.app.viewers["volume"].tr("Open Volume"), "", self.app.viewers["volume"].tr(self.app.viewers["volume"].renderer.getSupportedLoadFileFormats())))
    self.ui.lineEditAtomScore.setText(self.volumeName)
    self.bringToFront()
 

  def preprocess(self):
      subprocess.call(['python','EMAN2/bin/e2proc3d.py',self.volumeName, 'EMAN2/bin/map.mrc','--process',str(self.ui.normalizeLine.text()),'--process',str(self.ui.zeroThresholdLine.text())])
      self.app.viewers["volume"].loadDataFromFile("EMAN2/bin/map.mrc")


        