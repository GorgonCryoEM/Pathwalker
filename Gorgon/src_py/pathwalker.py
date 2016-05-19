# Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
# Author:        Ross A. Coleman (racolema@bcm.edu)
# Description:   Widget that lets the user see what chain models are in memory and select the one to work with


import sys
import os
import subprocess
from PyQt4 import QtCore, QtGui
from ui_dialog_pathwalker import Ui_DialogPathwalker
from seq_model.Chain import Chain
from base_dock_widget import BaseDockWidget
from base_viewer import *
from libpyGORGON import VolumeRenderer
from volume_viewer import *
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from seq_model.Chain import Chain

from string import split, upper

class Pathwalker(BaseDockWidget):

  '''
  def __init__(self, main, parent=None):
    super(tabdemo, self).__init__(parent)
    self.tab = QWidget()
    self.tab1 = QWidget()
    self.tab_2 = QWidget()
    self.addTab(self.tab,"Tab 1")
    self.addTab(self.tab1,"Tab 2")
    self.addTab(self.tab_2,"Tab 3")
    self.tab1UI()
    self.tab2UI()
    self.tab3UI()
    self.setWindowTitle("tab demo")

    '''
  def __init__(self, main, parent=None):
    BaseDockWidget.__init__(self, main,"&Pathwalker","Find path with pathwalker","perform_pathwalker","actions-calpha-pathwalker","actions-calpha", QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea | QtCore.Qt.BottomDockWidgetArea, QtCore.Qt.RightDockWidgetArea,parent)
    self.app = main
    self.renderer = VolumeRenderer()
    self.volumeName = "volume name"
    #self.mapname = "map name"
    #self.viewer = VolumeViewer(parent, main)
    self.createUi()
    self.connect(self.ui.preprocessPushButton, QtCore.SIGNAL("clicked()"), self.preprocessButtonPress)
    self.connect(self.ui.pushButtonBrowseAtomScore, QtCore.SIGNAL("clicked (bool)"), self.loadVolume)
    self.connect(self.ui.pushButton_2, QtCore.SIGNAL("clicked()"), self.generateAtomsButtonPress)
    self.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.pathwalkButtonPress)
    #self.connect(self.ui.pushButtonMap, QtCore.SIGNAL("clicked (bool)"), self.loadMap)

    
        #self.connect(self.ui.refreshPushButton, QtCore.SIGNAL("clicked()"), self.refresh)
        #self.connect(self.ui.chainModelsListWidget, QtCore.SIGNAL("currentTextChanged (const QString&)"), self.modelHighlighted)
    

  def preprocessButtonPress(self):
      print 'Preprocessing..'
      self.preprocess()
      print 'Finished Preprocessing'

  def generateAtomsButtonPress(self):
      print 'Generating Pseudo-Atoms...'
      self.findPseudoAtoms()
      print 'Finished Generating Pseudo-Atoms'

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
      subprocess.call(['python','EMAN2/bin/e2pathwalker.py','pseudoatoms.pdb', '--mapfile=EMAN2/bin/map.mrc','--output=path0.pdb','--solver=lkh','--overwrite',dmin,dmax,threshold,mapweight])
      self.generateAtoms("path0.pdb")

  def generateAtoms(self, fileName):
      def setupChain(mychain):            
            self.app.viewers['calpha'].main_chain = mychain
            self.app.viewers['calpha'].loadedChains.append(mychain)
            mychain.setViewer(self.app.viewers['calpha'])
            #Chain.setSelectedChainKey(mychain.getIDs())
            mychain.addCalphaBonds()
            mychain.addSideChainBonds()
            renderer = self.app.viewers['calpha'].renderer
            for i in mychain.residueRange():
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
                  mychainKeys = Chain.loadAllChains(str(self.atomFileName), qparent=self.app)
                  for chainKey in mychainKeys:
                    setupChain(Chain.getChain(chainKey))
              else:
                  mychain = Chain.load(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  setupChain(mychain)
        
              if not self.app.viewers['calpha'].loaded:
                  self.app.viewers['calpha'].dirty = False
                  self.app.viewers['calpha'].loaded = True
                  self.app.viewers['calpha'].setAtomColorsAndVisibility(self.app.viewers['calpha'].displayStyle)                        
                  self.app.viewers['calpha'].emitModelLoadedPreDraw()
                  self.app.viewers['calpha'].emitModelLoaded()
                  self.app.viewers['calpha'].emitViewerSetCenter()    
           
  def createUi(self):
      self.ui = Ui_DialogPathwalker()
      self.ui.setupUi(self)
      #self.ui.sequenceTextEdit.setReadOnly(True)        
      self.ui.preprocessPushButton.setEnabled(True)

  '''
  def loadMap(self, temp):
    self.mapName = str(QtGui.QFileDialog.getOpenFileName(self.app.viewers["volume"], self.app.viewers["volume"].tr("Open Volume"), "", self.app.viewers["volume"].tr(self.app.viewers["volume"].renderer.getSupportedLoadFileFormats())))
    self.ui.lineEditMap.setText(self.mapName)
    #self.app.viewers["volume"].loadDataFromFile(fileName)
        #self.app.actions.getAction("load_Volume").trigger()
    self.bringToFront()
  '''

  def loadVolume(self, temp):
    self.volumeName = str(QtGui.QFileDialog.getOpenFileName(self.app.viewers["volume"], self.app.viewers["volume"].tr("Open Volume"), "", self.app.viewers["volume"].tr(self.app.viewers["volume"].renderer.getSupportedLoadFileFormats())))
    self.ui.lineEditAtomScore.setText(self.volumeName)
    #self.app.viewers["volume"].loadDataFromFile(fileName)
        #self.app.actions.getAction("load_Volume").trigger()
    self.bringToFront()
  '''
  def loadData(self):
      fileName = str(QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Volume"), "", self.tr("MRC Files (*.mrc)")))
      self.loadDataFromFile(fileName)

  def loadDataFromFile(self, fileName):
        self.fileName = fileName
                
        if not self.fileName=="":  
            self.setCursor(QtCore.Qt.WaitCursor)
            
            tokens = split(str(self.fileName), '.')            
            extension = upper(tokens[len(tokens)-1])
            if(extension == "RAW"):
                if(self.rawLoader.exec_() == QtGui.QDialog.Accepted):                
                    self.renderer.loadFileRAW(str(self.fileName), self.rawLoader.bitsPerCell(), self.rawLoader.sizeX(), self.rawLoader.sizeY(), self.rawLoader.sizeZ())
                else:
                    return;
                    
            else:
                self.renderer.loadFile(str(self.fileName))
            self.app.viewers["volume"].setScaleNoEmit(self.renderer.getSpacingX(), self.renderer.getSpacingY(), self.renderer.getSpacingZ())       
            self.loaded = True
            self.dirty = False
            self.setCursor(QtCore.Qt.ArrowCursor)
            self.emitModelLoadedPreDraw()
            self.emitModelLoaded()            
            self.emitViewerSetCenter()
  '''

        
    #def loadWidget(self):
     #   BaseDockWidget.loadWidget(self)
      #  if self.app.actions.getAction("perform_pathwalker").isChecked():
       #     self.refresh()

    #def modelHighlighted(self, chainQString):
     #   currText = str(chainQString)
      #  currChainKey = tuple( currText.split(' - ') )
       # currChain = Chain.getChain(currChainKey)
        #text = str(currChain)
        #self.ui.sequenceTextEdit.setText(text)
        
    #def refresh(self):
     #   self.ui.chainModelsListWidget.clear()
      #  chainKeys = Chain.getChainKeys()
       # for key in chainKeys:
        #    item = str(key[0]) + ' - ' + str(key[1])
         #   self.ui.chainModelsListWidget.addItem(item)

  def preprocess(self):
      print str(self.ui.normalizeLine.text())
      subprocess.call(['python','EMAN2/bin/e2proc3d.py',self.volumeName, 'EMAN2/bin/map.mrc','--process',str(self.ui.normalizeLine.text()),'--process',str(self.ui.zeroThresholdLine.text())])
      self.app.viewers["volume"].loadDataFromFile("EMAN2/bin/map.mrc")
        #os.system('python /users/danzeng/Desktop/Gorgon/Gorgon/src_py/EMAN2/bin/e2proc3d.py /users/danzeng/Desktop/Gorgon/Gorgon/src_py/EMAN2/bin/sub-A.mrc /Users/danzeng/Desktop/Gorgon/Gorgon/src_py/EMAN2/bin/map.mrc --process normalize.edgemean --process threshold.belowtozero')
        #subprocess.call(['python','EMAN2/bin/e2proc3d.py','EMAN2/bin/Untitled.mrc','map.mrc','--process','normalize.edgemean','--process','threshold.belowtozero'])