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

#todo: if pathwalker mode, allow for calpha_viewer to delete or add using ctrl + click
class Pathwalker(BaseDockWidget):

  def __init__(self, main, parent=None):
    BaseDockWidget.__init__(self, main,"&Pathwalker","Find path with pathwalker","perform_pathwalker","actions-calpha-pathwalker","actions-calpha", QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea | QtCore.Qt.BottomDockWidgetArea, QtCore.Qt.RightDockWidgetArea,parent)
    self.app = main
    self.renderer = VolumeRenderer()
    self.volumeName = "volume name"
    self.createUi()
    self.pathWalkerMode = True
    self.connect(self.ui.preprocessPushButton, QtCore.SIGNAL("clicked()"), self.preprocessButtonPress)
    self.connect(self.ui.pushButtonBrowseAtomScore, QtCore.SIGNAL("clicked (bool)"), self.loadVolume)
    self.connect(self.ui.pushButton_2, QtCore.SIGNAL("clicked()"), self.generateAtomsButtonPress)
    self.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.pathwalkButtonPress)
    self.connect(self.ui.pushButton_5, QtCore.SIGNAL("clicked()"), self.deleteBonds)      
    self.connect(self.ui.pushButton_6, QtCore.SIGNAL("clicked()"), self.printDeleted)     

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

  def printDeleted(self):
       self.app.viewers['calpha'].main_chain.printDeletedBonds()
           
  def deleteBonds(self):
      self.app.viewers['calpha'].deleteSelectedBonds()
      self.app.viewers['calpha'].main_chain.unsetBonds()
      print self.app.viewers['calpha'].printDeletedBondAtoms()

  def createUi(self):
      self.ui = Ui_DialogPathwalker()
      self.ui.setupUi(self)      
      self.ui.preprocessPushButton.setEnabled(True)


  def loadVolume(self, temp):
    self.volumeName = str(QtGui.QFileDialog.getOpenFileName(self.app.viewers["volume"], self.app.viewers["volume"].tr("Open Volume"), "", self.app.viewers["volume"].tr(self.app.viewers["volume"].renderer.getSupportedLoadFileFormats())))
    self.ui.lineEditAtomScore.setText(self.volumeName)
    self.bringToFront()
 

  def preprocess(self):
      print str(self.ui.normalizeLine.text())
      subprocess.call(['python','EMAN2/bin/e2proc3d.py',self.volumeName, 'EMAN2/bin/map.mrc','--process',str(self.ui.normalizeLine.text()),'--process',str(self.ui.zeroThresholdLine.text())])
      self.app.viewers["volume"].loadDataFromFile("EMAN2/bin/map.mrc")



        