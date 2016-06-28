# Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
# Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
# Description:   A widget used ot perform binary skeletonization on a volume 


from PyQt4 import QtCore, QtGui
from extremalCurve import Ui_Extremal_Curve
from delayed_filter import DelayedFilter
from base_dialog_widget import BaseDialogWidget
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from seq_model.Chain import Chain

class ExtremeCurveForm(BaseDialogWidget):
    def __init__(self, main, volumeViewer, parent=None):
        BaseDialogWidget.__init__(self, main, 
                                  "&Extremal Curve", 
                                  "Apply extremal curve skeletonization to density map", 
                                  "perform_VolumeBinarySkeletonization", 
                                  "actions-volume-skeletonization-extremal", 
                                  "actions-volume-skeletonization", 
                                  False, parent)
        self.app = main
        self.viewer = volumeViewer
        #self.connect(self.viewer, QtCore.SIGNAL("modelLoaded()"), self.modelLoaded)
        #self.connect(self.viewer, QtCore.SIGNAL("modelUnloaded()"), self.modelUnloaded)
        self.createUI()
        self.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.findExtremalCurve)

        #self.createActions()

    def createUI(self):
      self.ui = Ui_Extremal_Curve()
      self.ui.setupUi(self)

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

    def findExtremalCurve(self):
      skeleton = self.viewer.renderer.performExtremalCurve2016();
      self.app.viewers["skeleton"].loadVolume(skeleton)
      self.generateAtoms("extremal.pdb")
      bonds1 = skeleton.getExtremalBonds1()
      bonds2 = skeleton.getExtremalBonds2()
      #print len(bonds1)
      for i in range(len(bonds1)):
                    #print str(bonds1[i]) + " " + str(bonds2[i])
        self.app.viewers['calpha'].renderer.drawAddedBond(bonds1[i], bonds2[i])
      self.app.viewers['calpha'].emitModelChanged()
