# Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
# Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
# Description:   A widget used ot perform binary skeletonization on a volume 


from PyQt4 import QtCore, QtGui
from ui_dialog_volume_binary_skeletonization import Ui_DialogVolumeBinarySkeletonization
from delayed_filter import DelayedFilter
from base_dialog_widget import BaseDialogWidget
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from seq_model.Chain import Chain

class VolumeBinarySkeletonizationForm(BaseDialogWidget):
    def __init__(self, main, volumeViewer, parent=None):
        BaseDialogWidget.__init__(self, main, 
                                  "&Binary Skeletonization", 
                                  "Apply binary skeletonization on the volume", 
                                  "perform_VolumeBinarySkeletonization", 
                                  "actions-volume-skeletonization-binary", 
                                  "actions-volume-skeletonization", 
                                  False, parent)
        self.app = main
        self.viewer = volumeViewer
        self.connect(self.viewer, QtCore.SIGNAL("modelLoaded()"), self.modelLoaded)
        self.connect(self.viewer, QtCore.SIGNAL("modelUnloaded()"), self.modelUnloaded)
        self.createUI()
        self.createActions()

    def createUI(self):
        self.ui = Ui_DialogVolumeBinarySkeletonization()
        self.ui.setupUi(self)
        
        self.connect(self.ui.horizontalSliderIsoLevel,QtCore.SIGNAL("valueChanged(int)"),self.isoValueChanged)
        self.connect(self.ui.comboBoxMethod, QtCore.SIGNAL("currentIndexChanged (int)"), self.methodChanged)                                            
        self.methodChanged(0)
        
    def createActions(self):               
        self.displayAct.setEnabled(False)  
        self.connect(self.displayAct, QtCore.SIGNAL("triggered()"), self.loadSelf)  
        
    def modelLoaded(self):
        maxDensity = self.viewer.renderer.getMaxDensity()
        minDensity = self.viewer.renderer.getMinDensity()
        self.ui.horizontalSliderIsoLevel.setMinimum(int(minDensity*100))
        self.ui.horizontalSliderIsoLevel.setMaximum(int(maxDensity*100))
        defaultDensity = (int(minDensity*100) + int(maxDensity*100.0)) / 2
        self.ui.horizontalSliderIsoLevel.setValue(defaultDensity)
        self.displayAct.setEnabled(True)  
    
    def modelUnloaded(self):
        self.displayAct.setEnabled(False)  
        self.close()

    def isoValueChanged(self, newLevel):
        self.ui.labelIsoValueDisplay.setNum(newLevel/100.0)

    def getDensityThreshold(self):
        return self.ui.horizontalSliderIsoLevel.value()/100.0
    
    def getSkeletonizationMethod(self):
        return self.ui.comboBoxMethod.currentIndex()
    
    def getMinCurveLength(self):
        return self.ui.spinBoxMinCurve.value()

    def getMinSurfaceSize(self):
        return self.ui.spinBoxMinSurface.value()


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
          
    def accept(self):
        if(self.viewer.loaded):
            self.setCursor(QtCore.Qt.BusyCursor)
            method = self.getSkeletonizationMethod()
            if(method == 0):
                skeleton = self.viewer.renderer.performBinarySkeletonizationJu2007(self.getDensityThreshold(), self.getMinCurveLength(), self.getMinSurfaceSize())
                #self.viewer.renderer.loadExtremal(skeleton);
                #self.app.viewers["volume"].loadExtremal(skeleton);
                self.app.viewers["skeleton"].loadVolume(skeleton)
                #self.generateAtoms("extremal.pdb")
                #bonds1 = skeleton.getExtremalBonds1()
                #bonds2 = skeleton.getExtremalBonds2()
                #print len(bonds1)
                #for i in range(len(bonds1)):
                    #print str(bonds1[i]) + " " + str(bonds2[i])
                #    self.app.viewers['calpha'].renderer.drawAddedBond(bonds1[i], bonds2[i])
                #self.app.viewers['calpha'].emitModelChanged()

            self.setCursor(QtCore.Qt.ArrowCursor)
            self.close() 
        else:
            QtGui.QMessageBox.critical(None, "Source volume unloaded", "A volume must be loaded to perform skeletonization", QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
        BaseDialogWidget.accept(self)

    def loadSelf(self):
        self.ui.horizontalSliderIsoLevel.setValue(int(self.viewer.getIsoValue()*100))
        
    def getCitationHtml(self, title, author, journal):
        return "<b>" + title + "</b><br>" + author + "<br><i>" + journal + "</i>"
                  
    def methodChanged(self, id):
        if(id == 0):
            self.ui.textCitation.setHtml(self.getCitationHtml("Computing a family of skeletons of volumetric models for shape description", "Tao Ju, Matthew Baker and Wah Chiu", "Computer Aided Design 39(5) 2007"))
        else :
            self.ui.textCitation.setHtml("")
