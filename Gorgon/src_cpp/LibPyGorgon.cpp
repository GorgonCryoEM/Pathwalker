// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi Abeysinghe (sasakthi@gmail.com)
// Description:   Boosts the interfaces required for the python GORGON interface.




#ifndef GORGON_LIB_PY_GORGON
#define GORGON_LIB_PY_GORGON

#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1

#ifndef GL_GLEXT_PROTOTYPES
	#ifdef __linux
		#define GL_GLEXT_PROTOTYPES
	#endif
#endif

#include <Gorgon/SSECorrespondenceEngine.h>
#include <Gorgon/VolumeRenderer.h>
#include <Gorgon/MeshRenderer.h>
#include <Gorgon/SSERenderer.h>
#include <Gorgon/Renderer.h>
#include <Gorgon/InteractiveSkeletonEngine.h>
#include <Gorgon/InteractiveLoopBuilderEngine.h>
#include <Gorgon/FlexibleFittingEngine.h>
#include <Gorgon/CAlphaRenderer.h>
#include <MathTools/Vector3D.h>
#include <MathTools/LinearSolver.h>
#include <GraphMatch/PDBAtom.h>
#include <GraphMatch/LinkedNode.h>
#include <GraphMatch/PDBBond.h>
#include <GraphMatch/SEQReader.h>
#include <ProteinMorph/SSEHunter.h>
#include <ProteinMorph/SSECorrespondenceNode.h>
#include <boost/python.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/tuple.hpp>
#include <boost/tuple/tuple.hpp>


using namespace std;
using namespace boost::python;
using namespace wustl_mm::Visualization;
using namespace wustl_mm::GraphMatch;
using namespace wustl_mm::SkeletonMaker;
using namespace wustl_mm::Protein_Morph;
using namespace boost::tuples;
using wustl_mm::MathTools::LinearSolver;



// ********************** From EMAN2 typeconverter.h ************************************
#include <boost/python/to_python_converter.hpp>
namespace python = boost::python;



template <class T>
struct vector_from_python {
	vector_from_python() {
		python::converter::registry::push_back(&convertible, &construct, python::type_id<std::vector<T> >());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!(PyList_Check(obj_ptr) || PyTuple_Check(obj_ptr) || PyIter_Check(obj_ptr)  || PyRange_Check(obj_ptr))) {
			return 0;
		}
		return obj_ptr;
	}


	static void construct(PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data) {
		void* storage = ((python::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes;
		new (storage) std::vector<T>();
		data->convertible = storage;
		std::vector<T>& result = *((std::vector<T>*) storage);
		python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
		while(1) {
			python::handle<> py_elem_hdl(python::allow_null(PyIter_Next(obj_iter.get())));
			if (PyErr_Occurred()) {
				python::throw_error_already_set();
			}
			if (!py_elem_hdl.get()) {
				break;
			}

			python::object py_elem_obj(py_elem_hdl);
			python::extract<T> elem_proxy(py_elem_obj);
			result.push_back(elem_proxy());
		}
	}
};

template<class T, class T2>
struct tuple3_from_python {
	tuple3_from_python() {
		python::converter::registry::push_back(&convertible, &construct,
											   python::type_id<T>());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!(PyList_Check(obj_ptr) || PyTuple_Check(obj_ptr) || PyIter_Check(obj_ptr)  || PyRange_Check(obj_ptr))) {
			return 0;
		}
		return obj_ptr;
	}


	static void construct(PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data)	{
		void* storage = ((python::converter::rvalue_from_python_storage<T>*)data)->storage.bytes;
		new (storage) T();
		data->convertible = storage;
		T& result = *((T*) storage);
		python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
		int i = 0;

		while(1) {
			python::handle<> py_elem_hdl(python::allow_null(PyIter_Next(obj_iter.get())));
			if (PyErr_Occurred()) {
				python::throw_error_already_set();
			}

			if (!py_elem_hdl.get()) {
				break;
			}

			python::object py_elem_obj(py_elem_hdl);
			python::extract<T2> elem_proxy(py_elem_obj);
			result[i] = elem_proxy();
			i++;
		}
	}
};

template <class T>
struct tuple3_to_python : python::to_python_converter<T, tuple3_to_python<T> > {
	static PyObject* convert(T const& p) {
		python::tuple result = python::make_tuple(p[0], p[1], p[2]);
		return python::incref(python::tuple(result).ptr());
	}
};

// **************************************************************************************

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
BOOST_PYTHON_MODULE(libpyGORGON)
{
	
	vector_from_python<float>();
	class_<std::vector<float> >("std::vector<float>")
	        .def(vector_indexing_suite<std::vector<float> >() );
	class_<std::vector<std::vector<float> > >("std::vector<std::vector<float> >")
	        .def(vector_indexing_suite<std::vector<std::vector<float> > >() );

	class_<std::vector<int> >("std::vector<int>")
	        .def(vector_indexing_suite<std::vector<int> >() );
	class_<std::vector<std::vector<int> > >("std::vector<std::vector<int> >")
	        .def(vector_indexing_suite<std::vector<std::vector<int> > >() );

	vector_from_python< std::vector<float> >();

	vector_from_python<int>();
	vector_from_python< std::vector<int> >();

	vector_from_python<bool>();
	vector_from_python< std::vector<bool> >();

    vector_from_python<unsigned long long>();

	vector_from_python<Vector3DFloat>();
	vector_from_python< std::vector<Vector3DFloat> >();


//	tuple3_to_python<Vector3DInt>();
//	tuple3_to_python<Vector3DFloat>();
	tuple3_from_python<Vector3DInt, int>();
	tuple3_from_python<Vector3DFloat, float>();

	class_<Vector3DFloat>("Vector3DFloat", init<float, float, float>())
		.def("x", &Vector3DFloat::X)
		.def("y", &Vector3DFloat::Y)
		.def("z", &Vector3DFloat::Z)
		.def("length", &Vector3DFloat::Length)
		.def("Transform", &Vector3DFloat::Transform)
		.def("rotate", &Vector3DFloat::Rotate)
		.def("Print", &Vector3DFloat::Print)
		.def(self * double())
		.def(self + self)
		.def(self - self)
		.def(self ^ self)
		.def(self += self)
		.def(self -= self)
		//Vector3DTemplate<T> operator-();
		//Vector3DTemplate<T>& operator=(const Vector3DTemplate<T>& d);

	;

	class_<Vector3DInt>("Vector3DInt", init<int, int, int>())
		.def("x", &Vector3DInt::X)
		.def("y", &Vector3DInt::Y)
		.def("z", &Vector3DInt::Z)
		.def("length", &Vector3DInt::Length)
		.def(self * double())
		.def(self + self)
		.def(self - self)
		.def(self ^ self)
		.def(self += self)
		.def(self -= self)
		
	;

	class_<MatrixFloat>("MatrixFloat", init<int, int>())
		.def("getValue", &MatrixFloat::GetValue)
		.def("setValue", &MatrixFloat::SetValue)
		.def(self * float())
		.def(self * self)
		.def(self + self)
	;
	


	class_<PDBAtom>("PDBAtom", init<string, char, unsigned int, string>())
		.def("getPDBId", &PDBAtom::GetPDBId)
		.def("getSerial", &PDBAtom::GetSerial)
		.def("getName", &PDBAtom::GetName)
		.def("getAltLoc", &PDBAtom::GetAltLoc)
		.def("getResName", &PDBAtom::GetResName)
		.def("getChainId", &PDBAtom::GetChainId)
		.def("getResSeq", &PDBAtom::GetResSeq)
		.def("getICode", &PDBAtom::GetICode)
		.def("getPosition", &PDBAtom::GetPosition)
		.def("getOccupancy", &PDBAtom::GetOccupancy)
		.def("getTempFactor", &PDBAtom::GetTempFactor)
		.def("getElement", &PDBAtom::GetElement)
		.def("getCharge", &PDBAtom::GetCharge)
		.def("getAtomRadius", &PDBAtom::GetAtomRadius)
		.def("getColorA", &PDBAtom::GetColorA)
		.def("getColorR", &PDBAtom::GetColorR)
		.def("getColorG", &PDBAtom::GetColorG)
		.def("getColorB", &PDBAtom::GetColorB)
		.def("getSelected", &PDBAtom::GetSelected)
		.def("getCorrelationScore", &PDBAtom::GetCorrelationScore)
		.def("getSkeletonScore", &PDBAtom::GetSkeletonScore)
		.def("getGeometryScore", &PDBAtom::GetGeometryScore)
		.def("getHashKey", &PDBAtom::GetHashKey)
		.def("getVisible", &PDBAtom::GetVisible)
		.def("constructHashKey", &PDBAtom::ConstructHashKey)
		.def("setSerial", &PDBAtom::SetSerial)
		.def("setName", &PDBAtom::SetName)
		.def("setAltLoc", &PDBAtom::SetAltLoc)
		.def("setResName", &PDBAtom::SetResName)
		.def("setChainId", &PDBAtom::SetChainId)
		.def("setResSeq", &PDBAtom::SetResSeq)
		.def("setICode", &PDBAtom::SetICode)
		.def("setPosition", &PDBAtom::SetPosition)
		.def("setOccupancy", &PDBAtom::SetOccupancy)
		.def("setTempFactor", &PDBAtom::SetTempFactor)
		.def("setElement", &PDBAtom::SetElement)
		.def("setCharge", &PDBAtom::SetCharge)
		.def("setAtomRadius", &PDBAtom::SetAtomRadius)
		.def("setColor", &PDBAtom::SetColor)
		.def("setSelected", &PDBAtom::SetSelected)
		.def("setVisible", &PDBAtom::SetVisible)
		.def("setCorrelationScore", &PDBAtom::SetCorrelationScore)
		.def("setSkeletonScore", &PDBAtom::SetSkeletonScore)
		.def("setGeometryScore", &PDBAtom::SetGeometryScore)
		.def("transform", &PDBAtom::Transform)
		.def("interpolateTransform", &PDBAtom::InterpolateTransform)
		.def("getInterpolateTransformLocation", &PDBAtom::GetInterpolateTransformLocation)
		.def("setPrevCAHash", &PDBAtom::SetPrevCAHash)
		.def("setNextCAHash", &PDBAtom::SetNextCAHash)
		.def("getDeletedBondAtom", &PDBAtom::GetDeletedBondAtom)
	;
	

	class_<PDBBond>("PDBBond", init<>())
		.def("getAtom0Ix", &PDBBond::GetAtom0Ix)
		.def("getAtom1Ix", &PDBBond::GetAtom1Ix)
		.def("getSelected", &PDBBond::GetSelected)
		.def("setAtom0Ix", &PDBBond::SetAtom0Ix)
		.def("setAtom1Ix", &PDBBond::SetAtom1Ix)
		.def("setSelected", &PDBBond::SetSelected)
	;

	double (Volume::*getDataAt1)(int, int, int)	= &Volume::getDataAt;
	
	class_<Volume>("Volume", init<int,int,int>())
		.def("getOriginX", &Volume::getOriginX)
		.def("getOriginY", &Volume::getOriginY)
		.def("getOriginZ", &Volume::getOriginZ)
		.def("getSizeX", &Volume::getSizeX)
		.def("getSizeY", &Volume::getSizeY)
		.def("getSizeZ", &Volume::getSizeZ)
		.def("getSpacingX", &Volume::getSpacingX)
		.def("getSpacingY", &Volume::getSpacingY)
		.def("getSpacingZ", &Volume::getSpacingZ)
		.def("getMin", &Volume::getMin)
		.def("getMax", &Volume::getMax)
		.def("getDataAt", getDataAt1)
		.def("buildHistogram", &Volume::buildHistogram)
		.def("getHistogramBinValue", &Volume::getHistogramBinValue)
		.def("getExtremalBonds1", &Volume::getExtremalBonds1)
		.def("getExtremalBonds2", &Volume::getExtremalBonds2)
		.def("getMinCurveBonds1", &Volume::getMinCurveBonds1)
		.def("getMinCurveBonds2", &Volume::getMinCurveBonds2)
		.def("getSaddleCurveBonds1", &Volume::getSaddleCurveBonds1)
		.def("getSaddleCurveBonds2", &Volume::getSaddleCurveBonds2)
		.def("setExtremalParams", &Volume::setExtremalParams)
		.def("hideSegments", &Volume::hideSegments)
		.def("setExtremalWidth", &Volume::setExtremalWidth)
		.def("setExtremalHeight", &Volume::setExtremalHeight)
		.def("setExtremalSlices", &Volume::setExtremalSlices)
		.def("setExtremalGridResolution", &Volume::setExtremalGridResolution)
		.def("setExtremalSlicesThickness", &Volume::setExtremalSlicesThickness)
		.def("setExtremalIsovalue", &Volume::setExtremalIsovalue)
		.def("setExtremalChecks", &Volume::setExtremalChecks)
		.def("getDisplayQuads", &Volume::getDisplayQuads)
		.def("getExtremalNormals", &Volume::getExtremalNormals)
		.def("getQuadTypes", &Volume::getQuadTypes)
		.def("hideSurfaces", &Volume::hideSurfaces)
		.def("getQuadSaliencies", &Volume::getQuadSaliencies)
		.def("getQuadIntensities", &Volume::getQuadIntensities)
		.def("getQuadEigenvalues", &Volume::getQuadEigenvalues)
		.def("setMaxPointOn", &Volume::setMaxPointOn)
		.def("setMinPointOn", &Volume::setMinPointOn)
		.def("setSaddlePointOn", &Volume::setSaddlePointOn)
		.def("getExtremalMaxPoints", &Volume::getExtremalMaxPoints)
		.def("getExtremalMinPoints", &Volume::getExtremalMinPoints)
		.def("getExtremalSaddlePoints", &Volume::getExtremalSaddlePoints)
		.def("setXMin", &Volume::setXMin)
		.def("setYMin", &Volume::setYMin)
		.def("setZMin", &Volume::setZMin)
		.def("setXMax", &Volume::setXMax)
		.def("setYMax", &Volume::setYMax)
		.def("setZMax", &Volume::setZMax)
		.def("getMaxSaliencies", &Volume::getMaxSaliencies)
		.def("getMaxIntensities", &Volume::getMaxIntensities)
		.def("getMaxEigenvalues", &Volume::getMaxEigenvalues)
		.def("getMinSaliencies", &Volume::getMinSaliencies)
		.def("getMinIntensities", &Volume::getMinIntensities)
		.def("getMinEigenvalues", &Volume::getMinEigenvalues)
		.def("getSaddleSaliencies", &Volume::getSaddleSaliencies)
		.def("getSaddleIntensities", &Volume::getSaddleIntensities)
		.def("getSaddleEigenvalues", &Volume::getSaddleEigenvalues)
		.def("getMinI", &Volume::getMinI)
		.def("getMaxI", &Volume::getMaxI)
		.def("getMinE", &Volume::getMinE)
		.def("getMaxE", &Volume::getMaxE)
		.def("getMaxHashes", &Volume::getMaxHashes)
		.def("getMinHashes", &Volume::getMinHashes)
		.def("getSaddleHashes", &Volume::getSaddleHashes)
		.def("getVertexPos", &Volume::getVertexPos)
		.def("getQuadIndices", &Volume::getQuadIndices)
		.def("getMaxPointSaliencies", &Volume::getMaxPointSaliencies)
		.def("getMaxPointIntensities", &Volume::getMaxPointIntensities)
		.def("getMaxPointEigenvalues", &Volume::getMaxPointEigenvalues)
		.def("getMinPointSaliencies", &Volume::getMinPointSaliencies)
		.def("getMinPointIntensities", &Volume::getMinPointIntensities)
		.def("getMinPointEigenvalues", &Volume::getMinPointEigenvalues)
		.def("getSaddlePointSaliencies", &Volume::getSaddlePointSaliencies)
		.def("getSaddlePointIntensities", &Volume::getSaddlePointIntensities)
		.def("getSaddlePointEigenvalues", &Volume::getSaddlePointEigenvalues)
		.def("getMinPointEigenvalues0", &Volume::getMinPointEigenvalues0)
		.def("getMinPointEigenvalues1", &Volume::getMinPointEigenvalues1)
		.def("getMaxPointEigenvalues0", &Volume::getMaxPointEigenvalues0)
		.def("getMaxPointEigenvalues1", &Volume::getMaxPointEigenvalues1)
		.def("getSaddlePointEigenvalues0", &Volume::getSaddlePointEigenvalues0)
		.def("getSaddlePointEigenvalues1", &Volume::getSaddlePointEigenvalues1)
		.def("getMaxCurveEigenvalues0", &Volume::getMaxCurveEigenvalues0)
		.def("getMaxCurveEigenvalues1", &Volume::getMaxCurveEigenvalues1)
		.def("getMinCurveEigenvalues0", &Volume::getMinCurveEigenvalues0)
		.def("getMinCurveEigenvalues1", &Volume::getMinCurveEigenvalues1)
		.def("getSaddleCurveEigenvalues0", &Volume::getSaddleCurveEigenvalues0)
		.def("getSaddleCurveEigenvalues1", &Volume::getSaddleCurveEigenvalues1)
		.def("getQuadEigenvalues0", &Volume::getQuadEigenvalues0)
		.def("getQuadEigenvalues1", &Volume::getQuadEigenvalues1)
		.def("getMaxCurveEigenvectors", &Volume::getMaxCurveEigenvectors)
		.def("getMinCurveEigenvectors", &Volume::getMinCurveEigenvectors)
		.def("getSaddleCurveEigenvectors", &Volume::getSaddleCurveEigenvectors)
		.def("getQuadEigenvectors", &Volume::getQuadEigenvectors)
		.def("getMaxPointEigenvectors", &Volume::getMaxPointEigenvectors)
		.def("getMinPointEigenvectors", &Volume::getMinPointEigenvectors)
		.def("getSaddlePointEigenvectors", &Volume::getSaddlePointEigenvectors)
	;


	class_<NonManifoldMesh_Annotated>("NonManifoldMesh_Annotated", init<>())
	;
	
	class_<Renderer>("Renderer", init<>())
		.def("draw", &Renderer::Draw)
		.def("drawBoundingBox", &Renderer::DrawBoundingBox)
		.def("loadFile", &Renderer::LoadFile)
		.def("saveFile", &Renderer::SaveFile)
		.def("unload", &Renderer::Unload)
		.def("selectionRotate", &Renderer::SelectionRotate)
		.def("selectionObjectCount", &Renderer::SelectionObjectCount)
		.def("selectionCenterOfMass", &Renderer::SelectionCenterOfMass)
		.def("selectionMove", &Renderer::SelectionMove)
		.def("selectionClear", &Renderer::SelectionClear)		
		.def("selectionToggle", &Renderer::SelectionToggle)	
		.def("getSupportedLoadFileFormats", &Renderer::GetSupportedLoadFileFormats)
		.def("getSupportedSaveFileFormats", &Renderer::GetSupportedSaveFileFormats)
		.def("getMin", &Renderer::GetMin)
		.def("getMax", &Renderer::GetMax)
		.def("get3DCoordinates", &Renderer::Get3DCoordinates)
		.def("setCuttingPlane", &Renderer::SetCuttingPlane)	
		.def("setSpacing", &Renderer::SetSpacing)
		.def("getSpacingX", &Renderer::GetSpacingX)
		.def("getSpacingY", &Renderer::GetSpacingY)
		.def("getSpacingZ", &Renderer::GetSpacingZ)
		.def("setOrigin", &Renderer::SetOrigin)
		.def("getOriginX", &Renderer::GetOriginX)
		.def("getOriginY", &Renderer::GetOriginY)
		.def("getOriginZ", &Renderer::GetOriginZ)
		.def("setDisplayStyle", &Renderer::SetDisplayStyle)
		.def("setObjectSpecificColoring", &Renderer::SetObjectSpecificColoring)
		.def("updateBoundingBox", &Renderer::UpdateBoundingBox)
		.def("clearOtherHighlights", &Renderer::ClearOtherHighlights)
	;
	
	
	class_< VolumeRenderer, bases<Renderer> >("VolumeRenderer", init<>())		
		.def("draw", &VolumeRenderer::Draw)
		.def("drawBoundingBox", &VolumeRenderer::DrawBoundingBox)
		.def("enableDraw", &VolumeRenderer::EnableDraw)
		.def("loadFile", &VolumeRenderer::LoadFile)
		.def("loadFileRAW", &VolumeRenderer::LoadFileRAW)
		.def("saveFile", &VolumeRenderer::SaveFile)
		.def("unload", &VolumeRenderer::Unload)
		.def("selectionRotate", &VolumeRenderer::SelectionRotate)
		.def("selectionObjectCount", &VolumeRenderer::SelectionObjectCount)
		.def("selectionCenterOfMass", &VolumeRenderer::SelectionCenterOfMass)
		.def("selectionMove", &VolumeRenderer::SelectionMove)
		.def("selectionClear", &VolumeRenderer::SelectionClear)		
		.def("selectionToggle", &VolumeRenderer::SelectionToggle)	
		.def("getSupportedLoadFileFormats", &VolumeRenderer::GetSupportedLoadFileFormats)
		.def("getSupportedSaveFileFormats", &VolumeRenderer::GetSupportedSaveFileFormats)
		.def("getMin", &VolumeRenderer::GetMin)
		.def("getMax", &VolumeRenderer::GetMax)
		.def("getMinDensity", &VolumeRenderer::GetMinDensity)
		.def("getMaxDensity", &VolumeRenderer::GetMaxDensity)
		.def("getVolume", &VolumeRenderer::GetVolume, return_value_policy<reference_existing_object>())		
		.def("getSurfaceValue", &VolumeRenderer::GetSurfaceValue)
		.def("get3DCoordinates", &VolumeRenderer::Get3DCoordinates)
		.def("setSurfaceValue", &VolumeRenderer::SetSurfaceValue)		
		.def("setMaxSurfaceValue", &VolumeRenderer::SetMaxSurfaceValue)		
		.def("setSampleInterval", &VolumeRenderer::SetSampleInterval)
		.def("setDisplayRadius", &VolumeRenderer::SetDisplayRadius)
		.def("setDisplayRadiusOrigin", &VolumeRenderer::SetDisplayRadiusOrigin)		
		.def("setCuttingPlane", &VolumeRenderer::SetCuttingPlane)
		.def("setViewingType", &VolumeRenderer::SetViewingType)
		.def("normalizeVolume", &VolumeRenderer::NormalizeVolume)
		.def("downsampleVolume", &VolumeRenderer::DownsampleVolume)
		.def("cropVolume", &VolumeRenderer::CropVolume)
		.def("performBinarySkeletonizationJu2007", &VolumeRenderer::PerformBinarySkeletonizationJu2007, return_value_policy<reference_existing_object>())		
		.def("performGrayscaleSkeletonizationAbeysinghe2008", &VolumeRenderer::PerformGrayscaleSkeletonizationAbeysinghe2008, return_value_policy<reference_existing_object>())		
		.def("performPreservingGrayscaleSkeletonizationAbeysinghe2008", &VolumeRenderer::PerformPreservingGrayscaleSkeletonizationAbeysinghe2008, return_value_policy<reference_existing_object>())		
		.def("setSpacing", &VolumeRenderer::SetSpacing)
		.def("getSpacingX", &VolumeRenderer::GetSpacingX)
		.def("getSpacingY", &VolumeRenderer::GetSpacingY)
		.def("getSpacingZ", &VolumeRenderer::GetSpacingZ)
		.def("setOrigin", &VolumeRenderer::SetOrigin)
		.def("getOriginX", &VolumeRenderer::GetOriginX)
		.def("getOriginY", &VolumeRenderer::GetOriginY)
		.def("getOriginZ", &VolumeRenderer::GetOriginZ)
		.def("useDisplayRadius", &VolumeRenderer::UseDisplayRadius)
		.def("performSmoothLaplacian", &VolumeRenderer::PerformSmoothLaplacian)		
		.def("setDisplayStyle", &VolumeRenderer::SetDisplayStyle)
		.def("setObjectSpecificColoring", &VolumeRenderer::SetObjectSpecificColoring)
		.def("updateBoundingBox", &VolumeRenderer::UpdateBoundingBox)
		.def("loadExtremal", &VolumeRenderer::loadExtremal)
		.def("performExtremalCurve2016", &VolumeRenderer::PerformExtremalCurve2016, return_value_policy<reference_existing_object>())
		
	;

	class_< MeshRenderer, bases<Renderer> >("MeshRenderer", init<>())
		.def("draw", &MeshRenderer::Draw)
		.def("drawBoundingBox", &MeshRenderer::DrawBoundingBox)
		.def("loadFile", &MeshRenderer::LoadFile)
		.def("loadVolume", &MeshRenderer::LoadVolume)
		.def("saveFile", &MeshRenderer::SaveFile)
		.def("unload", &MeshRenderer::Unload)
		.def("selectionRotate", &MeshRenderer::SelectionRotate)
		.def("selectionObjectCount", &MeshRenderer::SelectionObjectCount)
		.def("selectionCenterOfMass", &MeshRenderer::SelectionCenterOfMass)
		.def("selectionMove", &MeshRenderer::SelectionMove)
		.def("selectionClear", &MeshRenderer::SelectionClear)		
		.def("selectionToggle", &MeshRenderer::SelectionToggle)	
		.def("performSmoothLaplacian", &MeshRenderer::PerformSmoothLaplacian)
		.def("getSupportedLoadFileFormats", &MeshRenderer::GetSupportedLoadFileFormats)
		.def("getSupportedSaveFileFormats", &MeshRenderer::GetSupportedSaveFileFormats)
		.def("getMin", &MeshRenderer::GetMin)
		.def("getMax", &MeshRenderer::GetMax)
		.def("getMesh", &MeshRenderer::GetMesh, return_value_policy<reference_existing_object>())
		.def("get3DCoordinates", &MeshRenderer::Get3DCoordinates)
		.def("intersectMeshAndSphere", &MeshRenderer::IntersectMeshAndSphere)
		.def("getIntersectionPoint", &MeshRenderer::GetIntersectionPoint)
		.def("setLineThickness", &MeshRenderer::SetLineThickness)
		.def("setSpacing", &MeshRenderer::SetSpacing)
		.def("getSpacingX", &MeshRenderer::GetSpacingX)
		.def("getSpacingY", &MeshRenderer::GetSpacingY)
		.def("getSpacingZ", &MeshRenderer::GetSpacingZ)
		.def("setOrigin", &MeshRenderer::SetOrigin)
		.def("getOriginX", &MeshRenderer::GetOriginX)
		.def("getOriginY", &MeshRenderer::GetOriginY)
		.def("getOriginZ", &MeshRenderer::GetOriginZ)
		.def("setDisplayStyle", &MeshRenderer::SetDisplayStyle)
		.def("setObjectSpecificColoring", &MeshRenderer::SetObjectSpecificColoring)
		.def("updateBoundingBox", &MeshRenderer::UpdateBoundingBox)
	;

	class_< SSERenderer, bases<Renderer> >("SSERenderer", init<>())
		.def("draw", &SSERenderer::Draw)
		.def("drawBoundingBox", &SSERenderer::DrawBoundingBox)
		.def("loadHelixFile", &SSERenderer::LoadHelixFile)
		.def("loadSheetFile", &SSERenderer::LoadSheetFile)
		.def("saveHelixFile", &SSERenderer::SaveHelixFile)
		.def("saveSheetFile", &SSERenderer::SaveSheetFile)
		.def("loadGraphSSE", &SSERenderer::LoadGraphSSE)
		.def("unloadGraphSSEs", &SSERenderer::UnloadGraphSSEs)
		.def("unload", &SSERenderer::Unload)
		.def("selectionRotate", &SSERenderer::SelectionRotate)
		.def("selectionObjectCount", &SSERenderer::SelectionObjectCount)
		.def("selectionCenterOfMass", &SSERenderer::SelectionCenterOfMass)
		.def("selectionMove", &SSERenderer::SelectionMove)
		.def("selectionClear", &SSERenderer::SelectionClear)		
		.def("selectionToggle", &SSERenderer::SelectionToggle)	
		.def("getSupportedHelixLoadFileFormats", &SSERenderer::GetSupportedHelixLoadFileFormats)
		.def("getSupportedHelixSaveFileFormats", &SSERenderer::GetSupportedHelixSaveFileFormats)
		.def("getSupportedSheetLoadFileFormats", &SSERenderer::GetSupportedSheetLoadFileFormats)
		.def("getSupportedSheetSaveFileFormats", &SSERenderer::GetSupportedSheetSaveFileFormats)
		.def("getMin", &SSERenderer::GetMin)
		.def("getMax", &SSERenderer::GetMax)
		.def("get3DCoordinates", &SSERenderer::Get3DCoordinates)
		.def("setHelixColor", &SSERenderer::SetHelixColor)		
		.def("setSheetColor", &SSERenderer::SetSheetColor)		
		.def("setSSEColor", &SSERenderer::SetSSEColor)		
		.def("setSpacing", &SSERenderer::SetSpacing)
		.def("getSpacingX", &SSERenderer::GetSpacingX)
		.def("getSpacingY", &SSERenderer::GetSpacingY)
		.def("getSpacingZ", &SSERenderer::GetSpacingZ)
		.def("setOrigin", &SSERenderer::SetOrigin)
		.def("getOriginX", &SSERenderer::GetOriginX)
		.def("getOriginY", &SSERenderer::GetOriginY)
		.def("getOriginZ", &SSERenderer::GetOriginZ)
		.def("fitSelectedSSEs", &SSERenderer::FitSelectedSSEs)		
		.def("startNewSSE", &SSERenderer::StartNewSSE)
		.def("addSSEPoint", &SSERenderer::AddSSEPoint)
		.def("finalizeSheet", &SSERenderer::FinalizeSheet)
		.def("finalizeHelix", &SSERenderer::FinalizeHelix)
		.def("addHelix", &SSERenderer::AddHelix)
		.def("setObjectSpecificColoring", &SSERenderer::SetObjectSpecificColoring)
		.def("removeSelectedSSEs", &SSERenderer::RemoveSelectedSSEs)
		.def("removeHelices", &SSERenderer::RemoveHelices)
		.def("getHelixCount", &SSERenderer::GetHelixCount)
		.def("getHelixCorner", &SSERenderer::GetHelixCorner)
		.def("setDisplayStyle", &SSERenderer::SetDisplayStyle)
		.def("updateBoundingBox", &SSERenderer::UpdateBoundingBox)
		.def("getSelectedHelixIndices", &SSERenderer::GetSelectedHelixIndices)
		.def("setSSEOrientationFlips", &SSERenderer::SetSSEOrientationFlips)
		.def("setHelixCorrs", &SSERenderer::SetHelixCorrs)
		.def("setSelectedPDBHelices", &SSERenderer::SetSelectedPDBHelices)
		.def("clearOtherHighlights", &SSERenderer::ClearOtherHighlights)
	;


	class_< CAlphaRenderer, bases<Renderer> >("CAlphaRenderer", init<>())
		.def("draw", &CAlphaRenderer::Draw)
		.def("drawBoundingBox", &CAlphaRenderer::DrawBoundingBox)
		.def("loadFile", &CAlphaRenderer::LoadFile)
		.def("loadSSEHunterFile", &CAlphaRenderer::LoadSSEHunterFile)
		.def("saveSSEHunterFile", &CAlphaRenderer::SaveSSEHunterFile)
//		.def("getSSEHunterAtoms", &CAlphaRenderer::GetSSEHunterAtoms)
		.def("updateTotalScoreSSEHunterAtoms", &CAlphaRenderer::UpdateTotalScoreSSEHunterAtoms)
		.def("colorSSEHunterAtoms", &CAlphaRenderer::ColorSSEHunterAtoms)
		.def("saveFile", &CAlphaRenderer::SaveFile)
		.def("unload", &CAlphaRenderer::Unload)
		.def("selectionRotate", &CAlphaRenderer::SelectionRotate)
		.def("selectionObjectCount", &CAlphaRenderer::SelectionObjectCount)
		.def("selectionAtomCount", &CAlphaRenderer::SelectionAtomCount)
		.def("selectionCenterOfMass", &CAlphaRenderer::SelectionCenterOfMass)
		.def("selectionMove", &CAlphaRenderer::SelectionMove)
		.def("selectionClear", &CAlphaRenderer::SelectionClear)		
		.def("selectionToggle", &CAlphaRenderer::SelectionToggle)	
		.def("getAtomFromHitStack", &CAlphaRenderer::GetAtomFromHitStack, return_value_policy<reference_existing_object>())
		.def("getSupportedLoadFileFormats", &CAlphaRenderer::GetSupportedLoadFileFormats)
		.def("getSupportedSaveFileFormats", &CAlphaRenderer::GetSupportedSaveFileFormats)
		.def("getMin", &CAlphaRenderer::GetMin)
		.def("getMax", &CAlphaRenderer::GetMax)
		.def("get3DCoordinates", &CAlphaRenderer::Get3DCoordinates)
		.def("addAtom", &CAlphaRenderer::AddAtom, return_value_policy<reference_existing_object>())
		.def("getAtom", &CAlphaRenderer::GetAtom, return_value_policy<reference_existing_object>())
		.def("getSelectedAtom", &CAlphaRenderer::GetSelectedAtom, return_value_policy<reference_existing_object>())
		.def("getAtomCount", &CAlphaRenderer::GetAtomCount)
    .def("getAtomHashes", &CAlphaRenderer::GetAtomHashes)
		.def("deleteAtom", &CAlphaRenderer::DeleteAtom)	
		.def("getDeletedBondAtoms", &CAlphaRenderer::getDeletedBondAtoms)
		.def("getDeletedBonds1Ix", &CAlphaRenderer::getDeletedBonds1Ix)		
		.def("addBond", &CAlphaRenderer::AddBond)
		.def("getBond", &CAlphaRenderer::GetBond, return_value_policy<reference_existing_object>())
		.def("getBondIndex", &CAlphaRenderer::GetBondIndex)		
		.def("getBondCount", &CAlphaRenderer::GetBondCount)
		.def("deleteBond", &CAlphaRenderer::DeleteBond)		
		.def("addSideChainBond", &CAlphaRenderer::AddSideChainBond)
		.def("getSideChainBond", &CAlphaRenderer::GetSideChainBond, return_value_policy<reference_existing_object>())
		.def("getSideChainBondIndex", &CAlphaRenderer::GetSideChainBondIndex)		
		.def("getSideChainBondCount", &CAlphaRenderer::GetSideChainBondCount)
		.def("deleteSideChainBond", &CAlphaRenderer::DeleteSideChainBond)
		.def("getDeletedBonds", &CAlphaRenderer::getDeletedBonds)			
		.def("setSpacing", &CAlphaRenderer::SetSpacing)
		.def("getSpacingX", &CAlphaRenderer::GetSpacingX)
		.def("getSpacingY", &CAlphaRenderer::GetSpacingY)
		.def("getSpacingZ", &CAlphaRenderer::GetSpacingZ)
		.def("setOrigin", &CAlphaRenderer::SetOrigin)
		.def("getOriginX", &CAlphaRenderer::GetOriginX)
		.def("getOriginY", &CAlphaRenderer::GetOriginY)
		.def("getOriginZ", &CAlphaRenderer::GetOriginZ)
		.def("transformAllAtomLocations", &CAlphaRenderer::TransformAllAtomLocations)		
		.def("setDisplayStyle", &CAlphaRenderer::SetDisplayStyle)
		.def("setObjectSpecificColoring", &CAlphaRenderer::SetObjectSpecificColoring)
		.def("updateBoundingBox", &CAlphaRenderer::UpdateBoundingBox)
		.def("startHelix", &CAlphaRenderer::StartHelix)
		.def("addHelixElement", &CAlphaRenderer::AddHelixElement)
		.def("startStrand", &CAlphaRenderer::StartStrand)
		.def("addStrandElement", &CAlphaRenderer::AddStrandElement)
		.def("startLoop", &CAlphaRenderer::StartLoop)
		.def("addLoopElement", &CAlphaRenderer::AddLoopElement)
		.def("cleanSecondaryStructures", &CAlphaRenderer::CleanSecondaryStructures)
		.def("setNumSegments", &CAlphaRenderer::SetNumSegments)
		.def("setNumSlices", &CAlphaRenderer::SetNumSlices)
		.def("getSelectedHelixIndices", &CAlphaRenderer::GetSelectedHelixIndices)
		.def("setHelixCorrs", &CAlphaRenderer::SetHelixCorrs)
		.def("setSelectedSSEHelices", &CAlphaRenderer::SetSelectedSSEHelices)
		.def("clearOtherHighlights", &CAlphaRenderer::ClearOtherHighlights)
		.def("setFeatureVecs", &CAlphaRenderer::SetFeatureVecs)
		.def("setHelixColor", &CAlphaRenderer::SetHelixColor)
		.def("startHelix", &CAlphaRenderer::StartHelix)
		.def("addHelixElement", &CAlphaRenderer::AddHelixElement)
		.def("startStrand", &CAlphaRenderer::StartStrand)
		.def("addStrandElement", &CAlphaRenderer::AddStrandElement)
		.def("startLoop", &CAlphaRenderer::StartLoop)
		.def("addLoopElement", &CAlphaRenderer::AddLoopElement)
		.def("cleanSecondaryStructures", &CAlphaRenderer::CleanSecondaryStructures)
		.def("setNumSegments", &CAlphaRenderer::SetNumSegments)
		.def("setNumSlices", &CAlphaRenderer::SetNumSlices)
		.def("setHltRValue", &CAlphaRenderer::SetHltRValue)
		.def("setHltGValue", &CAlphaRenderer::SetHltGValue)
		.def("setHltBValue", &CAlphaRenderer::SetHltBValue)
		.def("setHltAValue", &CAlphaRenderer::SetHltAValue)
		.def("drawBackboneModelPathwalker", &CAlphaRenderer::DrawBackboneModelPathwalker)
		.def("drawDeletedBond", &CAlphaRenderer::DrawDeletedBond)
		.def("undrawDeletedBond", &CAlphaRenderer::UndrawDeletedBond)
		.def("drawAddedBond", &CAlphaRenderer::DrawAddedBond)
		.def("undrawAddedBond", &CAlphaRenderer::UndrawAddedBond)
		.def("addMaxBond", &CAlphaRenderer::AddMaxBond)
		.def("addMinBond", &CAlphaRenderer::AddMinBond)
		.def("addSaddleBond", &CAlphaRenderer::AddSaddleBond)
		.def("findDistance", &CAlphaRenderer::FindDistance)
		.def("deleteAtomFromVisualization", &CAlphaRenderer::DeleteAtomFromVisualization)
		.def("setMaxOn", &CAlphaRenderer::setMaxOn)
		.def("setMinOn", &CAlphaRenderer::setMinOn)
		.def("setSaddleOn", &CAlphaRenderer::setSaddleOn)
		.def("toggleMaxOn", &CAlphaRenderer::toggleMaxOn)
		.def("toggleMinOn", &CAlphaRenderer::toggleMinOn)
		.def("toggleSaddleOn", &CAlphaRenderer::toggleSaddleOn)
		.def("DrawSurface", &CAlphaRenderer::DrawSurface)
		.def("addQuad", &CAlphaRenderer::addQuad)
		.def("clearQuads", &CAlphaRenderer::clearQuads)
		.def("getQuadPts", &CAlphaRenderer::getQuadPts)
		.def("getQuadNormals", &CAlphaRenderer::getQuadNormals)
		.def("makeMax", &CAlphaRenderer::makeMax)
		.def("makeMin", &CAlphaRenderer::makeMin)
		.def("makeSaddle", &CAlphaRenderer::makeSaddle)
		.def("setExtremalMode", &CAlphaRenderer::SetExtremalMode)
		.def("setExtremalParams", &CAlphaRenderer::setExtremalParams)
		.def("setExtremalChecks", &CAlphaRenderer::setExtremalChecks)
		.def("setExtremalHide", &CAlphaRenderer::setExtremalHide)
		.def("setExtremalSurfaceHide", &CAlphaRenderer::setExtremalSurfaceHide)
		.def("setMinMaxIntensitiesEigenvalues", &CAlphaRenderer::setMinMaxIntensitiesEigenvalues)
		.def("addQuadSurface", &CAlphaRenderer::AddQuadSurface)
		.def("addQuadEigenvectors", &CAlphaRenderer::addQuadEigenvectors)
		.def("setMaxSurfaceOn", &CAlphaRenderer::setMaxSurfaceOn)
		.def("setMinSurfaceOn", &CAlphaRenderer::setMinSurfaceOn)
		.def("setMinDisplay", &CAlphaRenderer::setMinDisplay)
		.def("setMaxDisplay", &CAlphaRenderer::setMaxDisplay)
		.def("setSaddleDisplay", &CAlphaRenderer::setSaddleDisplay)
		.def("addVertexPos", &CAlphaRenderer::addVertexPos)
		.def("findHelices", &CAlphaRenderer::findHelices)
		.def("writeMaxCurveToPDB", &CAlphaRenderer::writeMaxCurveToPDB)
		.def("findMaxCurves", &CAlphaRenderer::findMaxCurves)
		.def("findMaxHelices", &CAlphaRenderer::findMaxHelices)
		.def("findMaxCurveStrand", &CAlphaRenderer::findMaxCurveStrand)
		.def("addCurveCoeff", &CAlphaRenderer::addCurveCoeff)
		.def("findCurveMinsMaxes", &CAlphaRenderer::findCurveMinsMaxes)
		.def("addMaxEigenVector", &CAlphaRenderer::addMaxEigenVector)
		.def("addMinEigenVector", &CAlphaRenderer::addMinEigenVector)
		.def("addSaddleEigenVector", &CAlphaRenderer::addSaddleEigenVector)
		.def("setMaxCurveEigenMaxesAndMins", &CAlphaRenderer::setMaxCurveEigenMaxesAndMins)
		.def("sortMaxCurveDimensions", &CAlphaRenderer::sortMaxCurveDimensions)
		.def("addMaxPtEigenVector", &CAlphaRenderer::addMaxPtEigenVector)
		.def("addMinPtEigenVector", &CAlphaRenderer::addMinPtEigenVector)
		.def("addSaddlePtEigenVector", &CAlphaRenderer::addSaddlePtEigenVector)
		.def("interpolateHelixPoints", &CAlphaRenderer::interpolateHelixPoints)
		.def("normalizedHelicesMinsAndMaxes", &CAlphaRenderer::normalizedHelicesMinsAndMaxes)
		.def("addOrd", &CAlphaRenderer::addOrd)
		.def("findEigenMinMaxes", &CAlphaRenderer::findEigenMinMaxes)
		.def("addQuadEigenvectors", &CAlphaRenderer::addQuadEigenvectors)
		.def("setHelixDisplay", &CAlphaRenderer::setHelixDisplay)
		.def("setHelixDisplayPoints", &CAlphaRenderer::setHelixDisplayPoints)
		.def("setHelixDisplayThreshold", &CAlphaRenderer::setHelixDisplayThreshold)
		.def("setEllipsoidDisplay", &CAlphaRenderer::setEllipsoidDisplay)
		.def("setScaleEllipsoid", &CAlphaRenderer::setScaleEllipsoid)
		.def("setHelixCoilHeight", &CAlphaRenderer::setHelixCoilHeight)
		.def("setSegmentThreshold", &CAlphaRenderer::setSegmentThreshold)
		.def("setOrthoDistMin", &CAlphaRenderer::setOrthoDistMin)
		.def("setOrthoDistMax", &CAlphaRenderer::setOrthoDistMax)
		.def("averageEllipsoidScale", &CAlphaRenderer::averageEllipsoidScale)
	;

	class_<InteractiveSkeletonEngine>("InteractiveSkeletonEngine", init<Volume *, NonManifoldMesh_Annotated *, float, int, int, int, unsigned int>())		
		.def("startEndPolyLineMode", &InteractiveSkeletonEngine::StartEndPolyLineMode)
		.def("startEndSingleRootMode", &InteractiveSkeletonEngine::StartEndSingleRootMode)
		.def("browseStartSeedRay", &InteractiveSkeletonEngine::BrowseStartSeedRay)
		.def("selectStartSeedRay", &InteractiveSkeletonEngine::SelectStartSeedRay)
		.def("selectEndSeed", &InteractiveSkeletonEngine::SelectEndSeed)
		.def("selectRootRay", &InteractiveSkeletonEngine::SelectRootRay)
		.def("analyzePathRay", &InteractiveSkeletonEngine::AnalyzePathRay)
		.def("setIsoValue", &InteractiveSkeletonEngine::SetIsoValue)		
		.def("clearSkeleton", &InteractiveSkeletonEngine::ClearSkeleton)		
		.def("clearCurrentPath", &InteractiveSkeletonEngine::ClearCurrentPath)		
		.def("finalizeSkeleton", &InteractiveSkeletonEngine::FinalizeSkeleton)		
		.def("draw", &InteractiveSkeletonEngine::Draw)		
		.def("clearSketch2D", &InteractiveSkeletonEngine::ClearSketch2D)		
		.def("setSketch2D", &InteractiveSkeletonEngine::SetSketch2D)		
		.def("setSketchRay", &InteractiveSkeletonEngine::SetSketchRay)		
		.def("clearSketchRay", &InteractiveSkeletonEngine::ClearSketchRay)
		.def("endSketchRay", &InteractiveSkeletonEngine::EndSketchRay)
		.def("addSelectionPoint", &InteractiveSkeletonEngine::AddSelectionPoint)
		.def("selectSelection", &InteractiveSkeletonEngine::SelectSelection)
		.def("deleteSelection", &InteractiveSkeletonEngine::DeleteSelection)
		.def("cancelSelection", &InteractiveSkeletonEngine::CancelSelection)
		.def("setLineThickness", &InteractiveSkeletonEngine::SetLineThickness)
	;

	class_<InteractiveLoopBuilderEngine>("InteractiveLoopBuilderEngine", init<Volume *, NonManifoldMesh_Annotated *, CAlphaRenderer *, float, int, int, int, unsigned int>())		
		.def("startEndPolyLineMode", &InteractiveLoopBuilderEngine::StartEndPolyLineMode)
		.def("startEndSingleRootMode", &InteractiveLoopBuilderEngine::StartEndSingleRootMode)
		.def("browseStartSeedRay", &InteractiveLoopBuilderEngine::BrowseStartSeedRay)
		.def("selectStartSeedRay", &InteractiveLoopBuilderEngine::SelectStartSeedRay)
		.def("selectEndSeed", &InteractiveLoopBuilderEngine::SelectEndSeed)
		.def("selectRootRay", &InteractiveLoopBuilderEngine::SelectRootRay)
		.def("analyzePathRay", &InteractiveLoopBuilderEngine::AnalyzePathRay)
		.def("setIsoValue", &InteractiveLoopBuilderEngine::SetIsoValue)		
		.def("clearSkeleton", &InteractiveLoopBuilderEngine::ClearSkeleton)		
		.def("clearCurrentPath", &InteractiveLoopBuilderEngine::ClearCurrentPath)		
		.def("finalizeSkeleton", &InteractiveLoopBuilderEngine::FinalizeSkeleton)		
		.def("draw", &InteractiveLoopBuilderEngine::Draw)		
		.def("clearSketch2D", &InteractiveLoopBuilderEngine::ClearSketch2D)		
		.def("setSketch2D", &InteractiveLoopBuilderEngine::SetSketch2D)		
		.def("setSketchRay", &InteractiveLoopBuilderEngine::SetSketchRay)		
		.def("clearSketchRay", &InteractiveLoopBuilderEngine::ClearSketchRay)
		.def("endSketchRay", &InteractiveLoopBuilderEngine::EndSketchRay)
		.def("addSelectionPoint", &InteractiveLoopBuilderEngine::AddSelectionPoint)
		.def("selectSelection", &InteractiveLoopBuilderEngine::SelectSelection)
		.def("deleteSelection", &InteractiveLoopBuilderEngine::DeleteSelection)
		.def("cancelSelection", &InteractiveLoopBuilderEngine::CancelSelection)
		.def("setLineThickness", &InteractiveLoopBuilderEngine::SetLineThickness)
		.def("addAtom", &InteractiveLoopBuilderEngine::AddAtom)
		.def("commitAtomList", &InteractiveLoopBuilderEngine::CommitAtomList)
		.def("clearAtomList", &InteractiveLoopBuilderEngine::ClearAtomList)
		.def("setTranslation", &InteractiveLoopBuilderEngine::SetTranslation)
		.def("setScaling", &InteractiveLoopBuilderEngine::SetScaling)
		.def("finishLoopBuilding", &InteractiveLoopBuilderEngine::FinishLoopBuilding)
		
	;

	class_<SSECorrespondenceResult>("SSECorrespondenceResult", init<>())
		.def("getNodeString", &SSECorrespondenceResult::GetNodeString)
		.def("getCost", &SSECorrespondenceResult::GetCost)
		.def("getNodeCount", &SSECorrespondenceResult::GetNodeCount)
		.def("getHelixCount", &SSECorrespondenceResult::GetHelixCount)
		.def("getSheetCount", &SSECorrespondenceResult::GetSheetCount)
		.def("getSkeletonNode", &SSECorrespondenceResult::GetSkeletonNode)
		.def("nodeToHelix", &SSECorrespondenceResult::NodeToHelix)		
	;

	class_<SecondaryStructure>("SecondaryStructure", init<>())
		.def("isHelix", &SecondaryStructure::IsHelix)
		.def("isSheet", &SecondaryStructure::IsSheet)
		.def("getSerialNumber", &SecondaryStructure::GetSerialNumber)
		.def("getStartPosition", &SecondaryStructure::GetStartPosition)
		.def("getEndPosition", &SecondaryStructure::GetEndPosition)
		.def("getSecondaryStructureID", &SecondaryStructure::GetSecondaryStructureID)
	;

	class_<GeometricShape>("GeometricShape", init<>())
		.def("isHelix", &GeometricShape::IsHelix)
		.def("isSheet", &GeometricShape::IsSheet)
		.def("getCornerCell2", &GeometricShape::GetCornerCell2)
		.def("getCornerCell3", &GeometricShape::GetCornerCell3)
		.def("setColor", &GeometricShape::SetColor)
		.def("getSelected", &GeometricShape::GetSelected)
	;

	bool (SSECorrespondenceEngine::*SetConstant1)(char *, char *)	= &SSECorrespondenceEngine::SetConstant;
	bool (SSECorrespondenceEngine::*SetConstant2)(char *, double)	= &SSECorrespondenceEngine::SetConstant;
	bool (SSECorrespondenceEngine::*SetConstant3)(char *, int)		= &SSECorrespondenceEngine::SetConstant;
	bool (SSECorrespondenceEngine::*SetConstant4)(char *, bool)		= &SSECorrespondenceEngine::SetConstant;

	class_<SSECorrespondenceEngine>("SSECorrespondenceEngine", init<>())
		.def("setConstant", SetConstant1)
		.def("setConstant", SetConstant2)
		.def("setConstantInt", SetConstant3)
		.def("setConstantBool", SetConstant4)
		.def("getConstant", &SSECorrespondenceEngine::GetConstant)
		.def("getConstantString", &SSECorrespondenceEngine::GetConstantString)
		.def("getConstantDouble", &SSECorrespondenceEngine::GetConstantDouble)
		.def("getConstantInt", &SSECorrespondenceEngine::GetConstantInt)
		.def("getConstantBool", &SSECorrespondenceEngine::GetConstantBool)
		.def("clearAllConstraints", &SSECorrespondenceEngine::ClearAllConstraints)			
		.def("setHelixConstraint", &SSECorrespondenceEngine::SetHelixConstraint)			
		.def("getHelixConstraintFwd", &SSECorrespondenceEngine::GetHelixConstraintFwd)			
		.def("getHelixConstraintRev", &SSECorrespondenceEngine::GetHelixConstraintRev)			
		.def("getHelixConstraintUnk", &SSECorrespondenceEngine::GetHelixConstraintUnk)
		.def("getStrandConstraint", &SSECorrespondenceEngine::GetStrandConstraint)			
		.def("setNodeConstraint", &SSECorrespondenceEngine::SetNodeConstraint)			
		.def("setSSEColor", &SSECorrespondenceEngine::SetSSEColor)			
        .def("setConstantsFromFile", &SSECorrespondenceEngine::SetConstantsFromFile)		
		.def("saveCorrespondenceToFile", &SSECorrespondenceEngine::SaveCorrespondenceToFile)
		.def("loadSkeletonGraph", &SSECorrespondenceEngine::LoadSkeletonGraph)
		.def("loadSequenceGraph", &SSECorrespondenceEngine::LoadSequenceGraph)
		.def("loadCorrespondenceFromFile", &SSECorrespondenceEngine::LoadCorrespondenceFromFile)
		.def("executeQuery", &SSECorrespondenceEngine::ExecuteQuery)
		.def("cleanupMemory", &SSECorrespondenceEngine::CleanupMemory)
		.def("getResult", &SSECorrespondenceEngine::GetResult)
		.def("getResultCount", &SSECorrespondenceEngine::GetResultCount)
		.def("getSupportedLoadFileFormats", &SSECorrespondenceEngine::GetSupportedLoadFileFormats)
		.def("getSupportedSaveFileFormats", &SSECorrespondenceEngine::GetSupportedSaveFileFormats)
		.def("getSkeletonSSE", &SSECorrespondenceEngine::GetSkeletonSSE, return_value_policy<reference_existing_object>())
		.def("getSequenceSSE", &SSECorrespondenceEngine::GetSequenceSSE, return_value_policy<reference_existing_object>())
		.def("getSkeletonSSECount", &SSECorrespondenceEngine::GetSkeletonSSECount)
		.def("getSequenceSSECount", &SSECorrespondenceEngine::GetSequenceSSECount)
		.def("setVisibleCorrespondence", &SSECorrespondenceEngine::SetVisibleCorrespondence)
		.def("draw", &SSECorrespondenceEngine::Draw)
		.def("drawAllPaths", &SSECorrespondenceEngine::DrawAllPaths)
		.def("initializePathFinder", &SSECorrespondenceEngine::InitializePathFinder)
		.def("initializePathHelix", &SSECorrespondenceEngine::InitializePathHelix)
		.def("getPathSpace", &SSECorrespondenceEngine::GetPathSpace)
		.def("clearPathSpace", &SSECorrespondenceEngine::ClearPathSpace)
		.def("getPathVertexCount", &SSECorrespondenceEngine::GetPathVertexCount)
		.def("getPathVertex", &SSECorrespondenceEngine::GetPathVertex)
		.def("getPathEdgeCount", &SSECorrespondenceEngine::GetPathEdgeCount)
		.def("getEdgeVertexIndex", &SSECorrespondenceEngine::GetEdgeVertexIndex)
		.def("clearPathFinder", &SSECorrespondenceEngine::ClearPathFinder)
	;

	class_<SEQFileData>("SeqFileData", init<>())
		.def("getStartResNo", &SEQFileData::GetStartResNo)
		.def("getSequenceString", &SEQFileData::GetSequenceString)
		.def("getStructureString", &SEQFileData::GetStructureString)
		.def("getNumberOfStructures", &SEQFileData::GetNumberOfStructures)
		.def("getStructure", &SEQFileData::GetStructure, return_value_policy<reference_existing_object>())
	;
	
	class_<SEQReader>("SeqReader", init<>())
		.def("loadFile", &SEQReader::ReadSeqFileData)
		.staticmethod("loadFile")
	;

	enum_<RadialProfileType>("RadialProfileType")
		.value("gaussian", GAUSSIAN)
		.value("gaussianDip", GAUSSIAN_DIP)
		.value("polynomial", POLYNOMIAL)
		.export_values()
	;
	
	class_<SSEHunter>("SSEHunter")
		.def("createPseudoAtoms", &SSEHunter::CreatePseudoAtoms)
		.def("getNumberOfPseudoAtoms", &SSEHunter::GetNumberOfPseudoAtoms)
		.def("getPseudoAtom", &SSEHunter::GetPseudoAtom, return_value_policy<reference_existing_object>())
		.def("setCorrelationScores", &SSEHunter::SetCorrelationScores)
		.def("setSkeletonScores", &SSEHunter::SetSkeletonScores)
//		.def("setGeometryScores", &SSEHunter::SetGeometryScores)
		.def("getAtomDistances", &SSEHunter::GetAtomDistances)
//		.def("getNeighborhoodVoxels", &SSEHunter::GetNeighborhoodVoxels)
		.def("getLocalDirectionalityScores", &SSEHunter::GetLocalDirectionalityScores)
	;
	
	class_<FlexibleFittingEngine>("FlexibleFittingEngine")
		.def("startPDBHelix", &FlexibleFittingEngine::StartPDBHelix)
		.def("addPDBAtomLocation", &FlexibleFittingEngine::AddPDBAtomLocation)
		.def("endPDBHelix", &FlexibleFittingEngine::EndPDBHelix)
		.def("addSSEHelix", &FlexibleFittingEngine::AddSSEHelix)		
		.def("startSearch", &FlexibleFittingEngine::StartSearch)		
		.def("getRigidTransform", &FlexibleFittingEngine::GetRigidTransform)		
		.def("getHelixFlexibleTransform", &FlexibleFittingEngine::GetHelixFlexibleTransform)		
		.def("getCorrespondenceCount", &FlexibleFittingEngine::GetCorrespondenceCount)		
		.def("getClusterCount", &FlexibleFittingEngine::GetClusterCount)		
		.def("getHelixCount", &FlexibleFittingEngine::GetHelixCount)		
		.def("getCorrespondenceNode", &FlexibleFittingEngine::GetCorrespondenceNode)
		.def("resetEngine", &FlexibleFittingEngine::ResetEngine)
		.def("Deform", &FlexibleFittingEngine::Deform)
		.def("addAtomLocation", &FlexibleFittingEngine::addAtomLocation)
		.def("addHardHandleLocation", &FlexibleFittingEngine::addHardHandleLocation)
		.def("addSoftHandleLocation", &FlexibleFittingEngine::addSoftHandleLocation)
		.def("returnFlattenedLocations", &FlexibleFittingEngine::returnFlattenedLocations)
		.def("prepareDeform", &FlexibleFittingEngine::prepareDeform)
		.def("setLaplacainW", &FlexibleFittingEngine::setLaplacainW)
		.def("setHardHandleW", &FlexibleFittingEngine::setHardHandleW)
		.def("setSoftHandleW", &FlexibleFittingEngine::setSoftHandleW)
		.def("addCorrespondence", &FlexibleFittingEngine::AddCorrespondence)
		.def("getCorrespondenceFlips", &FlexibleFittingEngine::GetCorrespondenceFlips)
		.def("getAllCorrespondencesFlat", &FlexibleFittingEngine::GetAllCorrespondencesFlat)
		.def("flipCorrespondencePair", &FlexibleFittingEngine::FlipCorrespondencePair)
		.def("saveCorrs", &FlexibleFittingEngine::SaveCorrs)
		.def("loadSavedCorrs", &FlexibleFittingEngine::LoadSavedCorrs)
		.def("getRigidTransform2", &FlexibleFittingEngine::GetRigidTransform2)
		.def("mergeClusters", &FlexibleFittingEngine::MergeClusters)
		.def("getAllCAlphaFeatureVecsFlat", &FlexibleFittingEngine::GetAllCAlphaFeatureVecsFlat)	
		.def("getIsForward", &FlexibleFittingEngine::GetIsForward)	
		.def("getPairRigidTransform", &FlexibleFittingEngine::GetPairRigidTransform)
		.def("getSideChainTransform", &FlexibleFittingEngine::GetSideChainTransform)
	;




	class_<SSECorrespondenceNode>("SSECorrespondenceNode")
		.def("getPIndex", &SSECorrespondenceNode::GetPIndex)
		.def("getQIndex", &SSECorrespondenceNode::GetQIndex)
		.def("isForward", &SSECorrespondenceNode::IsForward)
	;

    class_<LinearSolver>("LinearSolver")
        .def("findBestFitLine", &LinearSolver::FindBestFitLine)
        .staticmethod("findBestFitLine")
        .def("sumDistSqrd", LinearSolver::SumDistSqrd)
        .staticmethod("sumDistSqrd")
    ;
}



#endif
