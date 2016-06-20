/**
 * General data definition
 *
 * Author: Minxin Cheng
 * Date: 06/12/2012
 */


#pragma once

#include <cstdio>
#include <iostream>
#include <iomanip>
#include "Eigen/Dense"
#include <vector>
#include <list>
#include <cmath>
#include <ctime>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>
#include <fstream>
#include <string>
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "Segment.h"
#include "Quad.h"
#include "Point.h"
#include "Vertex.h"
#include "DispCell.h"
#include "Interpolation.h"
#include "View.h"
#include "mvc.h"
#include "ParallelCompute.h"
#include <tbb/task_scheduler_init.h>
#include "tbb/concurrent_vector.h"

using namespace std;
using namespace Eigen;

#define _USE_MATH_DEFINES

class General_Data
{
public:
  General_Data(View *_view3D);
  ~General_Data(void);
  virtual void dataGeneration(char *filename) {};
  virtual void buildGrid() {};
  void edgePhase(float* scalars,float* tensors,float* gradients, float *edgeTable);
  float *getGridShowPositions();
  Face *getFaces();
  float *getGridScalars();
  float *getGridV1s();
  float *getGridV2s();
  float *getGridV3s();
  float *getGridGradients();
  float *getEdgePoints();
  float *getFacePoints();
  float *getCellPoints();
  vector<Segment> *getSegments();
  vector<Quad> *getQuads();
  vector<Point> *getPoints();
  vector<Vertex> *getVertices();
  int getMaxGridIndex();
  float getMaxGradientNorm();
  float getMaxGrid();
  void facePhase(float* scalars,float* tensors,float* gradients, float *faceTable);
  void cellPhase(float* scalars,float* tensors,float* gradients);
  void edgePhase(float *coef);
  void facePhase(float *coef);
  void cellPhase(float *coef);
  void buildCurve();
  void buildSurface();
  void moveCell(int axis, bool inc);
  DispCell getDispCell();
  void getShowPos(float position[3], float showPos[3]);
  int getCellX();
  int getCellY();
  int getCellZ();
  void setCellX(int value);
  void setCellY(int value);
  void setCellZ(int value);
  void saveDisplay(bool resultDebug);
  void loadDisplay(bool resultDebug);
  void generateCubicVolume();
  void saveOFF(ofstream &ofs);
  void loadOFF(ifstream &ifs);
  void saveVPT(ofstream &ofs);
  void loadVPT(ifstream &ifs);
  int getFaceIndex(int axis, int x, int y, int z);

  char dataName[256];
  int *totalEdgePoints, *totalFacePoints, *totalCellPoints;
  float maxIntensity, minIntensity, maxScalar, minScalar, *maxEigenvalue, *minEigenvalue;
  int kernelSize, height, width, slices, dataType;
  float reverse, thickness, resolution;
  bool allCubic, resultDebug;
  float halfDataSize, rightBackTop[3], isovalue;
  int gridx, gridy, gridz;
  int sizex, sizey, sizez, halfSize;
  static const int subdNum = 3;

protected:
  virtual void getScalarGP(int index, float *scalar) {};
  virtual void getTensorGP(int index, float tensor[6]) {};
  virtual void getGradientGP(int index, float gradient[3]) {};
  virtual void getScalar(float position[3], float *scalar) {};
  virtual void getTensor(float position[3], float tensor[6]) {};
  virtual void getGradient(float position[3], float gradient[3]) {};
  virtual void getScalarCubic(float position[3], float *scalar) {};
  virtual void getTensorCubic(float position[3], float tensor[6]) {};
  virtual void getGradientCubic(float position[3], float gradient[3]) {};
  virtual void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6]) {};
  virtual void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3]) {};
  virtual void getGradientBicubicTable(int axis, int i, int j, int k, float sGradients[subdNum][subdNum][3]) {};

  int maxGrid, maxGridIndex;
  float *gridPoints;
  float gridSize;
  float maxGradientNorm;

private:
  bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver);
  bool getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver);
  bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver);
  bool getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver);
  bool getEigensolver2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2);
  bool getEigensolverCubic2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2);
  void computeGridShowPositions();
  void computeGridScalars();
  void computeGridVectors();
  void computeGridGradients();
  void computeEdgePoints();
  void computeFacePoints();
  void computeCellPoints();
  void getGridPointPos(int index, float v[3]);
  void getGridPointPosD(int index, float v[3]);
  void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
  void getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
  void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
  float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3]);
  void generateCubicMRC();
  void generateCubicDTI();


  void edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k);
  int getEdgeIndex(int axis, int x, int y, int z);
  void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC);
  float getDihedral(float b1[3], float b2[3], float b3[3]);


  bool adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
    float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG);
  /*bool adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
    float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG);*/
  bool adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
    float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth);
  /*bool adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
    float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth);*/
  void extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum);
  void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum);
  void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum);
  void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum);
  void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum);


  void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4,
    int index1, int index2, int index3, int index4, int axis, int i, int j, int k);
  void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int index1, int index2, int index3, int index4, int axis);
  void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h);
  void rotate2D(float p[2], float theta);
  void getCurveType(Face *face);
  void get3DWindingNum(Face *face, int ii, int jj, int kk, int axis);
  float signedArea3D(float v1[3], float v2[3], float v3[3]);


  void cellPhaseIteration(int i, int j, int k);
  bool getCentroid(int i, int j, int k, Cell *cell);
  void centroidProjection(int i, int j, int k, Cell *cell);
  void getProjDirection(int type, float position[3], float q[3]);
  float smallAbsA(float A);
  void buildCurveIteration(Face *face, Cell *cell1, Cell *cell2);
  void buildSurfaceIteration(Edge *edge, Cell *cell1, Cell *cell2, Cell *cell3, Cell *cell4);
  void buildPointIteration(Cell *cell);


  bool dispAdaptiveSampling(float p1[3], float p2[3], vector<float> *sampleP,
    vector<float> *sampleG, vector<float> *sampleV1, vector<float> *sampleV3);
  void dispAdaptiveRecursion(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3], float v3_1[3],
    float v3_2[3], vector<float> *PL, vector<float> *GL, vector<float> *V1L, vector<float> *V3L, int depth);
  void dispOrientV3(vector<float> *sampleV3);
  void dispPropagateXY(vector<float> *sampleV3, vector<float> *sampleX, vector<float> *sampleY);
  void dispCombineAB(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, float *sumAB);


  friend class boost::serialization::access;
  template <class Archive> void save(Archive &ar, const unsigned int version) const
  {
    if (resultDebug)
    {
      ar & (*totalEdgePoints);
      ar & (*totalFacePoints);
      ar & (*totalCellPoints);
      ar & maxGrid;
      ar & maxGridIndex;
      for ( int i = 0 ; i < 3 * gridx * gridy * gridz ; i ++ )
      {
        ar & gridPoints[i];
      }
      ar & gridSize;
      ar & maxGradientNorm;
      for ( int i = 0 ; i < maxGridIndex ; i ++ )
      {
        ar & gridScalars[i];
      }
      for ( int i = 0 ; i < 3 * maxGridIndex ; i ++ )
      {
        ar & gridShowPositions[i];
        ar & gridV1s[i];
        ar & gridV2s[i];
        ar & gridV3s[i];
        ar & gridGradients[i];
      }
      for ( int i = 0 ; i < 3 * (*totalEdgePoints) ; i ++ )
      {
        ar & edgePoints[i];
      }
      for ( int i = 0 ; i < 3 * (*totalFacePoints) ; i ++ )
      {
        ar & facePoints[i];
      }
      for ( int i = 0 ; i < 3 * (*totalCellPoints) ; i ++ )
      {
        ar & cellPoints[i];
      }
      ar & storeGridShowPos;
      ar & storeGridVector;
      ar & storeGridGradient;
      ar & storeGridScalar;
      ar & storeEdgePoint;
      ar & storeFacePoint;
      ar & storeCellPoint;
      ar & storeDispCell;
      ar & maxEdgeIndex;
      ar & maxFaceIndex;
      ar & maxCellIndex;
      for ( int i = 0 ; i < maxEdgeIndex ; i ++ )
      {
        ar & edges[i];
      }
      for ( int i = 0 ; i < maxFaceIndex ; i ++ )
      {
        ar & faces[i];
      }
      for ( int i = 0 ; i < maxCellIndex ; i ++ )
      {
        ar & cells[i];
      }
      ar & dispCellX;
      ar & dispCellY;
      ar & dispCellZ;
      ar & dispCell;
      ar & showGridScalar;
      ar & showGridV1;
      ar & showGridV2;
      ar & showGridV3;
      ar & showGridGradient;
      ar & showEdgePoint;
      ar & showFacePoint;
      ar & showCellPoint;
      ar & showCell;
      ar & showSample;
      ar & showFaceSampleGradient;
    }
    else
    {
      ar & dataName;
      ar & maxIntensity;
      ar & minIntensity;
      ar & maxScalar;
      ar & minScalar;
      ar & (*maxEigenvalue);
      ar & (*minEigenvalue);
      ar & kernelSize;
      ar & height;
      ar & width;
      ar & slices;
      ar & thickness;
      ar & resolution;
      ar & dataType;
      ar & reverse;
      ar & allCubic;
      ar & halfDataSize;
      ar & rightBackTop;
      ar & isovalue;
      ar & gridx;
      ar & gridy;
      ar & gridz;
      ar & segments;
      ar & quads;
      ar & points;
      ar & vertices;
      ar & pointRatio;
      ar & curveRatio;
      ar & surfaceRatio;
      ar & localIntensityThreshMinG;
      ar & localIntensityThreshMaxG;
      ar & eigenvalueThresh;
      ar & saddleCurve;
      ar & hideCurve;
      ar & hideSurface;
      ar & minCurve;
      ar & minSurface;
      ar & minPoint;
      ar & maxCurve;
      ar & maxSurface;
      ar & maxPoint;
      ar & showIsosurface;
      ar & isosurfaceFace;
      ar & cylinderCurve;
      ar & shadingFace;
      ar & shadingWireframe;
      ar & hideSaliencyRatio;
      ar & hideLocalIntensity;
      ar & hideEigenvalue;
      ar & saddlePoint;
      ar & hidePoint;
    }
  }
  
  template <class Archive> void load(Archive &ar, const unsigned int version)
  {
    if (resultDebug)
    {
      ar & (*totalEdgePoints);
      ar & (*totalFacePoints);
      ar & (*totalCellPoints);
      ar & maxGrid;
      ar & maxGridIndex;
      gridPoints = new float [ 3 * gridx * gridy * gridz ];
      for ( int i = 0 ; i < 3 * gridx * gridy * gridz ; i ++ )
      {
        ar & gridPoints[i];
      }
      ar & gridSize;
      ar & maxGradientNorm;
      gridScalars = new float [ maxGridIndex ];
      gridShowPositions = new float [ 3 * maxGridIndex ];
      gridV1s = new float [ 3 * maxGridIndex ];
      gridV2s = new float [ 3 * maxGridIndex ];
      gridV3s = new float [ 3 * maxGridIndex ];
      gridGradients = new float [ 3 * maxGridIndex ];
      for ( int i = 0 ; i < maxGridIndex ; i ++ )
      {
        ar & gridScalars[i];
      }
      for ( int i = 0 ; i < 3 * maxGridIndex ; i ++ )
      {
        ar & gridShowPositions[i];
        ar & gridV1s[i];
        ar & gridV2s[i];
        ar & gridV3s[i];
        ar & gridGradients[i];
      }
      edgePoints = new float [ 3 * (*totalEdgePoints) ];
      for ( int i = 0 ; i < 3 * (*totalEdgePoints) ; i ++ )
      {
        ar & edgePoints[i];
      }
      facePoints = new float [ 3 * (*totalFacePoints) ];
      for ( int i = 0 ; i < 3 * (*totalFacePoints) ; i ++ )
      {
        ar & facePoints[i];
      }
      cellPoints = new float [ 3 * (*totalCellPoints) ];
      for ( int i = 0 ; i < 3 * (*totalCellPoints) ; i ++ )
      {
        ar & cellPoints[i];
      }
      ar & storeGridShowPos;
      ar & storeGridVector;
      ar & storeGridGradient;
      ar & storeGridScalar;
      ar & storeEdgePoint;
      ar & storeFacePoint;
      ar & storeCellPoint;
      ar & storeDispCell;
      ar & maxEdgeIndex;
      ar & maxFaceIndex;
      ar & maxCellIndex;
      edges = new Edge [ maxEdgeIndex ];
      for ( int i = 0 ; i < maxEdgeIndex ; i ++ )
      {
        ar & edges[i];
      }
      faces = new Face [ maxFaceIndex ];
      for ( int i = 0 ; i < maxFaceIndex ; i ++ )
      {
        ar & faces[i];
      }
      cells = new Cell [ maxCellIndex ];
      for ( int i = 0 ; i < maxCellIndex ; i ++ )
      {
        ar & cells[i];
      }
      ar & dispCellX;
      ar & dispCellY;
      ar & dispCellZ;
      ar & dispCell;
      ar & showGridScalar;
      ar & showGridV1;
      ar & showGridV2;
      ar & showGridV3;
      ar & showGridGradient;
      ar & showEdgePoint;
      ar & showFacePoint;
      ar & showCellPoint;
      ar & showCell;
      ar & showSample;
      ar & showFaceSampleGradient;
    }
    else
    {
      ar & dataName;
      ar & maxIntensity;
      ar & minIntensity;
      ar & maxScalar;
      ar & minScalar;
      ar & (*maxEigenvalue);
      ar & (*minEigenvalue);
      ar & kernelSize;
      ar & height;
      ar & width;
      ar & slices;
      ar & thickness;
      ar & resolution;
      ar & dataType;
      ar & reverse;
      ar & allCubic;
      ar & halfDataSize;
      ar & rightBackTop;
      ar & isovalue;
      ar & gridx;
      ar & gridy;
      ar & gridz;
      ar & segments;
      ar & quads;
      ar & points;
      ar & vertices;
      ar & pointRatio;
      ar & curveRatio;
      ar & surfaceRatio;
      ar & localIntensityThreshMinG;
      ar & localIntensityThreshMaxG;
      ar & eigenvalueThresh;
      ar & saddleCurve;
      ar & hideCurve;
      ar & hideSurface;
      ar & minCurve;
      ar & minSurface;
      ar & minPoint;
      ar & maxCurve;
      ar & maxSurface;
      ar & maxPoint;
      ar & showIsosurface;
      ar & isosurfaceFace;
      ar & cylinderCurve;
      ar & shadingFace;
      ar & shadingWireframe;
      ar & hideSaliencyRatio;
      ar & hideLocalIntensity;
      ar & hideEigenvalue;
      ar & saddlePoint;
      ar & hidePoint;
    }
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  float *gridScalars;
  float *gridShowPositions, *gridV1s, *gridV2s, *gridV3s, *gridGradients, 
    *edgePoints, *facePoints, *cellPoints;
  bool storeGridShowPos, storeGridVector, storeGridGradient, storeGridScalar, 
    storeEdgePoint, storeFacePoint, storeCellPoint, storeDispCell;
  int maxEdgeIndex, maxFaceIndex, maxCellIndex;
  Edge *edges;     
  Face *faces;
  Cell *cells;
  float change;
  float globalVec[3];
  float *grid_2DGradMag;
  bool *storeGrid_2DGradMag;
  vector<Segment> *segments;
  vector<Quad> *quads;
  vector<Point> *points;
  vector<Vertex> *vertices;
  concurrent_vector<Point> *conpoints;
  concurrent_vector<Vertex> *convertices;
  float thresh;
  int dispCellX, dispCellY, dispCellZ;
  DispCell dispCell;
  /*clock_t t_sampling, t_sampling_interp, t_sampling_other, t_extrEdge, t_extrEdge_orient, t_extrEdge_interp, t_extrEdge_other,
    t_orientV3, t_propagate, t_propagate_getXY, t_propagate_getEdgeAB, t_projectG, t_getEdgeW,
    t_facePhase_other, t_facePhase_extremalPoint, t_cellPhase_other, t_cellPhase_extremalPoint;*/
  //clock_t t_sampling_interp, t_extrEdge_interp, t_faceSampling_interp, t_faceOther_interp, t_cellPhase_interp;

  View *view3D;
  float pointRatio, curveRatio, surfaceRatio, localIntensityThreshMinG, localIntensityThreshMaxG, eigenvalueThresh;
  bool showGridScalar, showGridV1, showGridV2, showGridV3, showGridGradient, showEdgePoint, showFacePoint, showCellPoint,
    saddleCurve, saddlePoint, hideCurve, hideSurface, hidePoint, minCurve, minSurface, minPoint,
    maxCurve, maxSurface, maxPoint, showCell, showSample, showFaceSampleGradient,
    showIsosurface, hideSaliencyRatio, hideLocalIntensity, hideEigenvalue,
    isosurfaceFace, cylinderCurve, shadingFace, shadingWireframe;
};
