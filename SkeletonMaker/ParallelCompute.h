#include <iostream>
#include <vector>
#include <list>
#include <ctime>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range3d.h"
#include "tbb/concurrent_vector.h"
#include <Eigen/Dense>

#include "Interpolation.h"
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "mvc.h"
#include "Point.h"
#include "Vertex.h"

using namespace std;
using namespace tbb;
using namespace Eigen;

class ParallelEdge
{
private:

	float thresh;
	int gridx,gridy,gridz;
	int *totalEdgePoints;
	int dataType;
	float* globalVec;
	Edge *edges;
	float *gridPoints;
	
	float change;
	int sizex,sizey,sizez, halfSize;
	float*scalars,*tensors,*gradients;
	float *edgeTable;


	void getGridPointPosD(int index, float v[3]) const;

	bool adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
		float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG, float *edgeTable) const;
	bool adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
		float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth, float *edgeTable) const;
	void extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints) const;
	void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum) const;
	void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum) const;
	void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum) const;
	void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum) const;

	void edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints, float *edgeTable) const;

	
	bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	bool getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	bool getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver, float *edgeTable) const;
	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;
	void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC) const;
	float getDihedral(float b1[3], float b2[3], float b3[3]) const;
	int getEdgeIndex(int axis, int x, int y, int z) const;

public: 
	bool allCubic;
	bool *storeGrid_2DGradMag;
	float *grid_2DGradMag;
	ParallelEdge(int lgridx,int lgridy,int lgridz,int *ltotalEdgePoints,Edge *ledges,bool lallCubic,float lglobalVec[3],int ldataType,float *lgridPoints,float lchange,bool *lstoreGrid_2DGradMag,float *lgrid_2DGradMag,int lsizex,int lsizey,int lsizez,int lhalfSize,float* lscalars,float* ltensors,float* lgradients,float *ledgeTable)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalEdgePoints=ltotalEdgePoints;
		edges = ledges;
		allCubic=lallCubic;
		globalVec = lglobalVec;
		dataType = ldataType;
		gridPoints = lgridPoints;
		change=lchange;
		storeGrid_2DGradMag = lstoreGrid_2DGradMag;
		grid_2DGradMag = lgrid_2DGradMag;
		thresh = cos(20 * (float)M_PI / 180);
		sizex = lsizex;
		sizey = lsizey;
		sizez = lsizez;
		halfSize = lhalfSize;
		scalars = lscalars;
		tensors = ltensors;
		gradients = lgradients;
		edgeTable = ledgeTable;
		

	};

	void getScalarGP(int index, float *scalar) const
	{
		*scalar = scalars[index];
	}

	void getGradientGP(int index, float gradient[3]) const
	{
		gradient[0] = gradients[index];
		gradient[1] = gradients[index+1];
		gradient[2] = gradients[index+2];
	}

	void getTensorGP(int index, float tensor[6]) const
	{
		tensor[0] = tensors[index];
		tensor[1] = tensors[index+1];
		tensor[2] = tensors[index+2];
		tensor[3] = tensors[index+3];
		tensor[4] = tensors[index+4];
		tensor[5] = tensors[index+5];
	}

	void getScalar(float position[3], float *scalar) const
	{
		trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradient(float position[3], float gradient[3]) const
	{
		trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensor(float position[3], float tensor[6]) const
	{
		trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
	};

	void getScalarCubic(float position[3], float *scalar) const
	{
		tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradientCubic(float position[3], float gradient[3]) const
	{
		tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensorCubic(float position[3], float tensor[6]) const
	{
		tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
	};

	void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3], float *edgeTable) const
	{
		tricubic_3f_table(sizex, sizey, sizez, edgeTable, axis, i+2, j+2, k+2, index, gradients, gradient);
	};

	void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6], float *edgeTable) const
	{
		tricubic_6f_table(sizex, sizey, sizez, edgeTable, axis, i+2, j+2, k+2, index, tensors, tensor);
	};
	// for tbb
	void operator()(const blocked_range3d<int>& r)	const {    

			int *ltotalEdgePoints= totalEdgePoints;			
			Edge *ledges = edges;

		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if (i+1==gridx)
						continue;
					int sIndex1 = getIndex(1, 0, i+2, j+2, k+2, sizex, sizey, sizez);
					int sIndex2 = getIndex(1, 0, i+3, j+2, k+2, sizex, sizey, sizez);
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, sIndex1, sIndex2, 1, i, j, k, ledges, ltotalEdgePoints, edgeTable);
				}
			}
		}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if(j+1==gridy)
						continue;
					int sIndex1 = getIndex(1, 0, i+2, j+2, k+2, sizex, sizey, sizez);
					int sIndex2 = getIndex(1, 0, i+2, j+3, k+2, sizex, sizey, sizez);
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, sIndex1, sIndex2, 2, i, j, k, ledges, ltotalEdgePoints, edgeTable);
				}
			}
		}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if(k+1==gridz)
						continue;
					int sIndex1 = getIndex(1, 0, i+2, j+2, k+2, sizex, sizey, sizez);
					int sIndex2 = getIndex(1, 0, i+2, j+2, k+3, sizex, sizey, sizez);
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, sIndex1, sIndex2, 3, i, j, k, ledges, ltotalEdgePoints, edgeTable);
				}
			}
		}

					
	};




};

class ParallelFace
{
	int gridx,gridy,gridz;
	Edge *edges;
	Face *faces;
	int *totalFacePoints;
	int dataType;
	float change;
	static const int subdNum = 3;

	float *gridPoints;

	int sizex,sizey,sizez, halfSize;
	float*scalars,*tensors,*gradients;
	float *faceTable;

	int getEdgeIndex(int axis, int x, int y, int z) const;
	int getFaceIndex(int axis, int x, int y, int z) const;
	void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int axis) const;
	void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h, int *ltotalFacePoints) const;
	void getCurveType(Face *face) const;
	void get3DWindingNum(Face *face, int ii, int jj, int kk, int axis, float *faceTable) const;
	float smallAbsA(float A) const;
	void getGridPointPos(int index, float v[3]) const;
	void getGridPointPosD(int index, float v[3]) const;
	void rotate2D(float p[2], float theta) const;
	bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	float signedArea3D(float v1[3], float v2[3], float v3[3]) const;
	float getDihedral(float b1[3], float b2[3], float b3[3]) const;
	void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4,int index1, int index2, int index3, int index4, int axis, int i, int j, int k,Face *lfaces,int *ltotalFacePoints, float *faceTable) const;

	void getScalarGP(int index, float *scalar) const
	{
		*scalar = scalars[index];
	}

	void getGradientGP(int index, float gradient[3]) const
	{
		gradient[0] = gradients[index];
		gradient[1] = gradients[index+1];
		gradient[2] = gradients[index+2];
	}

	void getTensorGP(int index, float tensor[6]) const
	{
		tensor[0] = tensors[index];
		tensor[1] = tensors[index+1];
		tensor[2] = tensors[index+2];
		tensor[3] = tensors[index+3];
		tensor[4] = tensors[index+4];
		tensor[5] = tensors[index+5];
	}

	void getScalar(float position[3], float *scalar) const
	{
		trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradient(float position[3], float gradient[3]) const
	{
		trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensor(float position[3], float tensor[6]) const
	{
		trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
	};

	void getScalarCubic(float position[3], float *scalar) const
	{
		tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradientCubic(float position[3], float gradient[3]) const
	{
		tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensorCubic(float position[3], float tensor[6]) const
	{
		tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
	};

	void getGradientBicubicTable(int axis, int ii, int jj, int kk, float sGradients[subdNum][subdNum][3], float *faceTable) const
	{
		float validData[16][3];
		int count = 0;
		switch(axis)
		{
		case 1:
			for ( int i = 0 ; i < 4 ; i ++ )
			{
				for ( int j = 0 ; j < 4 ; j ++ )
				{
					int idx = getIndex(3, 0, ii+2, jj+i+1, kk+j+1, sizex, sizey, sizez);
					validData[count][0] = gradients[idx];
					validData[count][1] = gradients[idx+1];
					validData[count][2] = gradients[idx+2];
					count++;
				}
			}
			break;
		case 2:
			for ( int i = 0 ; i < 4 ; i ++ )
			{
				for ( int j = 0 ; j < 4 ; j ++ )
				{
					int idx = getIndex(3, 0, ii+j+1, jj+2, kk+i+1, sizex, sizey, sizez);
					validData[count][0] = gradients[idx];
					validData[count][1] = gradients[idx+1];
					validData[count][2] = gradients[idx+2];
					count++;
				}
			}
			break;
		case 3:
			for ( int i = 0 ; i < 4 ; i ++ )
			{
				for ( int j = 0 ; j < 4 ; j ++ )
				{
					int idx = getIndex(3, 0, ii+i+1, jj+j+1, kk+2, sizex, sizey, sizez);
					validData[count][0] = gradients[idx];
					validData[count][1] = gradients[idx+1];
					validData[count][2] = gradients[idx+2];
					count++;
				}
			}
			break;
		}
		count = 0;
		for ( int si = 0 ; si < subdNum ; si ++ )
		{
			for ( int sj = 0 ; sj < subdNum ; sj ++ )
			{
				int tableIndex = 16 * count;
				count++;
				sGradients[si][sj][0] = 0;
				sGradients[si][sj][1] = 0;
				sGradients[si][sj][2] = 0;
				for ( int i = 0 ; i < 16 ; i ++ )
				{
					for ( int j = 0 ; j < 3 ; j ++ )
					{
						sGradients[si][sj][j] += faceTable[tableIndex + i] * validData[i][j];
					}
				}
			}
		}
	}

public: 
	bool *storeGrid_2DGradMag;
	float *grid_2DGradMag;

	ParallelFace(int lgridx,int lgridy,int lgridz,int *ltotalFacePoints,Edge *ledges,Face *lfaces,int ldataType,float *lgridPoints,float lchange,bool *lstoreGrid_2DGradMag,float *lgrid_2DGradMag,int lsizex,int lsizey,int lsizez,int lhalfSize,float* lscalars,float* ltensors,float* lgradients, float *lfaceTable)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalFacePoints=ltotalFacePoints;
		edges = ledges;
		faces  = lfaces;
		dataType = ldataType;
		gridPoints = lgridPoints;
		change=lchange;
		storeGrid_2DGradMag = lstoreGrid_2DGradMag;
		grid_2DGradMag = lgrid_2DGradMag;
		sizex = lsizex;
		sizey = lsizey;
		sizez = lsizez;
		halfSize = lhalfSize;
		scalars = lscalars;
		tensors = ltensors;
		gradients = lgradients;
		faceTable = lfaceTable;

	};
	void operator()(const blocked_range3d<int>& r)	const {
		int *ltotalFacePoints= totalFacePoints;			
		//Edge *ledges = edges;
		Face *lfaces = faces;

		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( j == gridy - 1)
						continue;
					if ( k == gridz - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i, j + 1, k + 1, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);

					Edge *edge1 = &(edges[getEdgeIndex(2, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(3, i, j + 1, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(2, i, j, k + 1)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
					if (!edge4->valid) continue;


					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 1, i, j, k,lfaces,ltotalFacePoints, faceTable);
				}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( i == gridx - 1)
						continue;
					if ( k == gridz - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i + 1, j, k + 1, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);


					Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(3, i + 1, j, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(1, i, j, k + 1)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
					if (!edge4->valid) continue;


					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 2, i, j, k,lfaces,ltotalFacePoints, faceTable);
				}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( i == gridx - 1)
						continue;
					if ( j == gridy - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i + 1, j + 1, k, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);

					Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(2, i + 1, j, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(1, i, j + 1, k)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(2, i, j, k)]);
					if (!edge4->valid) continue;
					

					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 3, i, j, k,lfaces,ltotalFacePoints, faceTable);
				}
	};
};


class ParallelCell
{


	int gridx,gridy,gridz;
	Edge *edges;
	Face *faces;
	Cell *cells;
	vector<Vertex> *vertices;
	vector<Point> *points;
	concurrent_vector<Point> *conpoints;
	concurrent_vector<Vertex> *convertices;

	float *maxEigenvalue, *minEigenvalue;
	float halfDataSize, *rightBackTop;
	float thickness;
	int *totalCellPoints;
	float gridSize;


	int dataType;
	float *gridPoints;
	float change;
	int sizex,sizey,sizez,halfSize;
	float*scalars,*tensors,*gradients;

	int getEdgeIndex(int axis, int x, int y, int z) const;
	int getFaceIndex(int axis, int x, int y, int z) const;

	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;
	void buildPointIteration(Cell *cell,concurrent_vector<Point> *points,concurrent_vector<Vertex> *vertices,float* maxEigenvalue, float* minEigenvalue) const;
	void getGridPointPosD(int index, float v[3]) const;
	void centroidProjection(int i, int j, int k, Cell *cell) const;
	bool getCentroid(int i, int j, int k, Cell *cell,int *ltotalCellPoints) const;
	void getProjDirection(int type, float position[3], float q[3]) const;
	bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3],float* maxEigenvalue, float* minEigenvalue) const;

	void cellPhaseIteration(int i, int j, int k, Cell *cells,concurrent_vector<Point> *lpoints,concurrent_vector<Vertex> *lvertices, int *ltotalCellPoints,float* lmaxEigenvalue, float* lminEigenvalue) const;

	void getShowPos(float position[3], float showPos[3]) const;

	void getScalarGP(int index, float *scalar) const
	{
		*scalar = scalars[index];
	}

	void getGradientGP(int index, float gradient[3]) const
	{
		gradient[0] = gradients[index];
		gradient[1] = gradients[index+1];
		gradient[2] = gradients[index+2];
	}

	void getTensorGP(int index, float tensor[6]) const
	{
		tensor[0] = tensors[index];
		tensor[1] = tensors[index+1];
		tensor[2] = tensors[index+2];
		tensor[3] = tensors[index+3];
		tensor[4] = tensors[index+4];
		tensor[5] = tensors[index+5];
	}

	void getScalar(float position[3], float *scalar) const
	{
		trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradient(float position[3], float gradient[3]) const
	{
		trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensor(float position[3], float tensor[6]) const
	{
		trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
	};

	void getScalarCubic(float position[3], float *scalar) const
	{
		tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
	};

	void getGradientCubic(float position[3], float gradient[3]) const
	{
		tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
	};

	void getTensorCubic(float position[3], float tensor[6]) const
	{
		tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
	};
public:

	ParallelCell(int lgridx,int lgridy,int lgridz,int *ltotalCellPoints,Edge *ledges,Face *lfaces,Cell *lcells,	concurrent_vector<Point> *lpoints,int ldataType,float *lgridPoints,float lchange,int lsizex,int lsizey,int lsizez,int lhalfSize,float* lscalars,float* ltensors,float* lgradients,float lgridSize,float* lmaxEigenvalue, float* lminEigenvalue,concurrent_vector<Vertex> *lvertices,float lhalfDataSize,float *lrightBackTop,	float lthickness)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalCellPoints=ltotalCellPoints;
		edges = ledges;
		faces  = lfaces;
		cells = lcells;
		conpoints = lpoints;
		convertices = lvertices;
		dataType = ldataType;
		gridPoints = lgridPoints;


		change=lchange;
		sizex = lsizex;
		sizey = lsizey;
		sizez = lsizez;
		halfSize = lhalfSize;
		scalars = lscalars;
		tensors = ltensors;
		gradients = lgradients;
		gridSize = lgridSize;
		//
		maxEigenvalue=lmaxEigenvalue;
		minEigenvalue=lminEigenvalue;

		halfDataSize = lhalfDataSize;
		rightBackTop = lrightBackTop;
		thickness = lthickness;

	};

	void operator()(const blocked_range3d<int>& r)	const {
		int *ltotalCellPoints= totalCellPoints;			
		//Edge *ledges = edges;
		//Face *lfaces = faces;
		Cell *lcells = cells;
		concurrent_vector<Point> *lpoints = conpoints; 
		float *lmaxEigenvalue = maxEigenvalue;
		float *lminEigenvalue = minEigenvalue;
		concurrent_vector<Vertex> *lvertices = convertices;

	for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
		for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
			for ( int k=r.pages().begin();k!=r.pages().end();++k )
			{
				if(i == gridx - 1 || j == gridy - 1 || k == gridz -1)
					continue;
				cellPhaseIteration(i, j, k,lcells,lpoints,lvertices,ltotalCellPoints,lmaxEigenvalue,lminEigenvalue);
			}
	};

};

class CGNSParallelEdge
{
private:
	float thresh;
	int gridx,gridy,gridz;
	int *totalEdgePoints;
	int dataType;
	float* globalVec;
	Edge *edges;
	float *gridPoints;
	float change;
	float *coef;

	void getGridPointPosD(int index, float v[3]) const;
	void edgePhaseIteration(int index1, int index2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints) const;
	bool adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
		float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG) const;
	bool adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
		float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth) const;
	void extremalEdgeArray(float p1[3], float p2[3], Edge *edge, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints) const;
	void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum) const;
	void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum) const;
	void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum) const;
	void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum) const;
	bool getEigensolverPN(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;
	void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC) const;
	float getDihedral(float b1[3], float b2[3], float b3[3]) const;
	int getEdgeIndex(int axis, int x, int y, int z) const;

public: 
	bool allCubic;
	bool *storeGrid_2DGradMag;
	float *grid_2DGradMag;

	CGNSParallelEdge(int lgridx,int lgridy,int lgridz,int *ltotalEdgePoints,Edge *ledges,bool lallCubic,float lglobalVec[3],
		int ldataType,float *lgridPoints,float lchange,bool *lstoreGrid_2DGradMag,float *lgrid_2DGradMag, float *lcoef)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalEdgePoints=ltotalEdgePoints;
		edges = ledges;
		allCubic=lallCubic;
		globalVec = lglobalVec;
		dataType = ldataType;
		gridPoints = lgridPoints;
		change=lchange;
		storeGrid_2DGradMag = lstoreGrid_2DGradMag;
		grid_2DGradMag = lgrid_2DGradMag;
		thresh = cos(20 * (float)M_PI / 180);
		coef = lcoef;
	};

	void getScalarPN(float position[3], float *scalar) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		*scalar = coef[56]*pow(x,6) + coef[63]*pow(x,5)*y + coef[57]*pow(x,5)*z + coef[35]*pow(x,5) + coef[69]*pow(x,4)*pow(y,2) + coef[64]*pow(x,4)*y*z +
			coef[41]*pow(x,4)*y + coef[58]*pow(x,4)*pow(z,2) + coef[36]*pow(x,4)*z + coef[20]*pow(x,4) + coef[74]*pow(x,3)*pow(y,3) + coef[70]*pow(x,3)*pow(y,2)*z +
			coef[46]*pow(x,3)*pow(y,2) + coef[65]*pow(x,3)*y*pow(z,2) + coef[42]*pow(x,3)*y*z + coef[25]*pow(x,3)*y + coef[59]*pow(x,3)*pow(z,3) + coef[37]*pow(x,3)*pow(z,2) +
			coef[21]*pow(x,3)*z + coef[10]*pow(x,3) + coef[78]*pow(x,2)*pow(y,4) + coef[75]*pow(x,2)*pow(y,3)*z + coef[50]*pow(x,2)*pow(y,3) + coef[71]*pow(x,2)*pow(y,2)*pow(z,2) +
			coef[47]*pow(x,2)*pow(y,2)*z + coef[29]*pow(x,2)*pow(y,2) + coef[66]*pow(x,2)*y*pow(z,3) + coef[43]*pow(x,2)*y*pow(z,2) + coef[26]*pow(x,2)*y*z +
			coef[14]*pow(x,2)*y + coef[60]*pow(x,2)*pow(z,4) + coef[38]*pow(x,2)*pow(z,3) + coef[22]*pow(x,2)*pow(z,2) + coef[11]*pow(x,2)*z + coef[4]*pow(x,2) +
			coef[81]*x*pow(y,5) + coef[79]*x*pow(y,4)*z + coef[53]*x*pow(y,4) + coef[76]*x*pow(y,3)*pow(z,2) + coef[51]*x*pow(y,3)*z + coef[32]*x*pow(y,3) +
			coef[72]*x*pow(y,2)*pow(z,3) + coef[48]*x*pow(y,2)*pow(z,2) + coef[30]*x*pow(y,2)*z + coef[17]*x*pow(y,2) + coef[67]*x*y*pow(z,4) + coef[44]*x*y*pow(z,3) +
			coef[27]*x*y*pow(z,2) + coef[15]*x*y*z + coef[7]*x*y + coef[61]*x*pow(z,5) + coef[39]*x*pow(z,4) + coef[23]*x*pow(z,3) +
			coef[12]*x*pow(z,2) + coef[5]*x*z + coef[1]*x + coef[83]*pow(y,6) + coef[82]*pow(y,5)*z + coef[55]*pow(y,5) + coef[80]*pow(y,4)*pow(z,2) +
			coef[54]*pow(y,4)*z + coef[34]*pow(y,4) + coef[77]*pow(y,3)*pow(z,3) + coef[52]*pow(y,3)*pow(z,2) + coef[33]*pow(y,3)*z + coef[19]*pow(y,3) +
			coef[73]*pow(y,2)*pow(z,4) + coef[49]*pow(y,2)*pow(z,3) + coef[31]*pow(y,2)*pow(z,2) + coef[18]*pow(y,2)*z + coef[9]*pow(y,2) + coef[68]*y*pow(z,5) +
			coef[45]*y*pow(z,4) + coef[28]*y*pow(z,3) + coef[16]*y*pow(z,2) + coef[8]*y*z + coef[3]*y + coef[62]*pow(z,6) + coef[40]*pow(z,5) +
			coef[24]*pow(z,4) + coef[13]*pow(z,3) + coef[6]*pow(z,2) + coef[2]*z + coef[0];
	};

	void getGradientPN(float position[3], float gradient[3]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		gradient[0] = 6*coef[56]*pow(x,5) + 5*coef[63]*pow(x,4)*y + 5*coef[57]*pow(x,4)*z + 5*coef[35]*pow(x,4) + 4*coef[69]*pow(x,3)*pow(y,2) +
			4*coef[64]*pow(x,3)*y*z + 4*coef[41]*pow(x,3)*y + 4*coef[58]*pow(x,3)*pow(z,2) + 4*coef[36]*pow(x,3)*z + 4*coef[20]*pow(x,3) + 3*coef[74]*pow(x,2)*pow(y,3) +
			3*coef[70]*pow(x,2)*pow(y,2)*z + 3*coef[46]*pow(x,2)*pow(y,2) + 3*coef[65]*pow(x,2)*y*pow(z,2) + 3*coef[42]*pow(x,2)*y*z + 3*coef[25]*pow(x,2)*y +
			3*coef[59]*pow(x,2)*pow(z,3) + 3*coef[37]*pow(x,2)*pow(z,2) + 3*coef[21]*pow(x,2)*z + 3*coef[10]*pow(x,2) + 2*coef[78]*x*pow(y,4) + 2*coef[75]*x*pow(y,3)*z +
			2*coef[50]*x*pow(y,3) + 2*coef[71]*x*pow(y,2)*pow(z,2) + 2*coef[47]*x*pow(y,2)*z + 2*coef[29]*x*pow(y,2) + 2*coef[66]*x*y*pow(z,3) +
			2*coef[43]*x*y*pow(z,2) + 2*coef[26]*x*y*z + 2*coef[14]*x*y + 2*coef[60]*x*pow(z,4) + 2*coef[38]*x*pow(z,3) + 2*coef[22]*x*pow(z,2) +
			2*coef[11]*x*z + 2*coef[4]*x + coef[81]*pow(y,5) + coef[79]*pow(y,4)*z + coef[53]*pow(y,4) + coef[76]*pow(y,3)*pow(z,2) + coef[51]*pow(y,3)*z +
			coef[32]*pow(y,3) + coef[72]*pow(y,2)*pow(z,3) + coef[48]*pow(y,2)*pow(z,2) + coef[30]*pow(y,2)*z + coef[17]*pow(y,2) + coef[67]*y*pow(z,4) + coef[44]*y*pow(z,3) +
			coef[27]*y*pow(z,2) + coef[15]*y*z + coef[7]*y + coef[61]*pow(z,5) + coef[39]*pow(z,4) + coef[23]*pow(z,3) + coef[12]*pow(z,2) + coef[5]*z + coef[1];
		gradient[1] = coef[63]*pow(x,5) + 2*coef[69]*pow(x,4)*y + coef[64]*pow(x,4)*z + coef[41]*pow(x,4) + 3*coef[74]*pow(x,3)*pow(y,2) + 2*coef[70]*pow(x,3)*y*z +
			2*coef[46]*pow(x,3)*y + coef[65]*pow(x,3)*pow(z,2) + coef[42]*pow(x,3)*z + coef[25]*pow(x,3) + 4*coef[78]*pow(x,2)*pow(y,3) + 3*coef[75]*pow(x,2)*pow(y,2)*z +
			3*coef[50]*pow(x,2)*pow(y,2) + 2*coef[71]*pow(x,2)*y*pow(z,2) + 2*coef[47]*pow(x,2)*y*z + 2*coef[29]*pow(x,2)*y + coef[66]*pow(x,2)*pow(z,3) + coef[43]*pow(x,2)*pow(z,2) +
			coef[26]*pow(x,2)*z + coef[14]*pow(x,2) + 5*coef[81]*x*pow(y,4) + 4*coef[79]*x*pow(y,3)*z + 4*coef[53]*x*pow(y,3) + 3*coef[76]*x*pow(y,2)*pow(z,2) +
			3*coef[51]*x*pow(y,2)*z + 3*coef[32]*x*pow(y,2) + 2*coef[72]*x*y*pow(z,3) + 2*coef[48]*x*y*pow(z,2) + 2*coef[30]*x*y*z + 2*coef[17]*x*y +
			coef[67]*x*pow(z,4) + coef[44]*x*pow(z,3) + coef[27]*x*pow(z,2) + coef[15]*x*z + coef[7]*x + 6*coef[83]*pow(y,5) + 5*coef[82]*pow(y,4)*z +
			5*coef[55]*pow(y,4) + 4*coef[80]*pow(y,3)*pow(z,2) + 4*coef[54]*pow(y,3)*z + 4*coef[34]*pow(y,3) + 3*coef[77]*pow(y,2)*pow(z,3) + 3*coef[52]*pow(y,2)*pow(z,2) +
			3*coef[33]*pow(y,2)*z + 3*coef[19]*pow(y,2) + 2*coef[73]*y*pow(z,4) + 2*coef[49]*y*pow(z,3) + 2*coef[31]*y*pow(z,2) + 2*coef[18]*y*z +
			2*coef[9]*y + coef[68]*pow(z,5) + coef[45]*pow(z,4) + coef[28]*pow(z,3) + coef[16]*pow(z,2) + coef[8]*z + coef[3];
		gradient[2] = coef[57]*pow(x,5) + coef[64]*pow(x,4)*y + 2*coef[58]*pow(x,4)*z + coef[36]*pow(x,4) + coef[70]*pow(x,3)*pow(y,2) + 2*coef[65]*pow(x,3)*y*z +
			coef[42]*pow(x,3)*y + 3*coef[59]*pow(x,3)*pow(z,2) + 2*coef[37]*pow(x,3)*z + coef[21]*pow(x,3) + coef[75]*pow(x,2)*pow(y,3) + 2*coef[71]*pow(x,2)*pow(y,2)*z +
			coef[47]*pow(x,2)*pow(y,2) + 3*coef[66]*pow(x,2)*y*pow(z,2) + 2*coef[43]*pow(x,2)*y*z + coef[26]*pow(x,2)*y + 4*coef[60]*pow(x,2)*pow(z,3) + 3*coef[38]*pow(x,2)*pow(z,2) +
			2*coef[22]*pow(x,2)*z + coef[11]*pow(x,2) + coef[79]*x*pow(y,4) + 2*coef[76]*x*pow(y,3)*z + coef[51]*x*pow(y,3) + 3*coef[72]*x*pow(y,2)*pow(z,2) +
			2*coef[48]*x*pow(y,2)*z + coef[30]*x*pow(y,2) + 4*coef[67]*x*y*pow(z,3) + 3*coef[44]*x*y*pow(z,2) + 2*coef[27]*x*y*z + coef[15]*x*y +
			5*coef[61]*x*pow(z,4) + 4*coef[39]*x*pow(z,3) + 3*coef[23]*x*pow(z,2) + 2*coef[12]*x*z + coef[5]*x + coef[82]*pow(y,5) + 2*coef[80]*pow(y,4)*z +
			coef[54]*pow(y,4) + 3*coef[77]*pow(y,3)*pow(z,2) + 2*coef[52]*pow(y,3)*z + coef[33]*pow(y,3) + 4*coef[73]*pow(y,2)*pow(z,3) + 3*coef[49]*pow(y,2)*pow(z,2) +
			2*coef[31]*pow(y,2)*z + coef[18]*pow(y,2) + 5*coef[68]*y*pow(z,4) + 4*coef[45]*y*pow(z,3) + 3*coef[28]*y*pow(z,2) + 2*coef[16]*y*z + coef[8]*y +
			6*coef[62]*pow(z,5) + 5*coef[40]*pow(z,4) + 4*coef[24]*pow(z,3) + 3*coef[13]*pow(z,2) + 2*coef[6]*z + coef[2];
	};

	void getTensorPN(float position[3], float tensor[6]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		tensor[0] = 30*coef[56]*pow(x,4) + 20*coef[63]*pow(x,3)*y + 20*coef[57]*pow(x,3)*z + 20*coef[35]*pow(x,3) + 12*coef[69]*pow(x,2)*pow(y,2) +
			12*coef[64]*pow(x,2)*y*z + 12*coef[41]*pow(x,2)*y + 12*coef[58]*pow(x,2)*pow(z,2) + 12*coef[36]*pow(x,2)*z + 12*coef[20]*pow(x,2) + 6*coef[74]*x*pow(y,3) +
			6*coef[70]*x*pow(y,2)*z + 6*coef[46]*x*pow(y,2) + 6*coef[65]*x*y*pow(z,2) + 6*coef[42]*x*y*z + 6*coef[25]*x*y + 6*coef[59]*x*pow(z,3) +
			6*coef[37]*x*pow(z,2) + 6*coef[21]*x*z + 6*coef[10]*x + 2*coef[78]*pow(y,4) + 2*coef[75]*pow(y,3)*z + 2*coef[50]*pow(y,3) +
			2*coef[71]*pow(y,2)*pow(z,2) + 2*coef[47]*pow(y,2)*z + 2*coef[29]*pow(y,2) + 2*coef[66]*y*pow(z,3) + 2*coef[43]*y*pow(z,2) + 2*coef[26]*y*z +
			2*coef[14]*y + 2*coef[60]*pow(z,4) + 2*coef[38]*pow(z,3) + 2*coef[22]*pow(z,2) + 2*coef[11]*z + 2*coef[4];
		tensor[1] = 5*coef[63]*pow(x,4) + 8*coef[69]*pow(x,3)*y + 4*coef[64]*pow(x,3)*z + 4*coef[41]*pow(x,3) + 9*coef[74]*pow(x,2)*pow(y,2) +
			6*coef[70]*pow(x,2)*y*z + 6*coef[46]*pow(x,2)*y + 3*coef[65]*pow(x,2)*pow(z,2) + 3*coef[42]*pow(x,2)*z + 3*coef[25]*pow(x,2) + 8*coef[78]*x*pow(y,3) +
			6*coef[75]*x*pow(y,2)*z + 6*coef[50]*x*pow(y,2) + 4*coef[71]*x*y*pow(z,2) + 4*coef[47]*x*y*z + 4*coef[29]*x*y + 2*coef[66]*x*pow(z,3) +
			2*coef[43]*x*pow(z,2) + 2*coef[26]*x*z + 2*coef[14]*x + 5*coef[81]*pow(y,4) + 4*coef[79]*pow(y,3)*z + 4*coef[53]*pow(y,3) +
			3*coef[76]*pow(y,2)*pow(z,2) + 3*coef[51]*pow(y,2)*z + 3*coef[32]*pow(y,2) + 2*coef[72]*y*pow(z,3) + 2*coef[48]*y*pow(z,2) + 2*coef[30]*y*z +
			2*coef[17]*y + coef[67]*pow(z,4) + coef[44]*pow(z,3) + coef[27]*pow(z,2) + coef[15]*z + coef[7];
		tensor[2] = 5*coef[57]*pow(x,4) + 4*coef[64]*pow(x,3)*y + 8*coef[58]*pow(x,3)*z + 4*coef[36]*pow(x,3) + 3*coef[70]*pow(x,2)*pow(y,2) +
			6*coef[65]*pow(x,2)*y*z + 3*coef[42]*pow(x,2)*y + 9*coef[59]*pow(x,2)*pow(z,2) + 6*coef[37]*pow(x,2)*z + 3*coef[21]*pow(x,2) + 2*coef[75]*x*pow(y,3) +
			4*coef[71]*x*pow(y,2)*z + 2*coef[47]*x*pow(y,2) + 6*coef[66]*x*y*pow(z,2) + 4*coef[43]*x*y*z + 2*coef[26]*x*y + 8*coef[60]*x*pow(z,3) +
			6*coef[38]*x*pow(z,2) + 4*coef[22]*x*z + 2*coef[11]*x + coef[79]*pow(y,4) + 2*coef[76]*pow(y,3)*z + coef[51]*pow(y,3) +
			3*coef[72]*pow(y,2)*pow(z,2) + 2*coef[48]*pow(y,2)*z + coef[30]*pow(y,2) + 4*coef[67]*y*pow(z,3) + 3*coef[44]*y*pow(z,2) + 2*coef[27]*y*z +
			coef[15]*y + 5*coef[61]*pow(z,4) + 4*coef[39]*pow(z,3) + 3*coef[23]*pow(z,2) + 2*coef[12]*z + coef[5];
		tensor[3] = 2*coef[69]*pow(x,4) + 6*coef[74]*pow(x,3)*y + 2*coef[70]*pow(x,3)*z + 2*coef[46]*pow(x,3) + 12*coef[78]*pow(x,2)*pow(y,2) +
			6*coef[75]*pow(x,2)*y*z + 6*coef[50]*pow(x,2)*y + 2*coef[71]*pow(x,2)*pow(z,2) + 2*coef[47]*pow(x,2)*z + 2*coef[29]*pow(x,2) + 20*coef[81]*x*pow(y,3) +
			12*coef[79]*x*pow(y,2)*z + 12*coef[53]*x*pow(y,2) + 6*coef[76]*x*y*pow(z,2) + 6*coef[51]*x*y*z + 6*coef[32]*x*y + 2*coef[72]*x*pow(z,3) +
			2*coef[48]*x*pow(z,2) + 2*coef[30]*x*z + 2*coef[17]*x + 30*coef[83]*pow(y,4) + 20*coef[82]*pow(y,3)*z + 20*coef[55]*pow(y,3) +
			12*coef[80]*pow(y,2)*pow(z,2) + 12*coef[54]*pow(y,2)*z + 12*coef[34]*pow(y,2) + 6*coef[77]*y*pow(z,3) + 6*coef[52]*y*pow(z,2) + 6*coef[33]*y*z +
			6*coef[19]*y + 2*coef[73]*pow(z,4) + 2*coef[49]*pow(z,3) + 2*coef[31]*pow(z,2) + 2*coef[18]*z + 2*coef[9];
		tensor[4] = coef[64]*pow(x,4) + 2*coef[70]*pow(x,3)*y + 2*coef[65]*pow(x,3)*z + coef[42]*pow(x,3) + 3*coef[75]*pow(x,2)*pow(y,2) + 4*coef[71]*pow(x,2)*y*z +
			2*coef[47]*pow(x,2)*y + 3*coef[66]*pow(x,2)*pow(z,2) + 2*coef[43]*pow(x,2)*z + coef[26]*pow(x,2) + 4*coef[79]*x*pow(y,3) + 6*coef[76]*x*pow(y,2)*z +
			3*coef[51]*x*pow(y,2) + 6*coef[72]*x*y*pow(z,2) + 4*coef[48]*x*y*z + 2*coef[30]*x*y + 4*coef[67]*x*pow(z,3) + 3*coef[44]*x*pow(z,2) +
			2*coef[27]*x*z + coef[15]*x + 5*coef[82]*pow(y,4) + 8*coef[80]*pow(y,3)*z + 4*coef[54]*pow(y,3) + 9*coef[77]*pow(y,2)*pow(z,2) +
			6*coef[52]*pow(y,2)*z + 3*coef[33]*pow(y,2) + 8*coef[73]*y*pow(z,3) + 6*coef[49]*y*pow(z,2) + 4*coef[31]*y*z + 2*coef[18]*y +
			5*coef[68]*pow(z,4) + 4*coef[45]*pow(z,3) + 3*coef[28]*pow(z,2) + 2*coef[16]*z + coef[8];
		tensor[5] = 2*coef[58]*pow(x,4) + 2*coef[65]*pow(x,3)*y + 6*coef[59]*pow(x,3)*z + 2*coef[37]*pow(x,3) + 2*coef[71]*pow(x,2)*pow(y,2) +
			6*coef[66]*pow(x,2)*y*z + 2*coef[43]*pow(x,2)*y + 12*coef[60]*pow(x,2)*pow(z,2) + 6*coef[38]*pow(x,2)*z + 2*coef[22]*pow(x,2) + 2*coef[76]*x*pow(y,3) +
			6*coef[72]*x*pow(y,2)*z + 2*coef[48]*x*pow(y,2) + 12*coef[67]*x*y*pow(z,2) + 6*coef[44]*x*y*z + 2*coef[27]*x*y + 20*coef[61]*x*pow(z,3) +
			12*coef[39]*x*pow(z,2) + 6*coef[23]*x*z + 2*coef[12]*x + 2*coef[80]*pow(y,4) + 6*coef[77]*pow(y,3)*z + 2*coef[52]*pow(y,3) +
			12*coef[73]*pow(y,2)*pow(z,2) + 6*coef[49]*pow(y,2)*z + 2*coef[31]*pow(y,2) + 20*coef[68]*y*pow(z,3) + 12*coef[45]*y*pow(z,2) + 6*coef[28]*y*z +
			2*coef[16]*y + 30*coef[62]*pow(z,4) + 20*coef[40]*pow(z,3) + 12*coef[24]*pow(z,2) + 6*coef[13]*z + 2*coef[6];
	};

	void operator()(const blocked_range3d<int>& r)	const {
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if (i+1==gridx)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, 1, i, j, k, edges, totalEdgePoints);
				}
			}
		}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if(j+1==gridy)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, 2, i, j, k, edges, totalEdgePoints);
				}
			}
		}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( int k=r.pages().begin();k!=r.pages().end();++k){
					if(k+1==gridz)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
					edgePhaseIteration(index1, index2, 3, i, j, k, edges, totalEdgePoints);
				}
			}
		}
	};
};

class CGNSParallelFace
{
	int gridx,gridy,gridz;
	Edge *edges;
	Face *faces;
	int *totalFacePoints;
	int dataType;
	float change;
	static const int subdNum = 3;
	float *gridPoints;
	float *coef;

	void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, int index1, int index2, int index3, int index4,
		int axis, int i, int j, int k, Face *lfaces, int *ltotalFacePoints) const;
	int getFaceIndex(int axis, int x, int y, int z) const;
	int getEdgeIndex(int axis, int x, int y, int z) const;
	void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int axis) const;
	void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h, int *ltotalFacePoints) const;
	void getCurveType(Face *face) const;
	void get3DWindingNum(int index1, int index3, Face *face, int axis) const;
	float smallAbsA(float A) const;
	void getGridPointPosD(int index, float v[3]) const;
	void rotate2D(float p[2], float theta) const;
	bool getEigensolverPN(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	float signedArea3D(float v1[3], float v2[3], float v3[3]) const;
	void getGridPointPos(int index, float v[3]) const;
	float getDihedral(float b1[3], float b2[3], float b3[3]) const;

	void getScalarPN(float position[3], float *scalar) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		*scalar = coef[56]*pow(x,6) + coef[63]*pow(x,5)*y + coef[57]*pow(x,5)*z + coef[35]*pow(x,5) + coef[69]*pow(x,4)*pow(y,2) + coef[64]*pow(x,4)*y*z +
			coef[41]*pow(x,4)*y + coef[58]*pow(x,4)*pow(z,2) + coef[36]*pow(x,4)*z + coef[20]*pow(x,4) + coef[74]*pow(x,3)*pow(y,3) + coef[70]*pow(x,3)*pow(y,2)*z +
			coef[46]*pow(x,3)*pow(y,2) + coef[65]*pow(x,3)*y*pow(z,2) + coef[42]*pow(x,3)*y*z + coef[25]*pow(x,3)*y + coef[59]*pow(x,3)*pow(z,3) + coef[37]*pow(x,3)*pow(z,2) +
			coef[21]*pow(x,3)*z + coef[10]*pow(x,3) + coef[78]*pow(x,2)*pow(y,4) + coef[75]*pow(x,2)*pow(y,3)*z + coef[50]*pow(x,2)*pow(y,3) + coef[71]*pow(x,2)*pow(y,2)*pow(z,2) +
			coef[47]*pow(x,2)*pow(y,2)*z + coef[29]*pow(x,2)*pow(y,2) + coef[66]*pow(x,2)*y*pow(z,3) + coef[43]*pow(x,2)*y*pow(z,2) + coef[26]*pow(x,2)*y*z +
			coef[14]*pow(x,2)*y + coef[60]*pow(x,2)*pow(z,4) + coef[38]*pow(x,2)*pow(z,3) + coef[22]*pow(x,2)*pow(z,2) + coef[11]*pow(x,2)*z + coef[4]*pow(x,2) +
			coef[81]*x*pow(y,5) + coef[79]*x*pow(y,4)*z + coef[53]*x*pow(y,4) + coef[76]*x*pow(y,3)*pow(z,2) + coef[51]*x*pow(y,3)*z + coef[32]*x*pow(y,3) +
			coef[72]*x*pow(y,2)*pow(z,3) + coef[48]*x*pow(y,2)*pow(z,2) + coef[30]*x*pow(y,2)*z + coef[17]*x*pow(y,2) + coef[67]*x*y*pow(z,4) + coef[44]*x*y*pow(z,3) +
			coef[27]*x*y*pow(z,2) + coef[15]*x*y*z + coef[7]*x*y + coef[61]*x*pow(z,5) + coef[39]*x*pow(z,4) + coef[23]*x*pow(z,3) +
			coef[12]*x*pow(z,2) + coef[5]*x*z + coef[1]*x + coef[83]*pow(y,6) + coef[82]*pow(y,5)*z + coef[55]*pow(y,5) + coef[80]*pow(y,4)*pow(z,2) +
			coef[54]*pow(y,4)*z + coef[34]*pow(y,4) + coef[77]*pow(y,3)*pow(z,3) + coef[52]*pow(y,3)*pow(z,2) + coef[33]*pow(y,3)*z + coef[19]*pow(y,3) +
			coef[73]*pow(y,2)*pow(z,4) + coef[49]*pow(y,2)*pow(z,3) + coef[31]*pow(y,2)*pow(z,2) + coef[18]*pow(y,2)*z + coef[9]*pow(y,2) + coef[68]*y*pow(z,5) +
			coef[45]*y*pow(z,4) + coef[28]*y*pow(z,3) + coef[16]*y*pow(z,2) + coef[8]*y*z + coef[3]*y + coef[62]*pow(z,6) + coef[40]*pow(z,5) +
			coef[24]*pow(z,4) + coef[13]*pow(z,3) + coef[6]*pow(z,2) + coef[2]*z + coef[0];
	};

	void getGradientPN(float position[3], float gradient[3]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		gradient[0] = 6*coef[56]*pow(x,5) + 5*coef[63]*pow(x,4)*y + 5*coef[57]*pow(x,4)*z + 5*coef[35]*pow(x,4) + 4*coef[69]*pow(x,3)*pow(y,2) +
			4*coef[64]*pow(x,3)*y*z + 4*coef[41]*pow(x,3)*y + 4*coef[58]*pow(x,3)*pow(z,2) + 4*coef[36]*pow(x,3)*z + 4*coef[20]*pow(x,3) + 3*coef[74]*pow(x,2)*pow(y,3) +
			3*coef[70]*pow(x,2)*pow(y,2)*z + 3*coef[46]*pow(x,2)*pow(y,2) + 3*coef[65]*pow(x,2)*y*pow(z,2) + 3*coef[42]*pow(x,2)*y*z + 3*coef[25]*pow(x,2)*y +
			3*coef[59]*pow(x,2)*pow(z,3) + 3*coef[37]*pow(x,2)*pow(z,2) + 3*coef[21]*pow(x,2)*z + 3*coef[10]*pow(x,2) + 2*coef[78]*x*pow(y,4) + 2*coef[75]*x*pow(y,3)*z +
			2*coef[50]*x*pow(y,3) + 2*coef[71]*x*pow(y,2)*pow(z,2) + 2*coef[47]*x*pow(y,2)*z + 2*coef[29]*x*pow(y,2) + 2*coef[66]*x*y*pow(z,3) +
			2*coef[43]*x*y*pow(z,2) + 2*coef[26]*x*y*z + 2*coef[14]*x*y + 2*coef[60]*x*pow(z,4) + 2*coef[38]*x*pow(z,3) + 2*coef[22]*x*pow(z,2) +
			2*coef[11]*x*z + 2*coef[4]*x + coef[81]*pow(y,5) + coef[79]*pow(y,4)*z + coef[53]*pow(y,4) + coef[76]*pow(y,3)*pow(z,2) + coef[51]*pow(y,3)*z +
			coef[32]*pow(y,3) + coef[72]*pow(y,2)*pow(z,3) + coef[48]*pow(y,2)*pow(z,2) + coef[30]*pow(y,2)*z + coef[17]*pow(y,2) + coef[67]*y*pow(z,4) + coef[44]*y*pow(z,3) +
			coef[27]*y*pow(z,2) + coef[15]*y*z + coef[7]*y + coef[61]*pow(z,5) + coef[39]*pow(z,4) + coef[23]*pow(z,3) + coef[12]*pow(z,2) + coef[5]*z + coef[1];
		gradient[1] = coef[63]*pow(x,5) + 2*coef[69]*pow(x,4)*y + coef[64]*pow(x,4)*z + coef[41]*pow(x,4) + 3*coef[74]*pow(x,3)*pow(y,2) + 2*coef[70]*pow(x,3)*y*z +
			2*coef[46]*pow(x,3)*y + coef[65]*pow(x,3)*pow(z,2) + coef[42]*pow(x,3)*z + coef[25]*pow(x,3) + 4*coef[78]*pow(x,2)*pow(y,3) + 3*coef[75]*pow(x,2)*pow(y,2)*z +
			3*coef[50]*pow(x,2)*pow(y,2) + 2*coef[71]*pow(x,2)*y*pow(z,2) + 2*coef[47]*pow(x,2)*y*z + 2*coef[29]*pow(x,2)*y + coef[66]*pow(x,2)*pow(z,3) + coef[43]*pow(x,2)*pow(z,2) +
			coef[26]*pow(x,2)*z + coef[14]*pow(x,2) + 5*coef[81]*x*pow(y,4) + 4*coef[79]*x*pow(y,3)*z + 4*coef[53]*x*pow(y,3) + 3*coef[76]*x*pow(y,2)*pow(z,2) +
			3*coef[51]*x*pow(y,2)*z + 3*coef[32]*x*pow(y,2) + 2*coef[72]*x*y*pow(z,3) + 2*coef[48]*x*y*pow(z,2) + 2*coef[30]*x*y*z + 2*coef[17]*x*y +
			coef[67]*x*pow(z,4) + coef[44]*x*pow(z,3) + coef[27]*x*pow(z,2) + coef[15]*x*z + coef[7]*x + 6*coef[83]*pow(y,5) + 5*coef[82]*pow(y,4)*z +
			5*coef[55]*pow(y,4) + 4*coef[80]*pow(y,3)*pow(z,2) + 4*coef[54]*pow(y,3)*z + 4*coef[34]*pow(y,3) + 3*coef[77]*pow(y,2)*pow(z,3) + 3*coef[52]*pow(y,2)*pow(z,2) +
			3*coef[33]*pow(y,2)*z + 3*coef[19]*pow(y,2) + 2*coef[73]*y*pow(z,4) + 2*coef[49]*y*pow(z,3) + 2*coef[31]*y*pow(z,2) + 2*coef[18]*y*z +
			2*coef[9]*y + coef[68]*pow(z,5) + coef[45]*pow(z,4) + coef[28]*pow(z,3) + coef[16]*pow(z,2) + coef[8]*z + coef[3];
		gradient[2] = coef[57]*pow(x,5) + coef[64]*pow(x,4)*y + 2*coef[58]*pow(x,4)*z + coef[36]*pow(x,4) + coef[70]*pow(x,3)*pow(y,2) + 2*coef[65]*pow(x,3)*y*z +
			coef[42]*pow(x,3)*y + 3*coef[59]*pow(x,3)*pow(z,2) + 2*coef[37]*pow(x,3)*z + coef[21]*pow(x,3) + coef[75]*pow(x,2)*pow(y,3) + 2*coef[71]*pow(x,2)*pow(y,2)*z +
			coef[47]*pow(x,2)*pow(y,2) + 3*coef[66]*pow(x,2)*y*pow(z,2) + 2*coef[43]*pow(x,2)*y*z + coef[26]*pow(x,2)*y + 4*coef[60]*pow(x,2)*pow(z,3) + 3*coef[38]*pow(x,2)*pow(z,2) +
			2*coef[22]*pow(x,2)*z + coef[11]*pow(x,2) + coef[79]*x*pow(y,4) + 2*coef[76]*x*pow(y,3)*z + coef[51]*x*pow(y,3) + 3*coef[72]*x*pow(y,2)*pow(z,2) +
			2*coef[48]*x*pow(y,2)*z + coef[30]*x*pow(y,2) + 4*coef[67]*x*y*pow(z,3) + 3*coef[44]*x*y*pow(z,2) + 2*coef[27]*x*y*z + coef[15]*x*y +
			5*coef[61]*x*pow(z,4) + 4*coef[39]*x*pow(z,3) + 3*coef[23]*x*pow(z,2) + 2*coef[12]*x*z + coef[5]*x + coef[82]*pow(y,5) + 2*coef[80]*pow(y,4)*z +
			coef[54]*pow(y,4) + 3*coef[77]*pow(y,3)*pow(z,2) + 2*coef[52]*pow(y,3)*z + coef[33]*pow(y,3) + 4*coef[73]*pow(y,2)*pow(z,3) + 3*coef[49]*pow(y,2)*pow(z,2) +
			2*coef[31]*pow(y,2)*z + coef[18]*pow(y,2) + 5*coef[68]*y*pow(z,4) + 4*coef[45]*y*pow(z,3) + 3*coef[28]*y*pow(z,2) + 2*coef[16]*y*z + coef[8]*y +
			6*coef[62]*pow(z,5) + 5*coef[40]*pow(z,4) + 4*coef[24]*pow(z,3) + 3*coef[13]*pow(z,2) + 2*coef[6]*z + coef[2];
	};

	void getTensorPN(float position[3], float tensor[6]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		tensor[0] = 30*coef[56]*pow(x,4) + 20*coef[63]*pow(x,3)*y + 20*coef[57]*pow(x,3)*z + 20*coef[35]*pow(x,3) + 12*coef[69]*pow(x,2)*pow(y,2) +
			12*coef[64]*pow(x,2)*y*z + 12*coef[41]*pow(x,2)*y + 12*coef[58]*pow(x,2)*pow(z,2) + 12*coef[36]*pow(x,2)*z + 12*coef[20]*pow(x,2) + 6*coef[74]*x*pow(y,3) +
			6*coef[70]*x*pow(y,2)*z + 6*coef[46]*x*pow(y,2) + 6*coef[65]*x*y*pow(z,2) + 6*coef[42]*x*y*z + 6*coef[25]*x*y + 6*coef[59]*x*pow(z,3) +
			6*coef[37]*x*pow(z,2) + 6*coef[21]*x*z + 6*coef[10]*x + 2*coef[78]*pow(y,4) + 2*coef[75]*pow(y,3)*z + 2*coef[50]*pow(y,3) +
			2*coef[71]*pow(y,2)*pow(z,2) + 2*coef[47]*pow(y,2)*z + 2*coef[29]*pow(y,2) + 2*coef[66]*y*pow(z,3) + 2*coef[43]*y*pow(z,2) + 2*coef[26]*y*z +
			2*coef[14]*y + 2*coef[60]*pow(z,4) + 2*coef[38]*pow(z,3) + 2*coef[22]*pow(z,2) + 2*coef[11]*z + 2*coef[4];
		tensor[1] = 5*coef[63]*pow(x,4) + 8*coef[69]*pow(x,3)*y + 4*coef[64]*pow(x,3)*z + 4*coef[41]*pow(x,3) + 9*coef[74]*pow(x,2)*pow(y,2) +
			6*coef[70]*pow(x,2)*y*z + 6*coef[46]*pow(x,2)*y + 3*coef[65]*pow(x,2)*pow(z,2) + 3*coef[42]*pow(x,2)*z + 3*coef[25]*pow(x,2) + 8*coef[78]*x*pow(y,3) +
			6*coef[75]*x*pow(y,2)*z + 6*coef[50]*x*pow(y,2) + 4*coef[71]*x*y*pow(z,2) + 4*coef[47]*x*y*z + 4*coef[29]*x*y + 2*coef[66]*x*pow(z,3) +
			2*coef[43]*x*pow(z,2) + 2*coef[26]*x*z + 2*coef[14]*x + 5*coef[81]*pow(y,4) + 4*coef[79]*pow(y,3)*z + 4*coef[53]*pow(y,3) +
			3*coef[76]*pow(y,2)*pow(z,2) + 3*coef[51]*pow(y,2)*z + 3*coef[32]*pow(y,2) + 2*coef[72]*y*pow(z,3) + 2*coef[48]*y*pow(z,2) + 2*coef[30]*y*z +
			2*coef[17]*y + coef[67]*pow(z,4) + coef[44]*pow(z,3) + coef[27]*pow(z,2) + coef[15]*z + coef[7];
		tensor[2] = 5*coef[57]*pow(x,4) + 4*coef[64]*pow(x,3)*y + 8*coef[58]*pow(x,3)*z + 4*coef[36]*pow(x,3) + 3*coef[70]*pow(x,2)*pow(y,2) +
			6*coef[65]*pow(x,2)*y*z + 3*coef[42]*pow(x,2)*y + 9*coef[59]*pow(x,2)*pow(z,2) + 6*coef[37]*pow(x,2)*z + 3*coef[21]*pow(x,2) + 2*coef[75]*x*pow(y,3) +
			4*coef[71]*x*pow(y,2)*z + 2*coef[47]*x*pow(y,2) + 6*coef[66]*x*y*pow(z,2) + 4*coef[43]*x*y*z + 2*coef[26]*x*y + 8*coef[60]*x*pow(z,3) +
			6*coef[38]*x*pow(z,2) + 4*coef[22]*x*z + 2*coef[11]*x + coef[79]*pow(y,4) + 2*coef[76]*pow(y,3)*z + coef[51]*pow(y,3) +
			3*coef[72]*pow(y,2)*pow(z,2) + 2*coef[48]*pow(y,2)*z + coef[30]*pow(y,2) + 4*coef[67]*y*pow(z,3) + 3*coef[44]*y*pow(z,2) + 2*coef[27]*y*z +
			coef[15]*y + 5*coef[61]*pow(z,4) + 4*coef[39]*pow(z,3) + 3*coef[23]*pow(z,2) + 2*coef[12]*z + coef[5];
		tensor[3] = 2*coef[69]*pow(x,4) + 6*coef[74]*pow(x,3)*y + 2*coef[70]*pow(x,3)*z + 2*coef[46]*pow(x,3) + 12*coef[78]*pow(x,2)*pow(y,2) +
			6*coef[75]*pow(x,2)*y*z + 6*coef[50]*pow(x,2)*y + 2*coef[71]*pow(x,2)*pow(z,2) + 2*coef[47]*pow(x,2)*z + 2*coef[29]*pow(x,2) + 20*coef[81]*x*pow(y,3) +
			12*coef[79]*x*pow(y,2)*z + 12*coef[53]*x*pow(y,2) + 6*coef[76]*x*y*pow(z,2) + 6*coef[51]*x*y*z + 6*coef[32]*x*y + 2*coef[72]*x*pow(z,3) +
			2*coef[48]*x*pow(z,2) + 2*coef[30]*x*z + 2*coef[17]*x + 30*coef[83]*pow(y,4) + 20*coef[82]*pow(y,3)*z + 20*coef[55]*pow(y,3) +
			12*coef[80]*pow(y,2)*pow(z,2) + 12*coef[54]*pow(y,2)*z + 12*coef[34]*pow(y,2) + 6*coef[77]*y*pow(z,3) + 6*coef[52]*y*pow(z,2) + 6*coef[33]*y*z +
			6*coef[19]*y + 2*coef[73]*pow(z,4) + 2*coef[49]*pow(z,3) + 2*coef[31]*pow(z,2) + 2*coef[18]*z + 2*coef[9];
		tensor[4] = coef[64]*pow(x,4) + 2*coef[70]*pow(x,3)*y + 2*coef[65]*pow(x,3)*z + coef[42]*pow(x,3) + 3*coef[75]*pow(x,2)*pow(y,2) + 4*coef[71]*pow(x,2)*y*z +
			2*coef[47]*pow(x,2)*y + 3*coef[66]*pow(x,2)*pow(z,2) + 2*coef[43]*pow(x,2)*z + coef[26]*pow(x,2) + 4*coef[79]*x*pow(y,3) + 6*coef[76]*x*pow(y,2)*z +
			3*coef[51]*x*pow(y,2) + 6*coef[72]*x*y*pow(z,2) + 4*coef[48]*x*y*z + 2*coef[30]*x*y + 4*coef[67]*x*pow(z,3) + 3*coef[44]*x*pow(z,2) +
			2*coef[27]*x*z + coef[15]*x + 5*coef[82]*pow(y,4) + 8*coef[80]*pow(y,3)*z + 4*coef[54]*pow(y,3) + 9*coef[77]*pow(y,2)*pow(z,2) +
			6*coef[52]*pow(y,2)*z + 3*coef[33]*pow(y,2) + 8*coef[73]*y*pow(z,3) + 6*coef[49]*y*pow(z,2) + 4*coef[31]*y*z + 2*coef[18]*y +
			5*coef[68]*pow(z,4) + 4*coef[45]*pow(z,3) + 3*coef[28]*pow(z,2) + 2*coef[16]*z + coef[8];
		tensor[5] = 2*coef[58]*pow(x,4) + 2*coef[65]*pow(x,3)*y + 6*coef[59]*pow(x,3)*z + 2*coef[37]*pow(x,3) + 2*coef[71]*pow(x,2)*pow(y,2) +
			6*coef[66]*pow(x,2)*y*z + 2*coef[43]*pow(x,2)*y + 12*coef[60]*pow(x,2)*pow(z,2) + 6*coef[38]*pow(x,2)*z + 2*coef[22]*pow(x,2) + 2*coef[76]*x*pow(y,3) +
			6*coef[72]*x*pow(y,2)*z + 2*coef[48]*x*pow(y,2) + 12*coef[67]*x*y*pow(z,2) + 6*coef[44]*x*y*z + 2*coef[27]*x*y + 20*coef[61]*x*pow(z,3) +
			12*coef[39]*x*pow(z,2) + 6*coef[23]*x*z + 2*coef[12]*x + 2*coef[80]*pow(y,4) + 6*coef[77]*pow(y,3)*z + 2*coef[52]*pow(y,3) +
			12*coef[73]*pow(y,2)*pow(z,2) + 6*coef[49]*pow(y,2)*z + 2*coef[31]*pow(y,2) + 20*coef[68]*y*pow(z,3) + 12*coef[45]*y*pow(z,2) + 6*coef[28]*y*z +
			2*coef[16]*y + 30*coef[62]*pow(z,4) + 20*coef[40]*pow(z,3) + 12*coef[24]*pow(z,2) + 6*coef[13]*z + 2*coef[6];
	};

public: 
	bool *storeGrid_2DGradMag;
	float *grid_2DGradMag;

	CGNSParallelFace(int lgridx,int lgridy,int lgridz,int *ltotalFacePoints,Edge *ledges,Face *lfaces,int ldataType,
		float *lgridPoints,float lchange,bool *lstoreGrid_2DGradMag,float *lgrid_2DGradMag, float *lcoef)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalFacePoints=ltotalFacePoints;
		edges = ledges;
		faces = lfaces;
		dataType = ldataType;
		gridPoints = lgridPoints;
		change=lchange;
		storeGrid_2DGradMag = lstoreGrid_2DGradMag;
		grid_2DGradMag = lgrid_2DGradMag;
		coef = lcoef;
	};

	void operator()(const blocked_range3d<int>& r)	const {
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( j == gridy - 1)
						continue;
					if ( k == gridz - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i, j + 1, k + 1, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
					Edge *edge1 = &(edges[getEdgeIndex(2, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(3, i, j + 1, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(2, i, j, k + 1)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
					if (!edge4->valid) continue;
					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 1, i, j, k, faces, totalFacePoints);
				}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( i == gridx - 1)
						continue;
					if ( k == gridz - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i + 1, j, k + 1, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
					Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(3, i + 1, j, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(1, i, j, k + 1)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
					if (!edge4->valid) continue;
					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 2, i, j, k, faces, totalFacePoints);
				}
		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if ( i == gridx - 1)
						continue;
					if ( j == gridy - 1)
						continue;
					int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
					int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
					int index3 = getIndex(1, 0, i + 1, j + 1, k, gridx, gridy, gridz);
					int index4 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
					Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
					if (!edge1->valid) continue;
					Edge *edge2 = &(edges[getEdgeIndex(2, i + 1, j, k)]);
					if (!edge2->valid) continue;
					Edge *edge3 = &(edges[getEdgeIndex(1, i, j + 1, k)]);
					if (!edge3->valid) continue;
					Edge *edge4 = &(edges[getEdgeIndex(2, i, j, k)]);
					if (!edge4->valid) continue;
					facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 3, i, j, k, faces, totalFacePoints);
				}
	};
};

class CGNSParallelCell
{
	int gridx,gridy,gridz;
	Edge *edges;
	Face *faces;
	Cell *cells;
	vector<Vertex> *vertices;
	vector<Point> *points;
	concurrent_vector<Point> *conpoints;
	concurrent_vector<Vertex> *convertices;
	float *maxEigenvalue, *minEigenvalue;
	float halfDataSize, *rightBackTop;
	float thickness;
	int *totalCellPoints;
	float gridSize;
	int dataType;
	float *gridPoints;
	float change;
	float *coef;

	void cellPhaseIteration(int i, int j, int k, Cell *cells, concurrent_vector<Point> *lpoints,
		concurrent_vector<Vertex> *lvertices, int *ltotalCellPoints, float* lmaxEigenvalue, float* lminEigenvalue) const;
	bool getCentroid(int i, int j, int k, Cell *cell,int *ltotalCellPoints) const;
	void centroidProjection(int i, int j, int k, Cell *cell) const;
	void buildPointIteration(Cell *cell,concurrent_vector<Point> *points,concurrent_vector<Vertex> *vertices,
		float* lmaxEigenvalue, float* lminEigenvalue) const;
	void getProjDirection(int type, float position[3], float q[3]) const;
	int getEdgeIndex(int axis, int x, int y, int z) const;
	int getFaceIndex(int axis, int x, int y, int z) const;
	void getGridPointPosD(int index, float v[3]) const;
	void getShowPos(float position[3], float showPos[3]) const;
	bool getEigensolverPN(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
	float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3],float* maxEigenvalue, float* minEigenvalue) const;
	void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
	void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;

	void getScalarPN(float position[3], float *scalar) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		*scalar = coef[56]*pow(x,6) + coef[63]*pow(x,5)*y + coef[57]*pow(x,5)*z + coef[35]*pow(x,5) + coef[69]*pow(x,4)*pow(y,2) + coef[64]*pow(x,4)*y*z +
			coef[41]*pow(x,4)*y + coef[58]*pow(x,4)*pow(z,2) + coef[36]*pow(x,4)*z + coef[20]*pow(x,4) + coef[74]*pow(x,3)*pow(y,3) + coef[70]*pow(x,3)*pow(y,2)*z +
			coef[46]*pow(x,3)*pow(y,2) + coef[65]*pow(x,3)*y*pow(z,2) + coef[42]*pow(x,3)*y*z + coef[25]*pow(x,3)*y + coef[59]*pow(x,3)*pow(z,3) + coef[37]*pow(x,3)*pow(z,2) +
			coef[21]*pow(x,3)*z + coef[10]*pow(x,3) + coef[78]*pow(x,2)*pow(y,4) + coef[75]*pow(x,2)*pow(y,3)*z + coef[50]*pow(x,2)*pow(y,3) + coef[71]*pow(x,2)*pow(y,2)*pow(z,2) +
			coef[47]*pow(x,2)*pow(y,2)*z + coef[29]*pow(x,2)*pow(y,2) + coef[66]*pow(x,2)*y*pow(z,3) + coef[43]*pow(x,2)*y*pow(z,2) + coef[26]*pow(x,2)*y*z +
			coef[14]*pow(x,2)*y + coef[60]*pow(x,2)*pow(z,4) + coef[38]*pow(x,2)*pow(z,3) + coef[22]*pow(x,2)*pow(z,2) + coef[11]*pow(x,2)*z + coef[4]*pow(x,2) +
			coef[81]*x*pow(y,5) + coef[79]*x*pow(y,4)*z + coef[53]*x*pow(y,4) + coef[76]*x*pow(y,3)*pow(z,2) + coef[51]*x*pow(y,3)*z + coef[32]*x*pow(y,3) +
			coef[72]*x*pow(y,2)*pow(z,3) + coef[48]*x*pow(y,2)*pow(z,2) + coef[30]*x*pow(y,2)*z + coef[17]*x*pow(y,2) + coef[67]*x*y*pow(z,4) + coef[44]*x*y*pow(z,3) +
			coef[27]*x*y*pow(z,2) + coef[15]*x*y*z + coef[7]*x*y + coef[61]*x*pow(z,5) + coef[39]*x*pow(z,4) + coef[23]*x*pow(z,3) +
			coef[12]*x*pow(z,2) + coef[5]*x*z + coef[1]*x + coef[83]*pow(y,6) + coef[82]*pow(y,5)*z + coef[55]*pow(y,5) + coef[80]*pow(y,4)*pow(z,2) +
			coef[54]*pow(y,4)*z + coef[34]*pow(y,4) + coef[77]*pow(y,3)*pow(z,3) + coef[52]*pow(y,3)*pow(z,2) + coef[33]*pow(y,3)*z + coef[19]*pow(y,3) +
			coef[73]*pow(y,2)*pow(z,4) + coef[49]*pow(y,2)*pow(z,3) + coef[31]*pow(y,2)*pow(z,2) + coef[18]*pow(y,2)*z + coef[9]*pow(y,2) + coef[68]*y*pow(z,5) +
			coef[45]*y*pow(z,4) + coef[28]*y*pow(z,3) + coef[16]*y*pow(z,2) + coef[8]*y*z + coef[3]*y + coef[62]*pow(z,6) + coef[40]*pow(z,5) +
			coef[24]*pow(z,4) + coef[13]*pow(z,3) + coef[6]*pow(z,2) + coef[2]*z + coef[0];
	};

	void getGradientPN(float position[3], float gradient[3]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		gradient[0] = 6*coef[56]*pow(x,5) + 5*coef[63]*pow(x,4)*y + 5*coef[57]*pow(x,4)*z + 5*coef[35]*pow(x,4) + 4*coef[69]*pow(x,3)*pow(y,2) +
			4*coef[64]*pow(x,3)*y*z + 4*coef[41]*pow(x,3)*y + 4*coef[58]*pow(x,3)*pow(z,2) + 4*coef[36]*pow(x,3)*z + 4*coef[20]*pow(x,3) + 3*coef[74]*pow(x,2)*pow(y,3) +
			3*coef[70]*pow(x,2)*pow(y,2)*z + 3*coef[46]*pow(x,2)*pow(y,2) + 3*coef[65]*pow(x,2)*y*pow(z,2) + 3*coef[42]*pow(x,2)*y*z + 3*coef[25]*pow(x,2)*y +
			3*coef[59]*pow(x,2)*pow(z,3) + 3*coef[37]*pow(x,2)*pow(z,2) + 3*coef[21]*pow(x,2)*z + 3*coef[10]*pow(x,2) + 2*coef[78]*x*pow(y,4) + 2*coef[75]*x*pow(y,3)*z +
			2*coef[50]*x*pow(y,3) + 2*coef[71]*x*pow(y,2)*pow(z,2) + 2*coef[47]*x*pow(y,2)*z + 2*coef[29]*x*pow(y,2) + 2*coef[66]*x*y*pow(z,3) +
			2*coef[43]*x*y*pow(z,2) + 2*coef[26]*x*y*z + 2*coef[14]*x*y + 2*coef[60]*x*pow(z,4) + 2*coef[38]*x*pow(z,3) + 2*coef[22]*x*pow(z,2) +
			2*coef[11]*x*z + 2*coef[4]*x + coef[81]*pow(y,5) + coef[79]*pow(y,4)*z + coef[53]*pow(y,4) + coef[76]*pow(y,3)*pow(z,2) + coef[51]*pow(y,3)*z +
			coef[32]*pow(y,3) + coef[72]*pow(y,2)*pow(z,3) + coef[48]*pow(y,2)*pow(z,2) + coef[30]*pow(y,2)*z + coef[17]*pow(y,2) + coef[67]*y*pow(z,4) + coef[44]*y*pow(z,3) +
			coef[27]*y*pow(z,2) + coef[15]*y*z + coef[7]*y + coef[61]*pow(z,5) + coef[39]*pow(z,4) + coef[23]*pow(z,3) + coef[12]*pow(z,2) + coef[5]*z + coef[1];
		gradient[1] = coef[63]*pow(x,5) + 2*coef[69]*pow(x,4)*y + coef[64]*pow(x,4)*z + coef[41]*pow(x,4) + 3*coef[74]*pow(x,3)*pow(y,2) + 2*coef[70]*pow(x,3)*y*z +
			2*coef[46]*pow(x,3)*y + coef[65]*pow(x,3)*pow(z,2) + coef[42]*pow(x,3)*z + coef[25]*pow(x,3) + 4*coef[78]*pow(x,2)*pow(y,3) + 3*coef[75]*pow(x,2)*pow(y,2)*z +
			3*coef[50]*pow(x,2)*pow(y,2) + 2*coef[71]*pow(x,2)*y*pow(z,2) + 2*coef[47]*pow(x,2)*y*z + 2*coef[29]*pow(x,2)*y + coef[66]*pow(x,2)*pow(z,3) + coef[43]*pow(x,2)*pow(z,2) +
			coef[26]*pow(x,2)*z + coef[14]*pow(x,2) + 5*coef[81]*x*pow(y,4) + 4*coef[79]*x*pow(y,3)*z + 4*coef[53]*x*pow(y,3) + 3*coef[76]*x*pow(y,2)*pow(z,2) +
			3*coef[51]*x*pow(y,2)*z + 3*coef[32]*x*pow(y,2) + 2*coef[72]*x*y*pow(z,3) + 2*coef[48]*x*y*pow(z,2) + 2*coef[30]*x*y*z + 2*coef[17]*x*y +
			coef[67]*x*pow(z,4) + coef[44]*x*pow(z,3) + coef[27]*x*pow(z,2) + coef[15]*x*z + coef[7]*x + 6*coef[83]*pow(y,5) + 5*coef[82]*pow(y,4)*z +
			5*coef[55]*pow(y,4) + 4*coef[80]*pow(y,3)*pow(z,2) + 4*coef[54]*pow(y,3)*z + 4*coef[34]*pow(y,3) + 3*coef[77]*pow(y,2)*pow(z,3) + 3*coef[52]*pow(y,2)*pow(z,2) +
			3*coef[33]*pow(y,2)*z + 3*coef[19]*pow(y,2) + 2*coef[73]*y*pow(z,4) + 2*coef[49]*y*pow(z,3) + 2*coef[31]*y*pow(z,2) + 2*coef[18]*y*z +
			2*coef[9]*y + coef[68]*pow(z,5) + coef[45]*pow(z,4) + coef[28]*pow(z,3) + coef[16]*pow(z,2) + coef[8]*z + coef[3];
		gradient[2] = coef[57]*pow(x,5) + coef[64]*pow(x,4)*y + 2*coef[58]*pow(x,4)*z + coef[36]*pow(x,4) + coef[70]*pow(x,3)*pow(y,2) + 2*coef[65]*pow(x,3)*y*z +
			coef[42]*pow(x,3)*y + 3*coef[59]*pow(x,3)*pow(z,2) + 2*coef[37]*pow(x,3)*z + coef[21]*pow(x,3) + coef[75]*pow(x,2)*pow(y,3) + 2*coef[71]*pow(x,2)*pow(y,2)*z +
			coef[47]*pow(x,2)*pow(y,2) + 3*coef[66]*pow(x,2)*y*pow(z,2) + 2*coef[43]*pow(x,2)*y*z + coef[26]*pow(x,2)*y + 4*coef[60]*pow(x,2)*pow(z,3) + 3*coef[38]*pow(x,2)*pow(z,2) +
			2*coef[22]*pow(x,2)*z + coef[11]*pow(x,2) + coef[79]*x*pow(y,4) + 2*coef[76]*x*pow(y,3)*z + coef[51]*x*pow(y,3) + 3*coef[72]*x*pow(y,2)*pow(z,2) +
			2*coef[48]*x*pow(y,2)*z + coef[30]*x*pow(y,2) + 4*coef[67]*x*y*pow(z,3) + 3*coef[44]*x*y*pow(z,2) + 2*coef[27]*x*y*z + coef[15]*x*y +
			5*coef[61]*x*pow(z,4) + 4*coef[39]*x*pow(z,3) + 3*coef[23]*x*pow(z,2) + 2*coef[12]*x*z + coef[5]*x + coef[82]*pow(y,5) + 2*coef[80]*pow(y,4)*z +
			coef[54]*pow(y,4) + 3*coef[77]*pow(y,3)*pow(z,2) + 2*coef[52]*pow(y,3)*z + coef[33]*pow(y,3) + 4*coef[73]*pow(y,2)*pow(z,3) + 3*coef[49]*pow(y,2)*pow(z,2) +
			2*coef[31]*pow(y,2)*z + coef[18]*pow(y,2) + 5*coef[68]*y*pow(z,4) + 4*coef[45]*y*pow(z,3) + 3*coef[28]*y*pow(z,2) + 2*coef[16]*y*z + coef[8]*y +
			6*coef[62]*pow(z,5) + 5*coef[40]*pow(z,4) + 4*coef[24]*pow(z,3) + 3*coef[13]*pow(z,2) + 2*coef[6]*z + coef[2];
	};

	void getTensorPN(float position[3], float tensor[6]) const
	{
		float x = position[0];
		float y = position[1];
		float z = position[2];
		tensor[0] = 30*coef[56]*pow(x,4) + 20*coef[63]*pow(x,3)*y + 20*coef[57]*pow(x,3)*z + 20*coef[35]*pow(x,3) + 12*coef[69]*pow(x,2)*pow(y,2) +
			12*coef[64]*pow(x,2)*y*z + 12*coef[41]*pow(x,2)*y + 12*coef[58]*pow(x,2)*pow(z,2) + 12*coef[36]*pow(x,2)*z + 12*coef[20]*pow(x,2) + 6*coef[74]*x*pow(y,3) +
			6*coef[70]*x*pow(y,2)*z + 6*coef[46]*x*pow(y,2) + 6*coef[65]*x*y*pow(z,2) + 6*coef[42]*x*y*z + 6*coef[25]*x*y + 6*coef[59]*x*pow(z,3) +
			6*coef[37]*x*pow(z,2) + 6*coef[21]*x*z + 6*coef[10]*x + 2*coef[78]*pow(y,4) + 2*coef[75]*pow(y,3)*z + 2*coef[50]*pow(y,3) +
			2*coef[71]*pow(y,2)*pow(z,2) + 2*coef[47]*pow(y,2)*z + 2*coef[29]*pow(y,2) + 2*coef[66]*y*pow(z,3) + 2*coef[43]*y*pow(z,2) + 2*coef[26]*y*z +
			2*coef[14]*y + 2*coef[60]*pow(z,4) + 2*coef[38]*pow(z,3) + 2*coef[22]*pow(z,2) + 2*coef[11]*z + 2*coef[4];
		tensor[1] = 5*coef[63]*pow(x,4) + 8*coef[69]*pow(x,3)*y + 4*coef[64]*pow(x,3)*z + 4*coef[41]*pow(x,3) + 9*coef[74]*pow(x,2)*pow(y,2) +
			6*coef[70]*pow(x,2)*y*z + 6*coef[46]*pow(x,2)*y + 3*coef[65]*pow(x,2)*pow(z,2) + 3*coef[42]*pow(x,2)*z + 3*coef[25]*pow(x,2) + 8*coef[78]*x*pow(y,3) +
			6*coef[75]*x*pow(y,2)*z + 6*coef[50]*x*pow(y,2) + 4*coef[71]*x*y*pow(z,2) + 4*coef[47]*x*y*z + 4*coef[29]*x*y + 2*coef[66]*x*pow(z,3) +
			2*coef[43]*x*pow(z,2) + 2*coef[26]*x*z + 2*coef[14]*x + 5*coef[81]*pow(y,4) + 4*coef[79]*pow(y,3)*z + 4*coef[53]*pow(y,3) +
			3*coef[76]*pow(y,2)*pow(z,2) + 3*coef[51]*pow(y,2)*z + 3*coef[32]*pow(y,2) + 2*coef[72]*y*pow(z,3) + 2*coef[48]*y*pow(z,2) + 2*coef[30]*y*z +
			2*coef[17]*y + coef[67]*pow(z,4) + coef[44]*pow(z,3) + coef[27]*pow(z,2) + coef[15]*z + coef[7];
		tensor[2] = 5*coef[57]*pow(x,4) + 4*coef[64]*pow(x,3)*y + 8*coef[58]*pow(x,3)*z + 4*coef[36]*pow(x,3) + 3*coef[70]*pow(x,2)*pow(y,2) +
			6*coef[65]*pow(x,2)*y*z + 3*coef[42]*pow(x,2)*y + 9*coef[59]*pow(x,2)*pow(z,2) + 6*coef[37]*pow(x,2)*z + 3*coef[21]*pow(x,2) + 2*coef[75]*x*pow(y,3) +
			4*coef[71]*x*pow(y,2)*z + 2*coef[47]*x*pow(y,2) + 6*coef[66]*x*y*pow(z,2) + 4*coef[43]*x*y*z + 2*coef[26]*x*y + 8*coef[60]*x*pow(z,3) +
			6*coef[38]*x*pow(z,2) + 4*coef[22]*x*z + 2*coef[11]*x + coef[79]*pow(y,4) + 2*coef[76]*pow(y,3)*z + coef[51]*pow(y,3) +
			3*coef[72]*pow(y,2)*pow(z,2) + 2*coef[48]*pow(y,2)*z + coef[30]*pow(y,2) + 4*coef[67]*y*pow(z,3) + 3*coef[44]*y*pow(z,2) + 2*coef[27]*y*z +
			coef[15]*y + 5*coef[61]*pow(z,4) + 4*coef[39]*pow(z,3) + 3*coef[23]*pow(z,2) + 2*coef[12]*z + coef[5];
		tensor[3] = 2*coef[69]*pow(x,4) + 6*coef[74]*pow(x,3)*y + 2*coef[70]*pow(x,3)*z + 2*coef[46]*pow(x,3) + 12*coef[78]*pow(x,2)*pow(y,2) +
			6*coef[75]*pow(x,2)*y*z + 6*coef[50]*pow(x,2)*y + 2*coef[71]*pow(x,2)*pow(z,2) + 2*coef[47]*pow(x,2)*z + 2*coef[29]*pow(x,2) + 20*coef[81]*x*pow(y,3) +
			12*coef[79]*x*pow(y,2)*z + 12*coef[53]*x*pow(y,2) + 6*coef[76]*x*y*pow(z,2) + 6*coef[51]*x*y*z + 6*coef[32]*x*y + 2*coef[72]*x*pow(z,3) +
			2*coef[48]*x*pow(z,2) + 2*coef[30]*x*z + 2*coef[17]*x + 30*coef[83]*pow(y,4) + 20*coef[82]*pow(y,3)*z + 20*coef[55]*pow(y,3) +
			12*coef[80]*pow(y,2)*pow(z,2) + 12*coef[54]*pow(y,2)*z + 12*coef[34]*pow(y,2) + 6*coef[77]*y*pow(z,3) + 6*coef[52]*y*pow(z,2) + 6*coef[33]*y*z +
			6*coef[19]*y + 2*coef[73]*pow(z,4) + 2*coef[49]*pow(z,3) + 2*coef[31]*pow(z,2) + 2*coef[18]*z + 2*coef[9];
		tensor[4] = coef[64]*pow(x,4) + 2*coef[70]*pow(x,3)*y + 2*coef[65]*pow(x,3)*z + coef[42]*pow(x,3) + 3*coef[75]*pow(x,2)*pow(y,2) + 4*coef[71]*pow(x,2)*y*z +
			2*coef[47]*pow(x,2)*y + 3*coef[66]*pow(x,2)*pow(z,2) + 2*coef[43]*pow(x,2)*z + coef[26]*pow(x,2) + 4*coef[79]*x*pow(y,3) + 6*coef[76]*x*pow(y,2)*z +
			3*coef[51]*x*pow(y,2) + 6*coef[72]*x*y*pow(z,2) + 4*coef[48]*x*y*z + 2*coef[30]*x*y + 4*coef[67]*x*pow(z,3) + 3*coef[44]*x*pow(z,2) +
			2*coef[27]*x*z + coef[15]*x + 5*coef[82]*pow(y,4) + 8*coef[80]*pow(y,3)*z + 4*coef[54]*pow(y,3) + 9*coef[77]*pow(y,2)*pow(z,2) +
			6*coef[52]*pow(y,2)*z + 3*coef[33]*pow(y,2) + 8*coef[73]*y*pow(z,3) + 6*coef[49]*y*pow(z,2) + 4*coef[31]*y*z + 2*coef[18]*y +
			5*coef[68]*pow(z,4) + 4*coef[45]*pow(z,3) + 3*coef[28]*pow(z,2) + 2*coef[16]*z + coef[8];
		tensor[5] = 2*coef[58]*pow(x,4) + 2*coef[65]*pow(x,3)*y + 6*coef[59]*pow(x,3)*z + 2*coef[37]*pow(x,3) + 2*coef[71]*pow(x,2)*pow(y,2) +
			6*coef[66]*pow(x,2)*y*z + 2*coef[43]*pow(x,2)*y + 12*coef[60]*pow(x,2)*pow(z,2) + 6*coef[38]*pow(x,2)*z + 2*coef[22]*pow(x,2) + 2*coef[76]*x*pow(y,3) +
			6*coef[72]*x*pow(y,2)*z + 2*coef[48]*x*pow(y,2) + 12*coef[67]*x*y*pow(z,2) + 6*coef[44]*x*y*z + 2*coef[27]*x*y + 20*coef[61]*x*pow(z,3) +
			12*coef[39]*x*pow(z,2) + 6*coef[23]*x*z + 2*coef[12]*x + 2*coef[80]*pow(y,4) + 6*coef[77]*pow(y,3)*z + 2*coef[52]*pow(y,3) +
			12*coef[73]*pow(y,2)*pow(z,2) + 6*coef[49]*pow(y,2)*z + 2*coef[31]*pow(y,2) + 20*coef[68]*y*pow(z,3) + 12*coef[45]*y*pow(z,2) + 6*coef[28]*y*z +
			2*coef[16]*y + 30*coef[62]*pow(z,4) + 20*coef[40]*pow(z,3) + 12*coef[24]*pow(z,2) + 6*coef[13]*z + 2*coef[6];
	};

public:
	CGNSParallelCell(int lgridx,int lgridy,int lgridz,int *ltotalCellPoints,Edge *ledges,Face *lfaces,Cell *lcells,
		concurrent_vector<Point> *lpoints,int ldataType,float *lgridPoints,float lchange,float lgridSize,
		float* lmaxEigenvalue, float* lminEigenvalue,concurrent_vector<Vertex> *lvertices,float lhalfDataSize,
		float *lrightBackTop,float lthickness, float *lcoef)
	{
		gridx=lgridx;
		gridy=lgridy;
		gridz=lgridz;
		totalCellPoints=ltotalCellPoints;
		edges = ledges;
		faces = lfaces;
		cells = lcells;
		conpoints = lpoints;
		convertices = lvertices;
		dataType = ldataType;
		gridPoints = lgridPoints;
		change=lchange;
		gridSize = lgridSize;
		maxEigenvalue=lmaxEigenvalue;
		minEigenvalue=lminEigenvalue;
		halfDataSize = lhalfDataSize;
		rightBackTop = lrightBackTop;
		thickness = lthickness;
		coef = lcoef;
	};

	void operator()(const blocked_range3d<int>& r)	const {
		int *ltotalCellPoints= totalCellPoints;
		Cell *lcells = cells;
		concurrent_vector<Point> *lpoints = conpoints; 
		float *lmaxEigenvalue = maxEigenvalue;
		float *lminEigenvalue = minEigenvalue;
		concurrent_vector<Vertex> *lvertices = convertices;

		for( int i=r.cols().begin(); i!=r.cols().end(); ++i )
			for( int j=r.rows().begin(); j!=r.rows().end(); ++j )
				for ( int k=r.pages().begin();k!=r.pages().end();++k )
				{
					if(i == gridx - 1 || j == gridy - 1 || k == gridz - 1)
						continue;
					cellPhaseIteration(i, j, k,lcells,lpoints,lvertices,ltotalCellPoints,lmaxEigenvalue,lminEigenvalue);
				}
	};
};
