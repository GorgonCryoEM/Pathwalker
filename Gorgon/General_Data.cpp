#include "General_Data.h"

General_Data::General_Data(View *_view3D)
{
	view3D = _view3D;
	allCubic = view3D->allCubic;
	reverse = view3D->reverse;
	kernelSize = view3D->kernelSize;
	height = view3D->height;
	width = view3D->width;
	slices = view3D->slices;
	thickness = view3D->thickness;
	resolution = view3D->resolution;
	dataType = view3D->dataType;

	storeGridShowPos = storeGridScalar = storeGridVector = storeGridGradient = 
		storeEdgePoint = storeFacePoint = storeCellPoint = storeDispCell = false;

	maxEigenvalue = new float;
	minEigenvalue = new float;
	maxIntensity = *maxEigenvalue = -1e12;
	minIntensity = *minEigenvalue = 1e12;

	globalVec[0] = exp(1.0f);
	globalVec[1] = (float)M_PI;
	globalVec[2] = exp(2.0f);
	float norm = getNorm(globalVec);
	globalVec[0] = globalVec[0]/norm;
	globalVec[1] = globalVec[1]/norm;
	globalVec[2] = globalVec[2]/norm;
	thresh = cos(20 * (float)M_PI / 180);

	gridPoints = NULL;
	gridScalars = NULL;
	gridShowPositions = NULL;
	gridV1s = NULL;
	gridV2s = NULL;
	gridV3s = NULL;
	gridGradients = NULL;
	edgePoints = NULL;
	facePoints = NULL;
	cellPoints = NULL;
	edges = NULL;
	faces = NULL;
	cells = NULL;
	grid_2DGradMag = NULL;
	storeGrid_2DGradMag = NULL;
	segments = NULL;
	quads = NULL;
	points = NULL;
	vertices = NULL;
}

General_Data::~General_Data(void)
{
	delete[] gridPoints;
	delete[] gridScalars;
	delete[] gridShowPositions;
	delete[] gridV1s;
	delete[] gridV2s;
	delete[] gridV3s;
	delete[] gridGradients;
	delete[] edgePoints;
	delete[] facePoints;
	delete[] cellPoints;
	delete[] edges;
	delete[] faces;
	delete[] cells;
	delete[] grid_2DGradMag;
	delete[] storeGrid_2DGradMag;
	delete segments;
	delete quads;
	delete points;
	delete vertices;
}

float *General_Data::getGridScalars()
{
	if (!storeGridScalar)
	{
		computeGridScalars();
		storeGridScalar = true;
	}
	return gridScalars;
}

float *General_Data::getGridV1s()
{
	if (!storeGridVector)
	{
		computeGridVectors();
		storeGridVector = true;
	}
	return gridV1s;
}

float *General_Data::getGridV2s()
{
	if (!storeGridVector)
	{
		computeGridVectors();
		storeGridVector = true;
	}
	return gridV2s;
}

float *General_Data::getGridV3s()
{
	if (!storeGridVector)
	{
		computeGridVectors();
		storeGridVector = true;
	}
	return gridV3s;
}

float *General_Data::getGridGradients()
{
	if (!storeGridGradient)
	{
		computeGridGradients();
		storeGridGradient = true;
	}
	return gridGradients;
}

float *General_Data::getGridShowPositions()
{
	if (!storeGridShowPos)
	{
		computeGridShowPositions();
		storeGridShowPos = true;
	}
	return gridShowPositions;
}

Face *General_Data::getFaces()
{
	return faces;
}

float *General_Data::getEdgePoints()
{
	if (!storeEdgePoint)
	{
		computeEdgePoints();
		storeEdgePoint = true;
	}
	return edgePoints;
}

float *General_Data::getFacePoints()
{
	if (!storeFacePoint)
	{
		computeFacePoints();
		storeFacePoint = true;
	}
	return facePoints;
}

float *General_Data::getCellPoints()
{
	if (!storeCellPoint)
	{
		computeCellPoints();
		storeCellPoint = true;
	}
	return cellPoints;
}

vector<Segment> *General_Data::getSegments()
{
	return segments;
}

vector<Quad> *General_Data::getQuads()
{
	return quads;
}

vector<Point> *General_Data::getPoints()
{
	return points;
}

vector<Vertex> *General_Data::getVertices()
{
	return vertices;
}

void General_Data::getGridPointPos(int index, float v[3])
{
	v[0] = gridPoints[3*index];
	v[1] = gridPoints[3*index+1];
	v[2] = gridPoints[3*index+2];
}

void General_Data::getGridPointPosD(int index, float v[3])
{
	v[0] = gridPoints[3*index];
	v[1] = gridPoints[3*index+1];
	v[2] = gridPoints[3*index+2];
}

void General_Data::computeGridShowPositions()
{
	gridShowPositions = new float[ 3 * maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		float position[3];
		float showPos[3];
		getGridPointPosD(index, position);
		getShowPos(position, showPos);
		gridShowPositions[3*index] = showPos[0];
		gridShowPositions[3*index+1] = showPos[1];
		gridShowPositions[3*index+2] = showPos[2];
	}
}

void General_Data::computeGridScalars()
{
	gridScalars = new float[ maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		float position[3];
		float scalar;
		getGridPointPosD(index, position);
		getScalarCubic(position, &scalar);
		gridScalars[index] = ((scalar - minScalar)/(maxScalar-minScalar));
	}
}

void General_Data::computeGridVectors()
{
	gridV1s = new float[ 3 * maxGridIndex ] ;
	gridV2s = new float[ 3 * maxGridIndex ] ;
	gridV3s = new float[ 3 * maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		float position[3];
		float v1[3], v2[3], v3[3];
		getGridPointPosD(index, position);
		SelfAdjointEigenSolver<Matrix3f> eigensolver;
		getEigensolverCubic(position, &eigensolver);
		getV1(&eigensolver, v1);
		getV2(&eigensolver, v2);
		getV3(&eigensolver, v3);
		gridV1s[3*index] = v1[0] / (float)maxGrid / 4;
		gridV1s[3*index+1] = v1[1] / (float)maxGrid / 4;
		gridV1s[3*index+2] = v1[2] / (float)maxGrid / 4;
		gridV2s[3*index] = v2[0] / (float)maxGrid / 4;
		gridV2s[3*index+1] = v2[1] / (float)maxGrid / 4;
		gridV2s[3*index+2] = v2[2] / (float)maxGrid / 4;
		gridV3s[3*index] = v3[0] / (float)maxGrid / 4;
		gridV3s[3*index+1] = v3[1] / (float)maxGrid / 4;
		gridV3s[3*index+2] = v3[2] / (float)maxGrid / 4;
	}
}

void General_Data::computeGridGradients()
{
	gridGradients = new float[ 3 * maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		float position[3];
		float gradient[3];
		getGridPointPosD(index, position);
		getGradient(position, gradient);
		gridGradients[3*index] = gradient[0] / maxGradientNorm / (float)maxGrid / 4;
		gridGradients[3*index+1] = gradient[1] / maxGradientNorm / (float)maxGrid / 4;
		gridGradients[3*index+2] = gradient[2] / maxGradientNorm / (float)maxGrid / 4;
	}
}

void General_Data::computeEdgePoints()
{
	edgePoints = new float[ 3 * (*totalEdgePoints) ] ;
	int count = 0;
	for ( int index = 0 ; index < maxEdgeIndex && count < (*totalEdgePoints) ; index ++ )
	{
		Edge *edge = &(edges[index]);
		if (edge->extremal)
		{
			float showPos[3];
			getShowPos(edge->edgePoint, showPos);
			edgePoints[3*count] = showPos[0];
			edgePoints[3*count+1] = showPos[1];
			edgePoints[3*count+2] = showPos[2];
			count ++;
		}
	}
}

void General_Data::computeFacePoints()
{
	facePoints = new float[ 3 * (*totalFacePoints) ] ;
	int count = 0;
	for ( int index = 0 ; index < maxFaceIndex && count < (*totalFacePoints) ; index ++ )
	{
		Face *face = &(faces[index]);
		if (face->extremal)
		{
			float showPos[3];
			getShowPos(face->facePoint, showPos);
			facePoints[3*count] = showPos[0];
			facePoints[3*count+1] = showPos[1];
			facePoints[3*count+2] = showPos[2];
			count ++;
		}
	}
}

void General_Data::computeCellPoints()
{
	cellPoints = new float[ 3 * (*totalCellPoints) ] ;
	int count = 0;
	for ( int index = 0 ; index < maxCellIndex && count < (*totalCellPoints) ; index ++ )
	{
		Cell *cell = &(cells[index]);
		if (cell->extremal)
		{
			float showPos[3];
			getShowPos(cell->centroid, showPos);
			cellPoints[3*count] = showPos[0];
			cellPoints[3*count+1] = showPos[1];
			cellPoints[3*count+2] = showPos[2];
			count ++;
		}
	}
}

void General_Data::getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])
{
	Vector3f ev;
	switch (dataType)
	{
	case 1:
		ev = eigensolver->eigenvectors().col(2).normalized();
		break;
	case 2:
		ev = eigensolver->eigenvectors().col(0).normalized();
		break;
	}
	v[0] = (float)ev(0);
	v[1] = (float)ev(1);
	v[2] = (float)ev(2);
}

void General_Data::getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])
{
	Vector3f ev = eigensolver->eigenvectors().col(1).normalized();
	v[0] = (float)ev(0);
	v[1] = (float)ev(1);
	v[2] = (float)ev(2);
}

void General_Data::getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])
{
	Vector3f ev;
	switch (dataType)
	{
	case 1:
		ev = eigensolver->eigenvectors().col(0).normalized();
		break;
	case 2:
		ev = eigensolver->eigenvectors().col(2).normalized();
		break;
	}
	v[0] = (float)ev(0);
	v[1] = (float)ev(1);
	v[2] = (float)ev(2);
}

float General_Data::getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3])
{
	Vector3f ev = eigensolver->eigenvalues();
	if ( ev(2) > 0 )
	{
		relativeSaliencies[0] = float((ev(2) - ev(1)) / ev(2));
		relativeSaliencies[1] = float((ev(1) - ev(0)) / ev(2));
		relativeSaliencies[2] = float(ev(0) / ev(2));
		if ((*minEigenvalue) > ev(2)) (*minEigenvalue) = ev(2);
		if ((*maxEigenvalue) < ev(2)) (*maxEigenvalue) = ev(2);
	}
	else
	{
		relativeSaliencies[0] = 0;
		relativeSaliencies[1] = 0;
		relativeSaliencies[2] = 0;
	}
	return float(ev(2));
}

bool General_Data::getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver)
{
	int count = 0;
	float epsilon = 1e-12;
	float tensorV[6];
	getTensorGP(index, tensorV);
	for ( int i = 0; i < 6; i ++ )
	{
		if ( fabs(tensorV[i]) < epsilon )
		{
			count ++;
			tensorV[i] = epsilon;
		}
	}
	Matrix3f tensorM;
	tensorM << tensorV[0], tensorV[1], tensorV[2],
			   tensorV[1], tensorV[3], tensorV[4],
			   tensorV[2], tensorV[4], tensorV[5];
	eigensolver->compute(tensorM);
	if ( count > 1 ) return false;
	return true;
}

bool General_Data::getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver)
{
	int count = 0;
	float epsilon = 1e-12;
	float tensorV[6];
	getTensor(position, tensorV);
	for ( int i = 0; i < 6; i ++ )
	{
		if ( fabs(tensorV[i]) < epsilon )
		{
			count ++;
			tensorV[i] = epsilon;
		}
	}
	Matrix3f tensorM;
	tensorM << tensorV[0], tensorV[1], tensorV[2],
			   tensorV[1], tensorV[3], tensorV[4],
			   tensorV[2], tensorV[4], tensorV[5];
	eigensolver->compute(tensorM);
	if ( count > 1 ) return false;
	return true;
}

bool General_Data::getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver)
{
	int count = 0;
	float epsilon = 1e-12;
	float tensorV[6];
	getTensorCubic(position, tensorV);
	for ( int i = 0; i < 6; i ++ )
	{
		if ( fabs(tensorV[i]) < epsilon )
		{
			count ++;
			tensorV[i] = epsilon;
		}
	}
	Matrix3f tensorM;
	tensorM << tensorV[0], tensorV[1], tensorV[2],
			   tensorV[1], tensorV[3], tensorV[4],
			   tensorV[2], tensorV[4], tensorV[5];
	eigensolver->compute(tensorM);
	if ( count > 1 ) return false;
	return true;
}

bool General_Data::getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver)
{
	int count = 0;
	float epsilon = 1e-12;
	float tensorV[6];
	getTensorCubicTable(axis, i, j, k, index, tensorV);
	for ( int i = 0; i < 6; i ++ )
	{
		if ( fabs(tensorV[i]) < epsilon )
		{
			count ++;
			tensorV[i] = epsilon;
		}
	}
	Matrix3f tensorM;
	tensorM << tensorV[0], tensorV[1], tensorV[2],
			   tensorV[1], tensorV[3], tensorV[4],
			   tensorV[2], tensorV[4], tensorV[5];
	eigensolver->compute(tensorM);
	if ( count > 1 ) return false;
	return true;
}

bool General_Data::getEigensolver2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2)
{
	float tensorV1[6], tensorV2[6];
	getTensor(p1, tensorV1);
	getTensor(p2, tensorV2);
	float epsilon = 1e-12;
	bool same = true;
	for ( int i = 0; i < 6; i ++ )
	{
		float abs1 = fabs(tensorV1[i]);
		float abs2 = fabs(tensorV2[i]);
		float absDiff = fabs(abs1-abs2);
		if (absDiff/abs1 > epsilon ||absDiff/abs2 > epsilon)
		{
			same = false;
			break;
		}
	}
	if (!same)
	{
		Matrix3f tensorM1, tensorM2;
		tensorM1 << tensorV1[0], tensorV1[1], tensorV1[2],
				   tensorV1[1], tensorV1[3], tensorV1[4],
				   tensorV1[2], tensorV1[4], tensorV1[5];
		tensorM2 << tensorV2[0], tensorV2[1], tensorV2[2],
				   tensorV2[1], tensorV2[3], tensorV2[4],
				   tensorV2[2], tensorV2[4], tensorV2[5];
		e1->compute(tensorM1);
		e2->compute(tensorM2);
	}
	return same;
}

bool General_Data::getEigensolverCubic2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2)
{
	float tensorV1[6], tensorV2[6];
	getTensorCubic(p1, tensorV1);
	getTensorCubic(p2, tensorV2);
	float epsilon = 1e-12;
	bool same = true;
	for ( int i = 0; i < 6; i ++ )
	{
		float abs1 = fabs(tensorV1[i]);
		float abs2 = fabs(tensorV2[i]);
		float absDiff = fabs(abs1-abs2);
		if (absDiff/abs1 > epsilon ||absDiff/abs2 > epsilon)
		{
			same = false;
			break;
		}
	}
	if (!same)
	{
		Matrix3f tensorM1, tensorM2;
		tensorM1 << tensorV1[0], tensorV1[1], tensorV1[2],
				   tensorV1[1], tensorV1[3], tensorV1[4],
				   tensorV1[2], tensorV1[4], tensorV1[5];
		tensorM2 << tensorV2[0], tensorV2[1], tensorV2[2],
				   tensorV2[1], tensorV2[3], tensorV2[4],
				   tensorV2[2], tensorV2[4], tensorV2[5];
		e1->compute(tensorM1);
		e2->compute(tensorM2);
	}
	return same;
}

int General_Data::getMaxGridIndex()
{
	return maxGridIndex;
}

float General_Data::getMaxGradientNorm()
{
	return maxGradientNorm;
}

float General_Data::getMaxGrid()
{
	return maxGrid;
}

void General_Data::getShowPos(float position[3], float showPos[3])
{
	showPos[0] = ((float)position[0] - 2 - rightBackTop[0]) / halfDataSize;
	showPos[1] = ((float)position[1] - 2 - rightBackTop[1]) / halfDataSize;
	showPos[2] = ((float)position[2] - 2 - rightBackTop[2]) / halfDataSize * thickness;
}

void General_Data::buildSurface()
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Build Surface Begin:" << endl;
	quads = new vector<Quad>;
	for ( int i = 0 ; i < gridx - 1 ; i ++ )
		for ( int j = 0 ; j < gridy ; j ++ )
			for ( int k = 0 ; k < gridz ; k ++ )
			{
				Edge *edge = &(edges[getEdgeIndex(1, i, j, k)]);
				if ( !edge->valid ) continue;
				if ( !edge->extremal || j == 0 || j == gridy-1 || k == 0 || k == gridz-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i, j-1, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				Cell *cell3 = &(cells[getIndex(1, 0, i, j-1, k-1, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell3->valid ) continue;
				Cell *cell4 = &(cells[getIndex(1, 0, i, j, k-1, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell4->valid ) continue;
				buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
			}
	for ( int i = 0 ; i < gridx ; i ++ )
		for ( int j = 0 ; j < gridy - 1 ; j ++ )
			for ( int k = 0 ; k < gridz ; k ++ )
			{
				Edge *edge = &(edges[getEdgeIndex(2, i, j, k)]);
				if ( !edge->valid ) continue;
				if ( !edge->extremal || i == 0 || i == gridx-1 || k == 0 || k == gridz-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i-1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				Cell *cell3 = &(cells[getIndex(1, 0, i-1, j, k-1, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell3->valid ) continue;
				Cell *cell4 = &(cells[getIndex(1, 0, i, j, k-1, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell4->valid ) continue;
				buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
			}
	for ( int i = 0 ; i < gridx ; i ++ )
		for ( int j = 0 ; j < gridy ; j ++ )
			for ( int k = 0 ; k < gridz - 1 ; k ++ )
			{
				Edge *edge = &(edges[getEdgeIndex(3, i, j, k)]);
				if ( !edge->valid ) continue;
				if ( !edge->extremal || i == 0 || i == gridx-1 || j == 0 || j == gridy-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i, j-1, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				Cell *cell3 = &(cells[getIndex(1, 0, i-1, j-1, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell3->valid ) continue;
				Cell *cell4 = &(cells[getIndex(1, 0, i-1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell4->valid ) continue;
				buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
			}
	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::buildSurfaceIteration(Edge *edge, Cell *cell1, Cell *cell2, Cell *cell3, Cell *cell4)
{
	Quad quad;
	if ( cell1->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell1->centroid, vertex.position);
		quad.vertIdxs[0] = vertices->size();
		cell1->vertIdx = quad.vertIdxs[0];
		vertices->push_back(vertex);
	}
	else
	{
		quad.vertIdxs[0] = cell1->vertIdx;
	}
	if ( cell2->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell2->centroid, vertex.position);
		quad.vertIdxs[1] = vertices->size();
		cell2->vertIdx = quad.vertIdxs[1];
		vertices->push_back(vertex);
	}
	else
	{
		quad.vertIdxs[1] = cell2->vertIdx;
	}
	if ( cell3->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell3->centroid, vertex.position);
		quad.vertIdxs[2] = vertices->size();
		cell3->vertIdx = quad.vertIdxs[2];
		vertices->push_back(vertex);
	}
	else
	{
		quad.vertIdxs[2] = cell3->vertIdx;
	}
	if ( cell4->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell4->centroid, vertex.position);
		quad.vertIdxs[3] = vertices->size();
		cell4->vertIdx = quad.vertIdxs[3];
		vertices->push_back(vertex);
	}
	else
	{
		quad.vertIdxs[3] = cell4->vertIdx;
	}
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	getEigensolverCubic(edge->edgePoint, &eigensolver);
	quad.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, quad.relativeSaliencies);
	getScalarCubic(edge->edgePoint, &(quad.localIntensity));
	quad.type = edge->edgeTag;
	quads->push_back(quad);
}

void General_Data::buildPointIteration(Cell *cell)
{
	Vertex vertex;
	Point point;
	getShowPos(cell->extremalPoint, vertex.position);
	point.vertIdx = vertices->size();
	vertices->push_back(vertex);
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	getEigensolverCubic(cell->extremalPoint, &eigensolver);
	point.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, point.relativeSaliencies);
	getScalarCubic(cell->extremalPoint, &(point.localIntensity));
	point.type = cell->extrPType;
	points->push_back(point);
}

void General_Data::buildCurve()
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Build Curve Begin:" << endl;
	segments = new vector<Segment>;
	for ( int i = 0 ; i < gridx ; i ++ )
		for ( int j = 0 ; j < gridy - 1 ; j ++ )
			for ( int k = 0 ; k < gridz - 1 ; k ++ )
			{
				Face *face = &(faces[getFaceIndex(1, i, j, k)]);
				if ( !face->valid ) continue;
				if ( !face->extremal || i == 0 || i == gridx-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i - 1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				buildCurveIteration(face, cell1, cell2);
			}
	for ( int i = 0 ; i < gridx - 1 ; i ++ )
		for ( int j = 0 ; j < gridy ; j ++ )
			for ( int k = 0 ; k < gridz - 1 ; k ++ )
			{
				Face *face = &(faces[getFaceIndex(2, i, j, k)]);
				if ( !face->valid ) continue;
				if ( !face->extremal || j == 0 || j == gridy-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i, j - 1, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				buildCurveIteration(face, cell1, cell2);
			}
	for ( int i = 0 ; i < gridx - 1 ; i ++ )
		for ( int j = 0 ; j < gridy - 1 ; j ++ )
			for ( int k = 0 ; k < gridz ; k ++ )
			{
				Face *face = &(faces[getFaceIndex(3, i, j, k)]);
				if ( !face->valid ) continue;
				if ( !face->extremal || k == 0 || k == gridz-1 ) continue;
				Cell *cell1 = &(cells[getIndex(1, 0, i, j, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell1->valid ) continue;
				Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
				if ( !cell2->valid ) continue;
				buildCurveIteration(face, cell1, cell2);
			}
	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::buildCurveIteration(Face *face, Cell *cell1, Cell *cell2)
{
	Segment segment;
	if ( cell1->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell1->centroid, vertex.position);
		segment.vertIdxs[0] = vertices->size();
		cell1->vertIdx = segment.vertIdxs[0];
		vertices->push_back(vertex);
	}
	else
	{
		segment.vertIdxs[0] = cell1->vertIdx;
	}
	if ( cell2->vertIdx < 0 )
	{
		Vertex vertex;
		getShowPos(cell2->centroid, vertex.position);
		segment.vertIdxs[1] = vertices->size();
		cell2->vertIdx = segment.vertIdxs[1];
		vertices->push_back(vertex);
	}
	else
	{
		segment.vertIdxs[1] = cell2->vertIdx;
	}
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	getEigensolverCubic(face->facePoint, &eigensolver);
	segment.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, segment.relativeSaliencies);
	float localI;
	getScalarCubic(face->facePoint, &localI);
	if (minIntensity > localI) minIntensity = localI;
	if (maxIntensity < localI) maxIntensity = localI;
	segment.localIntensity = localI;
	segment.type = face->faceTag;
	segments->push_back(segment);
}

void General_Data::cellPhase(float* scalars,float* tensors,float* gradients)
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Cell Phase Begin:" << endl;
	points = new vector<Point>;
	vertices = new vector<Vertex>;
	conpoints = new concurrent_vector<Point>;
	convertices = new concurrent_vector<Vertex>;
	change = 0.017f * gridSize;
	dispCellX = (int)floor( gridx / 2.0 ) - 1;
	dispCellY = (int)floor( gridy / 2.0 ) - 1;
	dispCellZ = (int)floor( gridz / 2.0 ) - 1;
	maxCellIndex = (gridx-1) * (gridy-1) * (gridz-1);
	cells = new Cell[ maxCellIndex ] ;
	totalCellPoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	ParallelCell parallel_cell(gridx,gridy,gridz,totalCellPoints,edges,faces,cells,conpoints,dataType,gridPoints,change,sizex,sizey,sizez,halfSize,scalars,tensors,gradients,gridSize,maxEigenvalue, minEigenvalue,convertices,halfDataSize,rightBackTop,thickness);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_cell,auto_partitioner());

	concurrent_vector<Point>::iterator iter = conpoints->begin();
	while( iter != conpoints->end())
	{
		Point p;
		p.vertIdx = (*iter).vertIdx;
	    p.relativeSaliencies[0] = (*iter).relativeSaliencies[0];
		p.relativeSaliencies[1] = (*iter).relativeSaliencies[1];
		p.relativeSaliencies[2] = (*iter).relativeSaliencies[2];
	    p.localIntensity = (*iter).localIntensity;
		p.firstEigenvalue = (*iter).firstEigenvalue;
		p.type=(*iter).type;		
		points->push_back(p);
		iter++;
	}

	concurrent_vector<Vertex>::iterator iterv = convertices->begin();
	while( iterv != convertices->end())
	{
		Vertex p;
		p.position[0]=(*iterv).position[0];
		p.position[1]=(*iterv).position[1];
		p.position[2]=(*iterv).position[2];
		vertices->push_back(p);
		iterv++;
	}
	timeEnd = clock();
	//t_cellPhase_other = timeEnd-timeStart - t_cellPhase_extremalPoint;
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	/*cout << "\tt_cellPhase_other: " << t_cellPhase_other << endl;
	cout << "\tt_cellPhase_extremalPoint: " << t_cellPhase_extremalPoint << endl;*/
	//cout << "\tt_cellPhase_interp: " << t_cellPhase_interp/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::cellPhaseIteration(int i, int j, int k)
{
	Cell *cell = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
	if (!getCentroid(i, j, k, cell)) return;
	cell->valid = true;
	centroidProjection(i, j, k, cell);
	if (cell->extrP) buildPointIteration(cell);
}

bool General_Data::getCentroid(int i, int j, int k, Cell *cell)
{
	Face *cellFaces[6];
	cellFaces[0] = &(faces[getFaceIndex(1, i, j, k)]);
	if (!cellFaces[0]->valid) return false;
	cellFaces[1] = &(faces[getFaceIndex(1, i + 1, j, k)]);
	if (!cellFaces[1]->valid) return false;
	cellFaces[2] = &(faces[getFaceIndex(2, i, j, k)]);
	if (!cellFaces[2]->valid) return false;
	cellFaces[3] = &(faces[getFaceIndex(2, i, j + 1, k)]);
	if (!cellFaces[3]->valid) return false;
	cellFaces[4] = &(faces[getFaceIndex(3, i, j, k)]);
	if (!cellFaces[4]->valid) return false;
	cellFaces[5] = &(faces[getFaceIndex(3, i, j, k + 1)]);
	if (!cellFaces[5]->valid) return false;
	Edge *cellEdges[12];
	cellEdges[0] = &(edges[getEdgeIndex(1, i, j, k)]);
	cellEdges[1] = &(edges[getEdgeIndex(1, i, j + 1, k)]);
	cellEdges[2] = &(edges[getEdgeIndex(1, i, j + 1, k + 1)]);
	cellEdges[3] = &(edges[getEdgeIndex(1, i, j, k + 1)]);
	cellEdges[4] = &(edges[getEdgeIndex(2, i, j, k)]);
	cellEdges[5] = &(edges[getEdgeIndex(2, i + 1, j, k)]);
	cellEdges[6] = &(edges[getEdgeIndex(2, i + 1, j, k + 1)]);
	cellEdges[7] = &(edges[getEdgeIndex(2, i, j, k + 1)]);
	cellEdges[8] = &(edges[getEdgeIndex(3, i, j, k)]);
	cellEdges[9] = &(edges[getEdgeIndex(3, i + 1, j, k)]);
	cellEdges[10] = &(edges[getEdgeIndex(3, i + 1, j + 1, k)]);
	cellEdges[11] = &(edges[getEdgeIndex(3, i, j + 1, k)]);

	int ci = i;
	int cj = j;
	int ck = k;

	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	float area = 0;
	for ( int i = 0; i < 3; i ++ )
	{
		area -= cellFaces[2*i]->spherArea;
		area += cellFaces[2*i+1]->spherArea;
	}
	//cell->area = area;
	float wn = area / ((float)(4*M_PI));
	wn = floor( wn + 0.5f );
	if ( (int)wn % 2 != 0 )
	{
		cell->extrP = true;

		float **grads;
		grads = new float* [8] ;
		for ( int i = 0 ; i < 8 ; i ++ )
		{
			grads[ i ] = new float[ 3 ] ;
		}
		float verts[8][3];
		for ( int i = 0 ; i < 2 ; i ++ )
			for ( int j = 0 ; j < 2 ; j ++ )
				for ( int k = 0 ; k < 2 ; k ++ )
				{
					int index = getIndex(1, 0, ci+i, cj+j, ck+k, gridx, gridy, gridz);
					getGridPointPosD(index, verts[i+2*j+4*k]);
					index = getIndex(3, 0, ci+i+2, cj+j+2, ck+k+2, sizex, sizey, sizez);
					getGradientGP(index, grads[i+2*j+4*k]);
				}
		int **tris;
		tris = new int* [12] ;
		for ( int i = 0 ; i < 12 ; i ++ )
		{
			tris[ i ] = new int[ 3 ] ;
		}
		int static_tris[12][3] = 
			{{0, 2, 3}, {0, 3, 1},
			{5, 7, 6}, {5, 6, 4},
			{4, 6, 2}, {4, 2, 0},
			{1, 3, 7}, {1, 7, 5},
			{4, 0, 1}, {4, 1, 5},
			{2, 6, 7}, {2, 7, 3}};
		for ( int i = 0 ; i < 12 ; i ++ )
			for ( int j = 0 ; j < 3 ; j ++ )
				tris[i][j] = static_tris[i][j];

		MVC mvc(8, 12, grads, tris);
		float weights[8];
		float origin[3] = {};
		mvc.getWeights(origin, weights);
		float p[3];
		p[0] = 0;
		p[1] = 0;
		p[2] = 0;
		for ( int i = 0 ; i < 8 ; i ++ )
			for ( int j = 0 ; j < 3 ; j ++ )
				p[j] += verts[i][j]*weights[i];
		if ( p[0] < verts[0][0] )
		{
			p[0] = verts[0][0];
		}
		if ( p[0] > verts[7][0] )
		{
			p[0] = verts[7][0];
		}
		if ( p[1] < verts[0][1] )
		{
			p[1] = verts[0][1];
		}
		if ( p[1] > verts[7][1] )
		{
			p[1] = verts[7][1];
		}
		if ( p[2] < verts[0][2] )
		{
			p[2] = verts[0][2];
		}
		if ( p[2] > verts[7][2] )
		{
			p[2] = verts[7][2];
		}
		for ( int i = 0 ; i < 3 ; i ++ )
			cell->extremalPoint[i] = p[i];

		float p10[3], p11[3], p13[3], p14[3];
		p10[0] = p[0]-2*change;
		p10[1] = p[1];
		p10[2] = p[2];
		p11[0] = p[0]-change;
		p11[1] = p[1];
		p11[2] = p[2];
		p13[0] = p[0]+change;
		p13[1] = p[1];
		p13[2] = p[2];
		p14[0] = p[0]+2*change;
		p14[1] = p[1];
		p14[2] = p[2];

		float p20[3], p21[3], p23[3], p24[3];
		p20[0] = p[0];
		p20[1] = p[1]-2*change;
		p20[2] = p[2];
		p21[0] = p[0];
		p21[1] = p[1]-change;
		p21[2] = p[2];
		p23[0] = p[0];
		p23[1] = p[1]+change;
		p23[2] = p[2];
		p24[0] = p[0];
		p24[1] = p[1]+2*change;
		p24[2] = p[2];

		float p30[3], p31[3], p33[3], p34[3];
		p30[0] = p[0];
		p30[1] = p[1];
		p30[2] = p[2]-2*change;
		p31[0] = p[0];
		p31[1] = p[1];
		p31[2] = p[2]-change;
		p33[0] = p[0];
		p33[1] = p[1];
		p33[2] = p[2]+change;
		p34[0] = p[0];
		p34[1] = p[1];
		p34[2] = p[2]+2*change;

		float f, f10, f11, f13, f14, f20, f21, f23, f24, f30, f31, f33, f34;

		clock_t timeStart_i, timeEnd_i;
		timeStart_i = clock();
		getScalarCubic(p, &f);
		getScalarCubic(p10, &f10);
		getScalarCubic(p11, &f11);
		getScalarCubic(p13, &f13);
		getScalarCubic(p14, &f14);
		getScalarCubic(p20, &f20);
		getScalarCubic(p21, &f21);
		getScalarCubic(p23, &f23);
		getScalarCubic(p24, &f24);
		getScalarCubic(p30, &f30);
		getScalarCubic(p31, &f31);
		getScalarCubic(p33, &f33);
		getScalarCubic(p34, &f34);
		timeEnd_i = clock();
		//t_cellPhase_interp = t_cellPhase_interp + timeEnd_i - timeStart_i;

		float secondDiff1 = -f10+16*f11-30*f+16*f13-f14;
		float secondDiff2 = -f20+16*f21-30*f+16*f23-f24;
		float secondDiff3 = -f30+16*f31-30*f+16*f33-f34;
		if (secondDiff1<0 && secondDiff2<0 && secondDiff3<0) cell->extrPType = 1;
		else if (secondDiff1>0 && secondDiff2>0 && secondDiff3>0) cell->extrPType = 2;
		else cell->extrPType = 3;
	}

	/*timeEnd = clock();
	t_cellPhase_extremalPoint = t_cellPhase_extremalPoint + timeEnd - timeStart;*/

	cell->centroid[0] = 0;
	cell->centroid[1] = 0;
	cell->centroid[2] = 0;
	int count = 0;
	bool minEP, maxEP, minFP, maxFP, other;
	minEP = maxEP = minFP = maxFP = other = false;

	for ( int i = 0; i < 12; i ++ )
	{
		Edge *edge = cellEdges[i];
		if ( edge->extremal )
		{
			cell->extremal = true;
			if ( edge->edgeTag == 1 ) maxEP = true;
			else minEP = true;
			cell->centroid[0] += edge->edgePoint[0];
			cell->centroid[1] += edge->edgePoint[1];
			cell->centroid[2] += edge->edgePoint[2];
			count ++;
		}
	}

	for ( int i = 0; i < 6; i ++ )
	{
		Face *face = cellFaces[i];
		if ( face->extremal )
		{
			cell->extremal = true;
			if ( face->faceTag == 1 ) maxFP = true;
			else if ( face->faceTag == 2 ) minFP = true;
			else other = true;
			cell->centroid[0] += face->facePoint[0];
			cell->centroid[1] += face->facePoint[1];
			cell->centroid[2] += face->facePoint[2];
			count ++;
		}
	}

	if ( cell->extremal )
	{
		cell->centroid[0] /= (float)count;
		cell->centroid[1] /= (float)count;
		cell->centroid[2] /= (float)count;
		(*totalCellPoints) ++;
	}
	if ( !other && !( (minEP || minFP) && (maxEP || maxFP) ) )
	{
		if (minFP) cell->projType = 3;
		else if (maxFP) cell->projType = 4;
		else if (minEP) cell->projType = 1;
		else cell->projType = 2;
	}
	return true;
}

void General_Data::centroidProjection(int i, int j, int k, Cell *cell)
{
	if ( cell->extremal && cell->projType > 0 )
	{
		float border1[3], border2[3];
		getGridPointPosD(getIndex(1, 0, i, j, k, gridx, gridy, gridz), border1);
		getGridPointPosD(getIndex(1, 0, i + 1, j + 1, k + 1, gridx, gridy, gridz), border2);
		float projThresh = 0.001f * gridSize;
		float stepSize = 0.1f * gridSize;
		float p[3];
		p[0] = cell->centroid[0];
		p[1] = cell->centroid[1];
		p[2] = cell->centroid[2];
		int type = cell->projType;
		float q[3];
		getProjDirection(type, p, q);
		float qNorm = getNorm(q);
		int count = 0;
		float p1[3];
		float q1[3];
		float q1Norm;
		while ( qNorm > projThresh && count < 10 )
		{
			count ++;
			p1[0] = p[0] - stepSize * ( q[0] / qNorm );
			p1[1] = p[1] - stepSize * ( q[1] / qNorm );
			p1[2] = p[2] - stepSize * ( q[2] / qNorm );
			int dCount = 0;
			if ( p1[0] < border1[0] )
			{
				p1[0] = border1[0];
				dCount ++;
			}
			if ( p1[0] > border2[0] )
			{
				p1[0] = border2[0];
				dCount ++;
			}
			if ( p1[1] < border1[1] )
			{
				p1[1] = border1[1];
				dCount ++;
			}
			if ( p1[1] > border2[1] )
			{
				p1[1] = border2[1];
				dCount ++;
			}
			if ( p1[2] < border1[2] )
			{
				p1[2] = border1[2];
				dCount ++;
			}
			if ( p1[2] > border2[2] )
			{
				p1[2] = border2[2];
				dCount ++;
			}
			getProjDirection(type, p1, q1);
			q1Norm = getNorm(q1);
			if ( q1Norm > qNorm )
			{
				stepSize /= 2;
				continue;
			}
			p[0] = p1[0];
			p[1] = p1[1];
			p[2] = p1[2];
			if ( dCount == 3 ) break;
			q[0] = q1[0];
			q[1] = q1[1];
			q[2] = q1[2];
			qNorm = q1Norm;
		}
		cell->centroid[0] = p[0];
		cell->centroid[1] = p[1];
		cell->centroid[2] = p[2];
	}
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	if ( cell->extrP && ( cell->extrPType == 1 || cell->extrPType == 2 ) )
	{
		float border1[3], border2[3];
		getGridPointPosD(getIndex(1, 0, i, j, k, gridx, gridy, gridz), border1);
		getGridPointPosD(getIndex(1, 0, i + 1, j + 1, k + 1, gridx, gridy, gridz), border2);
		float projThresh = 0.001f * gridSize;
		float stepSize = 0.01f * gridSize;
		if (cell->extrPType == 1) stepSize = -stepSize;
		float p[3];
		p[0] = cell->extremalPoint[0];
		p[1] = cell->extremalPoint[1];
		p[2] = cell->extremalPoint[2];
		float q[3];
		clock_t timeStart_i, timeEnd_i;
		timeStart_i = clock();
		getGradientCubic(p, q);
		timeEnd_i = clock();
		//t_cellPhase_interp = t_cellPhase_interp + timeEnd_i - timeStart_i;
		float qNorm = getNorm(q);
		int count = 0;
		float p1[3];
		float q1[3];
		float q1Norm;
		while ( qNorm > projThresh && count < 10 )
		{
			count ++;
			p1[0] = p[0] - stepSize * ( q[0] / qNorm );
			p1[1] = p[1] - stepSize * ( q[1] / qNorm );
			p1[2] = p[2] - stepSize * ( q[2] / qNorm );
			int dCount = 0;
			if ( p1[0] < border1[0] )
			{
				p1[0] = border1[0];
				dCount ++;
			}
			if ( p1[0] > border2[0] )
			{
				p1[0] = border2[0];
				dCount ++;
			}
			if ( p1[1] < border1[1] )
			{
				p1[1] = border1[1];
				dCount ++;
			}
			if ( p1[1] > border2[1] )
			{
				p1[1] = border2[1];
				dCount ++;
			}
			if ( p1[2] < border1[2] )
			{
				p1[2] = border1[2];
				dCount ++;
			}
			if ( p1[2] > border2[2] )
			{
				p1[2] = border2[2];
				dCount ++;
			}
			timeStart_i = clock();
			getGradientCubic(p1, q1);
			timeEnd_i = clock();
			//t_cellPhase_interp = t_cellPhase_interp + timeEnd_i - timeStart_i;
			q1Norm = getNorm(q1);
			if ( q1Norm > qNorm )
			{
				stepSize /= 10;
				continue;
			}
			p[0] = p1[0];
			p[1] = p1[1];
			p[2] = p1[2];
			if ( dCount == 3 ) break;
			q[0] = q1[0];
			q[1] = q1[1];
			q[2] = q1[2];
			qNorm = q1Norm;
		}
		cell->extremalPoint[0] = p[0];
		cell->extremalPoint[1] = p[1];
		cell->extremalPoint[2] = p[2];
	}
	/*timeEnd = clock();
	t_cellPhase_extremalPoint = t_cellPhase_extremalPoint + timeEnd - timeStart;*/
}

void General_Data::getProjDirection(int type, float position[3], float q[3])
{
	float g[3];
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	clock_t timeStart_i, timeEnd_i;
	timeStart_i = clock();
	getGradientCubic(position, g);
	getEigensolverCubic(position, &eigensolver);
	timeEnd_i = clock();
	//t_cellPhase_interp = t_cellPhase_interp + timeEnd_i - timeStart_i;
	float v1[3], v3[3];
	getV1(&eigensolver, v1);
	getV3(&eigensolver, v3);
	float dotP, crossP[3];
	switch (type)
	{
	case 1:
	case 2:
		dotP = getDotP(g, v1);
		q[0] = dotP*v1[0];
		q[1] = dotP*v1[1];
		q[2] = dotP*v1[2];
		break;
	case 3:
	case 4:
		getCrossP(v3, g, crossP);
		getCrossP(crossP, v3, q);
		break;
	}
	if ( type == 2 || type == 4 )
	{
		q[0] = -q[0];
		q[1] = -q[1];
		q[2] = -q[2];
	}
}

void General_Data::facePhase(float* scalars,float* tensors,float* gradients, float *faceTable)
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Face Phase Begin:" << endl;
	change = 0.017f * gridSize;
	maxFaceIndex = gridx * (gridy-1) * (gridz-1) + (gridx-1) * gridy * (gridz-1) + (gridx-1) * (gridy-1) * gridz;
	faces = new Face[ maxFaceIndex ] ;
	totalFacePoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	ParallelFace parallel_face(gridx,gridy,gridz,totalFacePoints,edges,faces,dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,sizex,sizey,sizez,halfSize,scalars,tensors,gradients,faceTable);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_face,auto_partitioner());

	delete[] grid_2DGradMag;
	delete[] storeGrid_2DGradMag;
	timeEnd = clock();
	//t_facePhase_other = timeEnd-timeStart - t_facePhase_extremalPoint;
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	/*cout << "\tt_facePhase_other: " << t_facePhase_other << endl;
	cout << "\tt_facePhase_extremalPoint: " << t_facePhase_extremalPoint << endl;*/
	/*cout << "\tt_faceSampling_interp: " << t_faceSampling_interp/CLOCKS_PER_SEC << endl;
	cout << "\tt_faceOther_interp: " << t_faceOther_interp/CLOCKS_PER_SEC << endl;*/
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

int General_Data::getFaceIndex(int axis, int x, int y, int z) {
	int index;
	switch (axis)
	{
	case 1:
		index = x * (gridy - 1) * (gridz - 1) + y * (gridz - 1) + z;
		return index;
		break;
	case 2:
		index = gridx * (gridy - 1) * (gridz - 1) +
			x * gridy * (gridz - 1) + y * (gridz - 1) + z;
		return index;
		break;
	case 3:
		index = gridx * (gridy - 1) * (gridz - 1) +
			(gridx - 1) * gridy * (gridz - 1) +
			x * (gridy - 1) * gridz + y * gridz + z;
		return index;
		break;
	}
	return 0;
}

void General_Data::facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4,
	int index1, int index2, int index3, int index4, int axis, int i, int j, int k)
{
	Face *face = &(faces[getFaceIndex(axis, i, j, k)]);
	face->valid = true;
	float h[4];
	combineWindingNum(edge1, edge2, edge3, edge4, face, h, index1, index2, index3, index4, axis);
	getFaceIntersection(index1, index2, index3, index4, face, h);
	getCurveType(face);
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	get3DWindingNum(face, i, j, k, axis);
	/*timeEnd = clock();
	t_facePhase_extremalPoint = t_facePhase_extremalPoint + timeEnd - timeStart;*/
}

void General_Data::combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int index1, int index2, int index3, int index4, int axis)
{
	int f[4];
	float a[4], b[4], w[4];
	f[0] = edge1->f;
	f[1] = edge2->f;
	f[2] = edge3->f;
	f[3] = edge4->f;
	a[0] = edge1->a;
	a[1] = edge2->a;
	a[2] = -edge3->a;
	a[3] = -edge4->a;
	b[0] = edge1->b;
	b[1] = edge2->b;
	b[2] = -edge3->b;
	b[3] = -edge4->b;
	w[0] = edge1->w;
	w[1] = edge2->w;
	w[2] = -edge3->w;
	w[3] = -edge4->w;
	int d[4] = {0, 0, 1, 1};
	int F = f[0] + f[1] + f[2] + f[3];
	if ( F % 2 != 0 ) return;
	F = 0;
	float A = 0;
	float W = 0;
	h[0] = 0;
	for (int i = 0; i < 4; i ++)
	{
		if ( ( F + d[i] * f[i] ) % 2 != 0 )
		{
			a[i] = 2 * b[i] - a[i];
			w[i] = -w[i];
		}
		F += f[i];
		A += a[i];
		W += w[i];
		if ( i < 3 )
		{
			h[i+1] = W;
		}
	}
	A = smallAbsA(A);
	float wn = (W + A) / (2 * (float)M_PI);

	/*float vertices[4][3];
	float gradients[4][3];
	clock_t timeStart_i, timeEnd_i;
	for ( int i = 0; i < 3; i ++ )
	{
		getGridPointPosD(index1, vertices[0]);
		getGridPointPosD(index2, vertices[1]);
		getGridPointPosD(index3, vertices[2]);
		getGridPointPosD(index4, vertices[3]);
		timeStart_i = clock();
		getGradientCubic(vertices[0], gradients[0]);
		getGradientCubic(vertices[1], gradients[1]);
		getGradientCubic(vertices[2], gradients[2]);
		getGradientCubic(vertices[3], gradients[3]);
		timeEnd_i = clock();
		t_faceOther_interp = t_faceOther_interp + timeEnd_i - timeStart_i;
	}*/

	wn = floor( wn + 0.5f );
	if ( (int)wn % 2 != 0 )
	{
		face->extremal = true;
		for (int i = 1; i < 4; i ++)
		{
			h[i] += (float)i / 4 * A;
		}
	}
}

float General_Data::smallAbsA(float A)
{
	float A1 = fmod(A, 4*(float)M_PI);
	float A2;
	if ( A1 > 0 ) A2 = A1 - 4*(float)M_PI;
	else A2 = A1 + 4*(float)M_PI;
	if ( fabs(A1) > fabs(A2) ) return A2;
	else return A1;
}

void General_Data::rotate2D(float p[2], float theta)
{
	float ox = p[0];
	float oy = p[1];
	float x = ox*cos(theta)-oy*sin(theta);
	float y = ox*sin(theta)+oy*cos(theta);
	p[0] = x;
	p[1] = y;
}

void General_Data::getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h)
{
	if (!face->extremal) return;
	(*totalFacePoints++);
	float vertices[4][3];
	getGridPointPos(index1, vertices[0]);
	getGridPointPos(index2, vertices[1]);
	getGridPointPos(index3, vertices[2]);
	getGridPointPos(index4, vertices[3]);
	float projVertices[4][2];
	projVertices[0][0] = grid_2DGradMag[index1];
	projVertices[1][0] = grid_2DGradMag[index2];
	projVertices[2][0] = grid_2DGradMag[index3];
	projVertices[3][0] = grid_2DGradMag[index4];
	projVertices[0][1] = 0;
	projVertices[1][1] = 0;
	projVertices[2][1] = 0;
	projVertices[3][1] = 0;
	for (int i = 0; i < 4; i ++)
	{
		rotate2D(projVertices[i], h[i]);
	}
	int facePointNum = 4;
	float spoke[4], rim[4], alpha[4], omega[4], lambda[4];
	for (int i = 0; i < facePointNum; i ++)
	{
		spoke[i] = getNorm2(projVertices[i]);
	}
	float diff[2];
	for (int i = 0; i < facePointNum-1; i ++)
	{
		diff[0] = projVertices[i+1][0]-projVertices[i][0];
		diff[1] = projVertices[i+1][1]-projVertices[i][1];
		rim[i] = getNorm2(diff);
	}
	diff[0] = projVertices[0][0]-projVertices[facePointNum-1][0];
	diff[1] = projVertices[0][1]-projVertices[facePointNum-1][1];
	rim[facePointNum-1] = getNorm2(diff);
	for (int i = 0; i < facePointNum-1; i ++)
	{
		if (2*spoke[i]*spoke[i+1] == 0)
		{
			alpha[i] = 0;
		}
		else
		{
			float angleCos = (pow(spoke[i],2)+pow(spoke[i+1],2)-pow(rim[i],2)) / (2*spoke[i]*spoke[i+1]);
			if (angleCos > 1) angleCos = 1;
			if (angleCos < -1) angleCos = -1;
			alpha[i] = acos(angleCos);
		}
	}
	if (2*spoke[facePointNum-1]*spoke[0] == 0)
	{
		alpha[facePointNum-1] = 0;
	}
	else
	{
		float angleCos = (pow(spoke[facePointNum-1],2)+pow(spoke[0],2)-pow(rim[facePointNum-1],2))
			/ (2*spoke[facePointNum-1]*spoke[0]);
		if (angleCos > 1) angleCos = 1;
		if (angleCos < -1) angleCos = -1;
		alpha[facePointNum-1] = acos(angleCos);
	}
	for (int i = 1; i < facePointNum; i ++)
	{
		if (spoke[i] == 0) omega[i] = 0;
		else omega[i] = (tan(alpha[i-1]/2)+tan(alpha[i]/2))/spoke[i];
	}
	if (spoke[0] == 0) omega[0] = 0;
	else omega[0] = (tan(alpha[facePointNum-1]/2)+tan(alpha[0]/2))/spoke[0];
	float omegaSum = 0;
	for (int i = 0; i < facePointNum; i ++)
	{
		omegaSum += omega[i];
	}
	if (omegaSum == 0)
	{
		for (int i = 0; i < facePointNum; i ++)
		{
			lambda[i] = 0.25;
		}
	}
	else
	{
		for (int i = 0; i < facePointNum; i ++)
		{
			lambda[i] = omega[i] / omegaSum;
		}
	}
	face->facePoint[0] = 0;
	face->facePoint[1] = 0;
	face->facePoint[2] = 0;
	for (int i = 0; i < facePointNum; i ++)
	{
		face->facePoint[0] += lambda[i] * vertices[i][0];
		face->facePoint[1] += lambda[i] * vertices[i][1];
		face->facePoint[2] += lambda[i] * vertices[i][2];
	}
}

void General_Data::getCurveType(Face *face)
{
	if ( !face->extremal ) return;

	float *p = face->facePoint;
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	clock_t timeStart_i, timeEnd_i;
	timeStart_i = clock();
	getEigensolverCubic(p, &eigensolver);
	timeEnd_i = clock();
	//t_faceOther_interp = t_faceOther_interp + timeEnd_i - timeStart_i;

	float v1[3];
	getV1(&eigensolver, v1);
	float p10[3], p11[3], p13[3], p14[3];
	p10[0] = p[0]-2*change*v1[0];
	p10[1] = p[1]-2*change*v1[1];
	p10[2] = p[2]-2*change*v1[2];
	p11[0] = p[0]-change*v1[0];
	p11[1] = p[1]-change*v1[1];
	p11[2] = p[2]-change*v1[2];
	p13[0] = p[0]+change*v1[0];
	p13[1] = p[1]+change*v1[1];
	p13[2] = p[2]+change*v1[2];
	p14[0] = p[0]+2*change*v1[0];
	p14[1] = p[1]+2*change*v1[1];
	p14[2] = p[2]+2*change*v1[2];

	float v2[3];
	getV2(&eigensolver, v2);
	float p20[3], p21[3], p23[3], p24[3];
	p20[0] = p[0]-2*change*v2[0];
	p20[1] = p[1]-2*change*v2[1];
	p20[2] = p[2]-2*change*v2[2];
	p21[0] = p[0]-change*v2[0];
	p21[1] = p[1]-change*v2[1];
	p21[2] = p[2]-change*v2[2];
	p23[0] = p[0]+change*v2[0];
	p23[1] = p[1]+change*v2[1];
	p23[2] = p[2]+change*v2[2];
	p24[0] = p[0]+2*change*v2[0];
	p24[1] = p[1]+2*change*v2[1];
	p24[2] = p[2]+2*change*v2[2];

	float f, f10, f11, f13, f14, f20, f21, f23, f24;
	timeStart_i = clock();
	getScalarCubic(p, &f);
	getScalarCubic(p10, &f10);
	getScalarCubic(p11, &f11);
	getScalarCubic(p13, &f13);
	getScalarCubic(p14, &f14);
	getScalarCubic(p20, &f20);
	getScalarCubic(p21, &f21);
	getScalarCubic(p23, &f23);
	getScalarCubic(p24, &f24);
	timeEnd_i = clock();
	//t_faceOther_interp = t_faceOther_interp + timeEnd_i - timeStart_i;
	float secondDiff1 = -f10+16*f11-30*f+16*f13-f14;
	float secondDiff2 = -f20+16*f21-30*f+16*f23-f24;
	if (secondDiff1 * secondDiff2 < 0) face->faceTag = 3;
	else if (secondDiff1 > 0 && secondDiff2 > 0) face->faceTag = 2;
	else face->faceTag = 1;
}

void General_Data::get3DWindingNum(Face *face, int ii, int jj, int kk, int axis)
{
	/*float vertices[2][3];
	getGridPointPos(index1, vertices[0]);
	getGridPointPos(index3, vertices[1]);
	float subdVerts[subdNum][subdNum][3];
	float iStep, jStep;
	switch (axis)
	{
	case 1:
		iStep = (vertices[1][1] - vertices[0][1])/(subdNum-1);
		jStep = (vertices[1][2] - vertices[0][2])/(subdNum-1);
		for ( int i = 0; i < subdNum; i ++ )
		{
			for ( int j = 0; j < subdNum; j ++ )
			{
				subdVerts[i][j][0] = vertices[0][0];
				subdVerts[i][j][1] = vertices[0][1] + i*iStep;
				subdVerts[i][j][2] = vertices[0][2] + j*jStep;
			}
		}
		break;
	case 2:
		iStep = (vertices[1][2] - vertices[0][2])/(subdNum-1);
		jStep = (vertices[1][0] - vertices[0][0])/(subdNum-1);
		for ( int i = 0; i < subdNum; i ++ )
		{
			for ( int j = 0; j < subdNum; j ++ )
			{
				subdVerts[i][j][1] = vertices[0][1];
				subdVerts[i][j][2] = vertices[0][2] + i*iStep;
				subdVerts[i][j][0] = vertices[0][0] + j*jStep;
			}
		}
		break;
	case 3:
		iStep = (vertices[1][0] - vertices[0][0])/(subdNum-1);
		jStep = (vertices[1][1] - vertices[0][1])/(subdNum-1);
		for ( int i = 0; i < subdNum; i ++ )
		{
			for ( int j = 0; j < subdNum; j ++ )
			{
				subdVerts[i][j][2] = vertices[0][2];
				subdVerts[i][j][0] = vertices[0][0] + i*iStep;
				subdVerts[i][j][1] = vertices[0][1] + j*jStep;
			}
		}
		break;
	}*/
	float sGradients[subdNum][subdNum][3];
	clock_t timeStart_i, timeEnd_i;
	timeStart_i = clock();
	getGradientBicubicTable(axis, ii, jj, kk, sGradients);
	timeEnd_i = clock();
	//t_faceSampling_interp = t_faceSampling_interp + timeEnd_i - timeStart_i;
	float normGrads[subdNum][subdNum][3];
	for ( int i = 0; i < subdNum; i ++ )
	{
		for ( int j = 0; j < subdNum; j ++ )
		{
			//getShowPos(subdVerts[i][j], face->subdVerts[i][j]);
			//getGradientCubic(subdVerts[i][j], gradients[i][j]);
			float norm = getNorm(sGradients[i][j]);
			if (norm != 0)
			{
				normGrads[i][j][0] = sGradients[i][j][0]/norm;
				normGrads[i][j][1] = sGradients[i][j][1]/norm;
				normGrads[i][j][2] = sGradients[i][j][2]/norm;
			}
		}
	}
	float area = 0, currentArea;
	int count = 0;
	for ( int i = 0; i < subdNum - 1; i ++ )
	{
		for ( int j = 0; j < subdNum - 1; j ++ )
		{
			currentArea = signedArea3D(normGrads[i][j], normGrads[i][j+1], normGrads[i+1][j]);
			/*face->areas[count] = currentArea;
			face->centers[count][0] = (face->subdVerts[i][j][0]+face->subdVerts[i][j+1][0]+face->subdVerts[i+1][j][0])/3;
			face->centers[count][1] = (face->subdVerts[i][j][1]+face->subdVerts[i][j+1][1]+face->subdVerts[i+1][j][1])/3;
			face->centers[count][2] = (face->subdVerts[i][j][2]+face->subdVerts[i][j+1][2]+face->subdVerts[i+1][j][2])/3;
			count ++;*/
			area += currentArea;
			currentArea = signedArea3D(normGrads[i][j+1], normGrads[i+1][j+1], normGrads[i+1][j]);
			/*face->areas[count] = currentArea;
			face->centers[count][0] = (face->subdVerts[i][j+1][0]+face->subdVerts[i+1][j+1][0]+face->subdVerts[i+1][j][0])/3;
			face->centers[count][1] = (face->subdVerts[i][j+1][1]+face->subdVerts[i+1][j+1][1]+face->subdVerts[i+1][j][1])/3;
			face->centers[count][2] = (face->subdVerts[i][j+1][2]+face->subdVerts[i+1][j+1][2]+face->subdVerts[i+1][j][2])/3;
			count ++;*/
			area += currentArea;
		}
	}
	face->spherArea = area;
}

float General_Data::signedArea3D(float v1[3], float v2[3], float v3[3])
{
	float phiC = getDihedral(v1, v3, v2);
	float phiP1 = getDihedral(v2, v1, v3);
	float phiP2 = getDihedral(v3, v2, v1);
	float area = phiC+phiP1+phiP2-(float)M_PI;
	float first[3], second[3], up[3];
	first[0] = v1[0]-v3[0];
	first[1] = v1[1]-v3[1];
	first[2] = v1[2]-v3[2];
	second[0] = v2[0]-v1[0];
	second[1] = v2[1]-v1[1];
	second[2] = v2[2]-v1[2];
	getCrossP(first, second, up);
	int aSign = sign(getDotP(up, v3));
	area *= aSign;
	return area;
}

void General_Data::edgePhase(float* scalars,float* tensors,float* gradients, float *edgeTable)
{
	/*t_sampling = t_sampling_interp = t_sampling_other = t_extrEdge = t_extrEdge_orient = t_extrEdge_interp = t_extrEdge_other
		= t_orientV3 = t_propagate = t_propagate_getXY = t_propagate_getEdgeAB = t_projectG = t_getEdgeW
		= t_facePhase_other = t_facePhase_extremalPoint = t_cellPhase_other = t_cellPhase_extremalPoint = 0;*/
	//t_sampling_interp = t_extrEdge_interp = t_faceSampling_interp = t_faceOther_interp = t_cellPhase_interp = 0;

	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Edge Phase Begin:" << endl;
	change = 0.017f * gridSize;
	maxEdgeIndex = (gridx-1) * gridy * gridz + gridx * (gridy-1) * gridz + gridx * gridy * (gridz-1);
	edges = new Edge[ maxEdgeIndex ] ;
	grid_2DGradMag = new float[ maxGridIndex ] ;
	storeGrid_2DGradMag = new bool[ maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		storeGrid_2DGradMag[index] = false;
	}
	totalEdgePoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	ParallelEdge parallel_edge(gridx,gridy,gridz,totalEdgePoints,edges,allCubic,globalVec,dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,sizex,sizey,sizez,halfSize,scalars,tensors,gradients,edgeTable);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_edge,auto_partitioner());

	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	/*cout << "t_sampling: " << t_sampling << endl;
	cout << "\tt_sampling_interp: " << t_sampling_interp << endl;
	cout << "\tt_sampling_other: " << t_sampling_other << endl;
	cout << "t_extrEdge: " << t_extrEdge << endl;
	cout << "\tt_extrEdge_orient: " << t_extrEdge_orient << endl;
	cout << "\tt_extrEdge_interp: " << t_extrEdge_interp << endl;
	cout << "\tt_extrEdge_other: " << t_extrEdge_other << endl;
	cout << "t_orientV3: " << t_orientV3 << endl;
	cout << "t_propagate: " << t_propagate << endl;
	cout << "\tt_propagate_getXY: " << t_propagate_getXY << endl;
	cout << "\tt_propagate_getEdgeAB: " << t_propagate_getEdgeAB << endl;
	cout << "t_projectG: " << t_projectG << endl;
	cout << "t_getEdgeW: " << t_getEdgeW << endl;*/
	/*cout << "\tt_sampling_interp: " << t_sampling_interp/CLOCKS_PER_SEC << endl;
	cout << "\tt_extrEdge_interp: " << t_extrEdge_interp/CLOCKS_PER_SEC << endl;*/
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k)
{
	float *sampleG, *sampleV1, *sampleV3, *X, *Y, *sampleProjG;
	int edgeSampleNum;
	if (!adaptiveSamplingArrayTable(si1, si2, axis, i, j, k, &edgeSampleNum, &sampleG, &sampleV1, &sampleV3, &X, &Y, &sampleProjG)) return;
	//if (!adaptiveSamplingArray(p1, p2, &edgeSampleNum, &sampleG, &sampleV1, &sampleV3, &X, &Y, &sampleProjG)) return;
	Edge *edge = &(edges[getEdgeIndex(axis, i, j, k)]);
	edge->valid = true;
	extremalEdgeArray(index1, index2, edge, si1, si2, sampleV1, edgeSampleNum);
	orientV3Array(edge, sampleV3, edgeSampleNum);
	propagateXYArray(edge, sampleV3, X, Y, edgeSampleNum);
	projectGArray(sampleG, X, Y, sampleProjG, edgeSampleNum);
	getWindingNumArray(edge, sampleProjG, index1, index2, edgeSampleNum);
	delete[] sampleG;
	delete[] sampleV1;
	delete[] sampleV3;
	delete[] X;
	delete[] Y;
	delete[] sampleProjG;
}

void General_Data::extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	for ( int i = 1 ; i < edgeSampleNum ; i ++ )
	{
		float fSign = (float)sign(getDotP(sampleV1+3*i, sampleV1+3*(i-1)));
		if (fSign == 0) fSign = 1;
		sampleV1[3*i] = fSign*sampleV1[3*i];
		sampleV1[3*i+1] = fSign*sampleV1[3*i+1];
		sampleV1[3*i+2] = fSign*sampleV1[3*i+2];
	}
	/*timeEnd = clock();
	t_extrEdge_orient = t_extrEdge_orient + timeEnd - timeStart;*/
	float p1[3], p2[3];
	getGridPointPosD(index1, p1);
	getGridPointPosD(index2, p2);
	float g1[3], g2[3];
	getGradientGP(si1*3, g1);
	getGradientGP(si2*3, g2);
	float d1 = getDotP(sampleV1, g1);
	float d2 = getDotP(sampleV1+3*(edgeSampleNum-1), g2);
	if (sign(d1*d2) == -1)
	{
		(*totalEdgePoints)++;
		float p[3];
		linearInterpolate(d1, d2, 0, p1, p2, p);
		edge->edgePoint[0] = p[0];
		edge->edgePoint[1] = p[1];
		edge->edgePoint[2] = p[2];
		SelfAdjointEigenSolver<Matrix3f> eigensolver;
		clock_t timeStart_i, timeEnd_i;
		timeStart_i = clock();
		getEigensolverCubic(p, &eigensolver);
		timeEnd_i = clock();
		//t_extrEdge_interp = t_extrEdge_interp + timeEnd_i - timeStart_i;

		float v1[3];
		float p10[3], p11[3], p13[3], p14[3];
		getV1(&eigensolver, v1);
		p10[0] = p[0]-2*change*v1[0];
		p10[1] = p[1]-2*change*v1[1];
		p10[2] = p[2]-2*change*v1[2];
		p11[0] = p[0]-change*v1[0];
		p11[1] = p[1]-change*v1[1];
		p11[2] = p[2]-change*v1[2];
		p13[0] = p[0]+change*v1[0];
		p13[1] = p[1]+change*v1[1];
		p13[2] = p[2]+change*v1[2];
		p14[0] = p[0]+2*change*v1[0];
		p14[1] = p[1]+2*change*v1[1];
		p14[2] = p[2]+2*change*v1[2];

		float f10, f11, f13, f14, f;
		timeStart_i = clock();
		getScalarCubic(p10, &f10);
		getScalarCubic(p11, &f11);
		getScalarCubic(p13, &f13);
		getScalarCubic(p14, &f14);
		getScalarCubic(p, &f);
		timeEnd_i = clock();
		//t_extrEdge_interp = t_extrEdge_interp + timeEnd_i - timeStart_i;
		edge->extremal = true;
		float secondDiff = -f10+16*f11-30*f+16*f13-f14;
		if (secondDiff < 0) edge->edgeTag = 1;
		else edge->edgeTag = 2;
	}
	/*timeEnd = clock();
	t_extrEdge = t_extrEdge + timeEnd - timeStart;
	t_extrEdge_other = t_extrEdge - t_extrEdge_orient - t_extrEdge_interp;*/
}

void General_Data::orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	edge->f = 0;
	float oldEnd[3];
	oldEnd[0] = sampleV3[3*(edgeSampleNum-1)];
	oldEnd[1] = sampleV3[3*(edgeSampleNum-1)+1];
	oldEnd[2] = sampleV3[3*(edgeSampleNum-1)+2];
	for ( int i = 1 ; i < edgeSampleNum ; i ++ )
	{
		float fSign = (float)sign(getDotP(sampleV3+3*i, sampleV3+3*(i-1)));
		if (fSign == 0) fSign = 1;
		sampleV3[3*i] = fSign*sampleV3[3*i];
		sampleV3[3*i+1] = fSign*sampleV3[3*i+1];
		sampleV3[3*i+2] = fSign*sampleV3[3*i+2];
	}
	if (sign(getDotP(oldEnd, sampleV3+3*(edgeSampleNum-1))) == -1) edge->f = 1;
	/*timeEnd = clock();
	t_orientV3 = t_orientV3 + timeEnd - timeStart;*/
}

void General_Data::propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	float epsilon = 1e-12;
	float ox[3], oz[3];
	ox[0] = 1;
	ox[1] = 0;
	ox[2] = 0;
	oz[0] = 0;
	oz[1] = 0;
	oz[2] = 1;
	for ( int i = 0 ; i < edgeSampleNum ; i ++ )
	{
		float *v = sampleV3+3*i;
		float axle[3];
		getCrossP(oz, v, axle);
		float axleNorm = getNorm(axle);
		if (axleNorm > epsilon)
		{
			axle[0] = axle[0] / axleNorm;
			axle[1] = axle[1] / axleNorm;
			axle[2] = axle[2] / axleNorm;
			float angleCos = getDotP(oz, v);
			if (angleCos > 1) angleCos = 1;
			if (angleCos < -1) angleCos = -1;
			float angle = acos(angleCos);
			float a = getDotP(ox, axle)*(1-angleCos);
			float crossP[3];
			getCrossP(axle, ox, crossP);
			float b = sin(angle);
			ox[0] = ox[0]*angleCos + axle[0]*a + crossP[0]*b;
			ox[1] = ox[1]*angleCos + axle[1]*a + crossP[1]*b;
			ox[2] = ox[2]*angleCos + axle[2]*a + crossP[2]*b;
			oz[0] = v[0];
			oz[1] = v[1];
			oz[2] = v[2];
		}
		X[3*i] = ox[0];
		X[3*i+1] = ox[1];
		X[3*i+2] = ox[2];
		getCrossP(v, ox, Y+3*i);
	}
	/*timeEnd = clock();
	t_propagate_getXY = t_propagate_getXY + timeEnd - timeStart;*/
	float area = 0;
	float phiC = 0;
	for ( int i = 0 ; i < edgeSampleNum-1 ; i ++ )
	{
		float *p1 = sampleV3+3*i;
		float *p2 = sampleV3+3*(i+1);
		float incArea, incPhiC;
		signedSphericalTriangleArea(p1, p2, &incArea, &incPhiC);
		area += incArea;
		phiC += incPhiC;
	}
	edge->a = area;
	edge->b = phiC;
	/*timeEnd = clock();
	t_propagate = t_propagate + timeEnd - timeStart;
	t_propagate_getEdgeAB = t_propagate - t_propagate_getEdgeAB;*/
}

void General_Data::signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC)
{
	*phiC = getDihedral(p1, globalVec, p2);
	float phiP1 = getDihedral(p2, p1, globalVec);
	float phiP2 = getDihedral(globalVec, p2, p1);
	*area = *phiC+phiP1+phiP2-(float)M_PI;
	float first[3], second[3], up[3];
	first[0] = p1[0]-globalVec[0];
	first[1] = p1[1]-globalVec[1];
	first[2] = p1[2]-globalVec[2];
	second[0] = p2[0]-p1[0];
	second[1] = p2[1]-p1[1];
	second[2] = p2[2]-p1[2];
	getCrossP(first, second, up);
	int aSign = sign(getDotP(up, globalVec));
	*area *= aSign;
	*phiC *= aSign;
}

float General_Data::getDihedral(float b1[3], float b2[3], float b3[3])
{
	float b21[3], b23[3];
	getCrossP(b2, b1, b21);
	getCrossP(b2, b3, b23);
	float norm = getNorm(b21);
	if (norm != 0)
	{
		b21[0] = b21[0]/norm;
		b21[1] = b21[1]/norm;
		b21[2] = b21[2]/norm;
	}
	norm = getNorm(b23);
	if (norm != 0)
	{
		b23[0] = b23[0]/norm;
		b23[1] = b23[1]/norm;
		b23[2] = b23[2]/norm;
	}
	float dotP = getDotP(b21, b23);
	if (dotP > 1) dotP = 1;
	if (dotP < -1) dotP = -1;
	float phi = acos(dotP);
	return phi;
}

void General_Data::projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	for ( int i = 0 ; i < edgeSampleNum ; i ++ )
	{
		float *g = sampleG+3*i;
		float *x = X+3*i;
		float *y = Y+3*i;
		sampleProjG[3*i] = getDotP(g, x);
		sampleProjG[3*i+1] = getDotP(g, y);
		sampleProjG[3*i+2] = 0;
	}
	/*timeEnd = clock();
	t_projectG = t_projectG + timeEnd - timeStart;*/
}

void General_Data::getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	float totalAngle = 0;
	float epsilon = 1e-12;
	float *end1 = sampleProjG;
	for ( int i = 1 ; i < edgeSampleNum ; i ++ )
	{
		float *end2 = sampleProjG+3*i;
		/*float ang1, ang2;
		if (end1[0] != 0)
		{
			ang1 = atan(end1[1]/end1[0]);
		}
		else
		{
			ang1 = (float)M_PI/2;
		}
		if (end2[0] != 0)
		{
			ang2 = atan(end2[1]/end2[0]);
		}
		else
		{
			ang2 = (float)M_PI/2;
		}
		if (end1[0] < 0) ang1 = ang1-(float)M_PI;
		if (end2[0] < 0) ang2 = ang2+(float)M_PI;
		float localAngle = ang2-ang1;
		if (localAngle > M_PI) localAngle = localAngle-2*(float)M_PI;
		if (localAngle < -M_PI) localAngle = localAngle+2*(float)M_PI;*/
		float end1Norm = getNorm2(end1);
		float end2Norm = getNorm2(end2);
		if (end1Norm > 0 && end2Norm > 0)
		{
			float angleCos = getDotP(end1, end2) / (end1Norm*end2Norm);
			if (angleCos > 1) angleCos = 1;
			if (angleCos < -1) angleCos = -1;
			float localAngle = acos( angleCos );
			float crossP[3];
			getCrossP(end1, end2, crossP);
			localAngle *= sign(crossP[2]);
			totalAngle += localAngle;
		}
		end1 = end2;
	}
	edge->w = totalAngle;
	if (!storeGrid_2DGradMag[index1])
	{
		storeGrid_2DGradMag[index1] = true;
		grid_2DGradMag[index1] = getNorm(sampleProjG);
	}
	if (!storeGrid_2DGradMag[index2])
	{
		storeGrid_2DGradMag[index2] = true;
		grid_2DGradMag[index2] = getNorm(sampleProjG+3*(edgeSampleNum-1));
	}
	/*timeEnd = clock();
	t_getEdgeW = t_getEdgeW + timeEnd - timeStart;*/
}

bool General_Data::adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
	float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG)
{
	/*clock_t timeStart, timeEnd;
	timeStart = clock();*/
	float epsilon = 1e-12;
	int depth = 0;
	list<float> sampleGList, sampleV1List, sampleV3List;
	list<float> *sampleGLP, *sampleV1LP, *sampleV3LP;
	sampleGLP = &sampleGList;
	sampleV1LP = &sampleV1List;
	sampleV3LP = &sampleV3List;
	SelfAdjointEigenSolver<Matrix3f> eigensolver1, eigensolver2;
	float g1[3], g2[3];
	int gi1 = si1*3;
	int gi2 = si2*3;
	int ti1 = si1*6;
	int ti2 = si2*6;
	getGradientGP(gi1, g1);
	getGradientGP(gi2, g2);
	if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
	if ( !getEigensolverGP(ti1, &eigensolver1) ) return false;
	if ( !getEigensolverGP(ti2, &eigensolver2) ) return false;
	float v1_1[3], v3_1[3], v1_2[3], v3_2[3];
	getV1(&eigensolver1, v1_1);
	getV3(&eigensolver1, v3_1);
	getV1(&eigensolver2, v1_2);
	getV3(&eigensolver2, v3_2);

	for ( int j = 0; j < 3; j ++ )
	{
		sampleGLP->push_back(g1[j]);
		sampleV1LP->push_back(v1_1[j]);
		sampleV3LP->push_back(v3_1[j]);
	}
	if ( !adaptiveRecursionArrayTable(axis, ii, jj, kk, 2, g1, g2, v1_1, v1_2, v3_1, v3_2, sampleGLP, sampleV1LP, sampleV3LP, depth) ) return false;
	for ( int j = 0; j < 3; j ++ )
	{
		sampleGLP->push_back(g2[j]);
		sampleV1LP->push_back(v1_2[j]);
		sampleV3LP->push_back(v3_2[j]);
	}

	*edgeSampleNum = (int)(sampleGLP->size()/3);
	*sampleG = new float[3*(*edgeSampleNum)];
	*sampleV1 = new float[3*(*edgeSampleNum)];
	*sampleV3 = new float[3*(*edgeSampleNum)];
	*X = new float[3*(*edgeSampleNum)];
	*Y = new float[3*(*edgeSampleNum)];
	*sampleProjG = new float[3*(*edgeSampleNum)];
	list<float>::iterator ig = sampleGLP->begin();
	list<float>::iterator iv1 = sampleV1LP->begin();
	list<float>::iterator iv3 = sampleV3LP->begin();
	for ( int i = 0 ; ig != sampleGLP->end() ; ig ++, iv1 ++, iv3 ++, i ++ )
	{
		(*sampleG)[i] = *ig;
		(*sampleV1)[i] = *iv1;
		(*sampleV3)[i] = *iv3;
	}
	/*timeEnd = clock();
	t_sampling = t_sampling + timeEnd - timeStart;
	t_sampling_other = t_sampling - t_sampling_interp;*/
	return true;
}

//bool General_Data::adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
//	float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG)
//{
//	/*clock_t timeStart, timeEnd;
//	timeStart = clock();*/
//	float epsilon = 1e-6;
//	int depth = 0;
//	list<float> sampleGList, sampleV1List, sampleV3List;
//	list<float> *sampleGLP, *sampleV1LP, *sampleV3LP;
//	sampleGLP = &sampleGList;
//	sampleV1LP = &sampleV1List;
//	sampleV3LP = &sampleV3List;
//	
//	clock_t timeStart_i, timeEnd_i;
//	timeStart_i = clock();
//	SelfAdjointEigenSolver<Matrix3f> eigensolver1, eigensolver2;
//	float g1[3], g2[3];
//	if ( allCubic )
//	{
//		getGradientCubic(p1, g1);
//		getGradientCubic(p2, g2);
//		if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
//		if ( !getEigensolverCubic(p1, &eigensolver1) ) return false;
//		if ( !getEigensolverCubic(p2, &eigensolver2) ) return false;
//	}
//	else
//	{
//		getGradient(p1, g1);
//		getGradient(p2, g2);
//		if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
//		if ( !getEigensolver(p1, &eigensolver1) ) return false;
//		if ( !getEigensolver(p2, &eigensolver2) ) return false;
//	}
//	timeEnd_i = clock();
//	t_sampling_interp = t_sampling_interp + timeEnd_i - timeStart_i;
//
//	float v1_1[3], v3_1[3], v1_2[3], v3_2[3];
//	getV1(&eigensolver1, v1_1);
//	getV3(&eigensolver1, v3_1);
//	getV1(&eigensolver2, v1_2);
//	getV3(&eigensolver2, v3_2);
//
//	for ( int j = 0; j < 3; j ++ )
//	{
//		sampleGLP->push_back(g1[j]);
//		sampleV1LP->push_back(v1_1[j]);
//		sampleV3LP->push_back(v3_1[j]);
//	}
//	if ( !adaptiveRecursionArray(p1, p2, g1, g2, v1_1, v1_2, v3_1, v3_2, sampleGLP, sampleV1LP, sampleV3LP, depth) ) return false;
//	for ( int j = 0; j < 3; j ++ )
//	{
//		sampleGLP->push_back(g2[j]);
//		sampleV1LP->push_back(v1_2[j]);
//		sampleV3LP->push_back(v3_2[j]);
//	}
//
//	*edgeSampleNum = (int)(sampleGLP->size()/3);
//	*sampleG = new float[3*(*edgeSampleNum)];
//	*sampleV1 = new float[3*(*edgeSampleNum)];
//	*sampleV3 = new float[3*(*edgeSampleNum)];
//	*X = new float[3*(*edgeSampleNum)];
//	*Y = new float[3*(*edgeSampleNum)];
//	*sampleProjG = new float[3*(*edgeSampleNum)];
//	list<float>::iterator ig = sampleGLP->begin();
//	list<float>::iterator iv1 = sampleV1LP->begin();
//	list<float>::iterator iv3 = sampleV3LP->begin();
//	for ( int i = 0 ; ig != sampleGLP->end() ; ig ++, iv1 ++, iv3 ++, i ++ )
//	{
//		(*sampleG)[i] = *ig;
//		(*sampleV1)[i] = *iv1;
//		(*sampleV3)[i] = *iv3;
//	}
//	/*timeEnd = clock();
//	t_sampling = t_sampling + timeEnd - timeStart;
//	t_sampling_other = t_sampling - t_sampling_interp;*/
//	return true;
//}

bool General_Data::adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
	float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth)
{
	depth ++;
	if ( depth > 10 ) return true;
	float cosG, cosV1, cosV3;
	float normG = getNorm(g1)*getNorm(g2);
	float normV1 = getNorm(v1_1)*getNorm(v1_2);
	float normV3 = getNorm(v3_1)*getNorm(v3_2);
	if (normG == 0) cosG = 1;
	else cosG = getDotP(g1, g2)/normG;
	if (normV1 == 0) cosV1 = 1;
	else cosV1 = getDotP(v1_1, v1_2)/normV1;
	if (normV3 == 0) cosV3 = 1;
	else cosV3 = getDotP(v3_1, v3_2)/normV3;
	if ( cosG > thresh && fabs(cosV1) > thresh && fabs(cosV3) > thresh ) return true;

	/*float p[3];
	p[0] = (p1[0]+p2[0])/2;
	p[1] = (p1[1]+p2[1])/2;
	p[2] = (p1[2]+p2[2])/2;*/
	float g[3], v1[3], v3[3];

	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	clock_t timeStart, timeEnd;
	timeStart = clock();
	getGradientCubicTable(axis, ii, jj, kk, index, g);
	if ( !getEigensolverCubicTable(axis, ii, jj, kk, index, &eigensolver) ) return false;
	/*if ( allCubic )
	{
		getGradientCubic(p, g);
		if ( !getEigensolverCubic(p, &eigensolver) ) return false;
	}
	else
	{
		getGradient(p, g);
		if ( !getEigensolver(p, &eigensolver) ) return false;
	}*/
	timeEnd = clock();
	//t_sampling_interp = t_sampling_interp + timeEnd - timeStart;
	getV1(&eigensolver, v1);
	getV3(&eigensolver, v3);

	if ( !adaptiveRecursionArrayTable(axis, ii, jj, kk, index*2-1, g1, g, v1_1, v1, v3_1, v3, GL, V1L, V3L, depth) ) return false;
	for ( int j = 0; j < 3; j ++ )
	{
		GL->push_back(g[j]);
		V1L->push_back(v1[j]);
		V3L->push_back(v3[j]);
	}
	if ( !adaptiveRecursionArrayTable(axis, ii, jj, kk, index*2, g, g2, v1, v1_2, v3, v3_2, GL, V1L, V3L, depth) ) return false;
	return true;
}

//bool General_Data::adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
//	float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth)
//{
//	depth ++;
//	float cosG, cosV1, cosV3;
//	float normG = getNorm(g1)*getNorm(g2);
//	float normV1 = getNorm(v1_1)*getNorm(v1_2);
//	float normV3 = getNorm(v3_1)*getNorm(v3_2);
//	if (normG == 0) cosG = 1;
//	else cosG = getDotP(g1, g2)/normG;
//	if (normV1 == 0) cosV1 = 1;
//	else cosV1 = getDotP(v1_1, v1_2)/normV1;
//	if (normV3 == 0) cosV3 = 1;
//	else cosV3 = getDotP(v3_1, v3_2)/normV3;
//	if ( ( cosG > thresh && fabs(cosV1) > thresh && fabs(cosV3) > thresh ) || depth > 10 ) return true;
//
//	float p[3];
//	float g[3], v1[3], v3[3];
//	p[0] = (p1[0]+p2[0])/2;
//	p[1] = (p1[1]+p2[1])/2;
//	p[2] = (p1[2]+p2[2])/2;
//
//	clock_t timeStart, timeEnd;
//	timeStart = clock();
//	SelfAdjointEigenSolver<Matrix3f> eigensolver;
//	if ( allCubic )
//	{
//		getGradientCubic(p, g);
//		if ( !getEigensolverCubic(p, &eigensolver) ) return false;
//	}
//	else
//	{
//		getGradient(p, g);
//		if ( !getEigensolver(p, &eigensolver) ) return false;
//	}
//	timeEnd = clock();
//	t_sampling_interp = t_sampling_interp + timeEnd - timeStart;
//	getV1(&eigensolver, v1);
//	getV3(&eigensolver, v3);
//
//	if ( !adaptiveRecursionArray(p1, p, g1, g, v1_1, v1, v3_1, v3, GL, V1L, V3L, depth) ) return false;
//	for ( int j = 0; j < 3; j ++ )
//	{
//		GL->push_back(g[j]);
//		V1L->push_back(v1[j]);
//		V3L->push_back(v3[j]);
//	}
//	if ( !adaptiveRecursionArray(p, p2, g, g2, v1, v1_2, v3, v3_2, GL, V1L, V3L, depth) ) return false;
//	return true;
//}

int General_Data::getEdgeIndex(int axis, int x, int y, int z) {
	int index;
	switch (axis)
	{
	case 1:
		index = x * gridy * gridz + y * gridz + z;
		return index;
		break;
	case 2:
		index = (gridx - 1) * gridy * gridz +
			x * (gridy - 1) * gridz + y * gridz + z;
		return index;
		break;
	case 3:
		index = (gridx - 1) * gridy * gridz +
			gridx * (gridy - 1) * gridz +
			x * gridy * (gridz - 1) + y * (gridz - 1) + z;
		return index;
		break;
	}
	return 0;
}

void General_Data::moveCell(int axis, bool inc)
{
	switch ( axis )
	{
	case 1:
		if ( !inc && dispCellX > 0 )
		{
			dispCellX --;
			storeDispCell = false;
		}
		else if ( inc && dispCellX < gridx - 2 )
		{
			dispCellX ++;
			storeDispCell = false;
		}
		break;
	case 2:
		if ( !inc && dispCellY > 0 )
		{
			dispCellY --;
			storeDispCell = false;
		}
		else if ( inc && dispCellY < gridy - 2 )
		{
			dispCellY ++;
			storeDispCell = false;
		}
		break;
	case 3:
		if ( !inc && dispCellZ > 0 )
		{
			dispCellZ --;
			storeDispCell = false;
		}
		else if ( inc && dispCellZ < gridz - 2 )
		{
			dispCellZ ++;
			storeDispCell = false;
		}
		break;
	}
}

DispCell General_Data::getDispCell()
{
	if (storeDispCell) return dispCell;

	if (!storeGridShowPos)
	{
		computeGridShowPositions();
		storeGridShowPos = true;
	}
	if (!storeGridGradient)
	{
		computeGridGradients();
		storeGridGradient = true;
	}
	if (!storeGridVector)
	{
		computeGridVectors();
		storeGridVector = true;
	}
	dispCell.cell = cells[getIndex(1, 0, dispCellX, dispCellY, dispCellZ, gridx - 1, gridy - 1, gridz - 1)];
	dispCell.faces[0] = faces[getFaceIndex(1, dispCellX, dispCellY, dispCellZ)];
	dispCell.faces[1] = faces[getFaceIndex(1, dispCellX + 1, dispCellY, dispCellZ)];
	dispCell.faces[2] = faces[getFaceIndex(2, dispCellX, dispCellY, dispCellZ)];
	dispCell.faces[3] = faces[getFaceIndex(2, dispCellX, dispCellY + 1, dispCellZ)];
	dispCell.faces[4] = faces[getFaceIndex(3, dispCellX, dispCellY, dispCellZ)];
	dispCell.faces[5] = faces[getFaceIndex(3, dispCellX, dispCellY, dispCellZ + 1)];
	dispCell.edges[0] = edges[getEdgeIndex(1, dispCellX, dispCellY, dispCellZ)];
	dispCell.edges[1] = edges[getEdgeIndex(1, dispCellX, dispCellY + 1, dispCellZ)];
	dispCell.edges[2] = edges[getEdgeIndex(1, dispCellX, dispCellY + 1, dispCellZ + 1)];
	dispCell.edges[3] = edges[getEdgeIndex(1, dispCellX, dispCellY, dispCellZ + 1)];
	dispCell.edges[4] = edges[getEdgeIndex(2, dispCellX, dispCellY, dispCellZ)];
	dispCell.edges[5] = edges[getEdgeIndex(2, dispCellX + 1, dispCellY, dispCellZ)];
	dispCell.edges[6] = edges[getEdgeIndex(2, dispCellX + 1, dispCellY, dispCellZ + 1)];
	dispCell.edges[7] = edges[getEdgeIndex(2, dispCellX, dispCellY, dispCellZ + 1)];
	dispCell.edges[8] = edges[getEdgeIndex(3, dispCellX, dispCellY, dispCellZ)];
	dispCell.edges[9] = edges[getEdgeIndex(3, dispCellX + 1, dispCellY, dispCellZ)];
	dispCell.edges[10] = edges[getEdgeIndex(3, dispCellX + 1, dispCellY + 1, dispCellZ)];
	dispCell.edges[11] = edges[getEdgeIndex(3, dispCellX, dispCellY + 1, dispCellZ)];

	int indices[8];
	indices[0] = getIndex(1, 0, dispCellX, dispCellY, dispCellZ, gridx, gridy, gridz);
	indices[1] = getIndex(1, 0, dispCellX, dispCellY, dispCellZ + 1, gridx, gridy, gridz);
	indices[2] = getIndex(1, 0, dispCellX, dispCellY + 1, dispCellZ, gridx, gridy, gridz);
	indices[3] = getIndex(1, 0, dispCellX, dispCellY + 1, dispCellZ + 1, gridx, gridy, gridz);
	indices[4] = getIndex(1, 0, dispCellX + 1, dispCellY, dispCellZ, gridx, gridy, gridz);
	indices[5] = getIndex(1, 0, dispCellX + 1, dispCellY, dispCellZ + 1, gridx, gridy, gridz);
	indices[6] = getIndex(1, 0, dispCellX + 1, dispCellY + 1, dispCellZ, gridx, gridy, gridz);
	indices[7] = getIndex(1, 0, dispCellX + 1, dispCellY + 1, dispCellZ + 1, gridx, gridy, gridz);

	float cornerPos[8][3];
	dispCell.center[0] = 0;
	dispCell.center[1] = 0;
	dispCell.center[2] = 0;
	for ( int i = 0; i < 8; i ++ )
	{
		for ( int j = 0; j < 3; j ++ )
		{
			dispCell.corners[i][j] = gridShowPositions[indices[i]*3+j];
			dispCell.center[j] += dispCell.corners[i][j];
			dispCell.gradients[i][j] = gridGradients[indices[i]*3+j];
			dispCell.v1s[i][j] = gridV1s[indices[i]*3+j];
			dispCell.v3s[i][j] = gridV3s[indices[i]*3+j];
			cornerPos[i][j] = gridPoints[indices[i]*3+j];
		}
	}
	dispCell.center[0] /= 8;
	dispCell.center[1] /= 8;
	dispCell.center[2] /= 8;

	int edgesNum[24] = {5, 12, 8, 9, 6, 11, 7, 10, 1, 10, 4, 9,
		2, 11, 3, 12, 1, 6, 2, 5, 4, 7, 3, 8};
	for ( int i = 0; i < 24; i ++ )
	{
		edgesNum[i] -= 1;
	}
	for ( int i = 0; i < 6; i ++ )
	{
		if ( dispCell.faces[i].faceTag == -1 ) continue;
		dispCombineAB(&(dispCell.edges[edgesNum[4*i]]), &(dispCell.edges[edgesNum[4*i+1]]),
			&(dispCell.edges[edgesNum[4*i+2]]), &(dispCell.edges[edgesNum[4*i+3]]), &(dispCell.sumAB[i]));
		if (dispCell.faces[i].extremal)
		{
			SelfAdjointEigenSolver<Matrix3f> eigensolver;
			getEigensolverCubic(dispCell.faces[i].facePoint, &eigensolver);
			Vector3f ev = eigensolver.eigenvalues();
			dispCell.eigenvalues[i][0] = float(ev(2));
			dispCell.eigenvalues[i][1] = float(ev(1));
			dispCell.eigenvalues[i][2] = float(ev(0));
		}
	}

	int edgeIdx = 0;
	for ( int i = 0; i < 4; i ++ )
	{
		float p1[3], p2[3];
		for ( int j = 0; j < 3; j ++ )
		{
			p1[j] = cornerPos[2 * i][j];
			p2[j] = cornerPos[2 * i + 1][j];
		}
		dispCell.sampleP[edgeIdx].clear();
		dispCell.sampleG[edgeIdx].clear();
		dispCell.sampleV1[edgeIdx].clear();
		dispCell.sampleV3[edgeIdx].clear();
		dispCell.sampleX[edgeIdx].clear();
		dispCell.sampleY[edgeIdx].clear();
		dispCell.sampleProjG[edgeIdx].clear();
		if ( dispAdaptiveSampling(p1, p2, &(dispCell.sampleP[edgeIdx]), &(dispCell.sampleG[edgeIdx]),
			&(dispCell.sampleV1[edgeIdx]), &(dispCell.sampleV3[edgeIdx])) )
		{
			dispOrientV3(&(dispCell.sampleV3[edgeIdx]));
			dispPropagateXY(&(dispCell.sampleV3[edgeIdx]), &(dispCell.sampleX[edgeIdx]), &(dispCell.sampleY[edgeIdx]));
		}
		edgeIdx ++;
	}
	for ( int i = 0; i < 4; i ++ )
	{
		float p1[3], p2[3];
		for ( int j = 0; j < 3; j ++ )
		{
			p1[j] = cornerPos[i][j];
			p2[j] = cornerPos[i + 4][j];
		}
		dispCell.sampleP[edgeIdx].clear();
		dispCell.sampleG[edgeIdx].clear();
		dispCell.sampleV1[edgeIdx].clear();
		dispCell.sampleV3[edgeIdx].clear();
		dispCell.sampleX[edgeIdx].clear();
		dispCell.sampleY[edgeIdx].clear();
		dispCell.sampleProjG[edgeIdx].clear();
		if ( dispAdaptiveSampling(p1, p2, &(dispCell.sampleP[edgeIdx]), &(dispCell.sampleG[edgeIdx]),
			&(dispCell.sampleV1[edgeIdx]), &(dispCell.sampleV3[edgeIdx])) )
		{
			dispOrientV3(&(dispCell.sampleV3[edgeIdx]));
			dispPropagateXY(&(dispCell.sampleV3[edgeIdx]), &(dispCell.sampleX[edgeIdx]), &(dispCell.sampleY[edgeIdx]));
		}
		edgeIdx ++;
	}
	for ( int i = 0; i < 2; i ++ )
	{
		float p1[3], p2[3];
		for ( int j = 0; j < 3; j ++ )
		{
			p1[j] = cornerPos[i][j];
			p2[j] = cornerPos[i + 2][j];
		}
		dispCell.sampleP[edgeIdx].clear();
		dispCell.sampleG[edgeIdx].clear();
		dispCell.sampleV1[edgeIdx].clear();
		dispCell.sampleV3[edgeIdx].clear();
		dispCell.sampleX[edgeIdx].clear();
		dispCell.sampleY[edgeIdx].clear();
		dispCell.sampleProjG[edgeIdx].clear();
		if ( dispAdaptiveSampling(p1, p2, &(dispCell.sampleP[edgeIdx]), &(dispCell.sampleG[edgeIdx]),
			&(dispCell.sampleV1[edgeIdx]), &(dispCell.sampleV3[edgeIdx])) )
		{
			dispOrientV3(&(dispCell.sampleV3[edgeIdx]));
			dispPropagateXY(&(dispCell.sampleV3[edgeIdx]), &(dispCell.sampleX[edgeIdx]), &(dispCell.sampleY[edgeIdx]));
		}
		edgeIdx ++;
	}
	for ( int i = 4; i < 6; i ++ )
	{
		float p1[3], p2[3];
		for ( int j = 0; j < 3; j ++ )
		{
			p1[j] = cornerPos[i][j];
			p2[j] = cornerPos[i + 2][j];
		}
		dispCell.sampleP[edgeIdx].clear();
		dispCell.sampleG[edgeIdx].clear();
		dispCell.sampleV1[edgeIdx].clear();
		dispCell.sampleV3[edgeIdx].clear();
		dispCell.sampleX[edgeIdx].clear();
		dispCell.sampleY[edgeIdx].clear();
		dispCell.sampleProjG[edgeIdx].clear();
		if ( dispAdaptiveSampling(p1, p2, &(dispCell.sampleP[edgeIdx]), &(dispCell.sampleG[edgeIdx]),
			&(dispCell.sampleV1[edgeIdx]), &(dispCell.sampleV3[edgeIdx])) )
		{
			dispOrientV3(&(dispCell.sampleV3[edgeIdx]));
			dispPropagateXY(&(dispCell.sampleV3[edgeIdx]), &(dispCell.sampleX[edgeIdx]), &(dispCell.sampleY[edgeIdx]));
		}
		edgeIdx ++;
	}

	Cell *cell = &(dispCell.cell);
	if ( cell->extrP )
	{
		float change = 0.017f * gridSize;
		float *p = cell->extremalPoint;
		float p10[3], p11[3], p13[3], p14[3];
		p10[0] = p[0]-2*change;
		p10[1] = p[1];
		p10[2] = p[2];
		p11[0] = p[0]-change;
		p11[1] = p[1];
		p11[2] = p[2];
		p13[0] = p[0]+change;
		p13[1] = p[1];
		p13[2] = p[2];
		p14[0] = p[0]+2*change;
		p14[1] = p[1];
		p14[2] = p[2];

		float p20[3], p21[3], p23[3], p24[3];
		p20[0] = p[0];
		p20[1] = p[1]-2*change;
		p20[2] = p[2];
		p21[0] = p[0];
		p21[1] = p[1]-change;
		p21[2] = p[2];
		p23[0] = p[0];
		p23[1] = p[1]+change;
		p23[2] = p[2];
		p24[0] = p[0];
		p24[1] = p[1]+2*change;
		p24[2] = p[2];

		float p30[3], p31[3], p33[3], p34[3];
		p30[0] = p[0];
		p30[1] = p[1];
		p30[2] = p[2]-2*change;
		p31[0] = p[0];
		p31[1] = p[1];
		p31[2] = p[2]-change;
		p33[0] = p[0];
		p33[1] = p[1];
		p33[2] = p[2]+change;
		p34[0] = p[0];
		p34[1] = p[1];
		p34[2] = p[2]+2*change;

		float f, f10, f11, f13, f14, f20, f21, f23, f24, f30, f31, f33, f34;
		getScalarCubic(p, &f);
		getScalarCubic(p10, &f10);
		getScalarCubic(p11, &f11);
		getScalarCubic(p13, &f13);
		getScalarCubic(p14, &f14);
		getScalarCubic(p20, &f20);
		getScalarCubic(p21, &f21);
		getScalarCubic(p23, &f23);
		getScalarCubic(p24, &f24);
		getScalarCubic(p30, &f30);
		getScalarCubic(p31, &f31);
		getScalarCubic(p33, &f33);
		getScalarCubic(p34, &f34);
		float secondDiff1 = -f10+16*f11-30*f+16*f13-f14;
		float secondDiff2 = -f20+16*f21-30*f+16*f23-f24;
		float secondDiff3 = -f30+16*f31-30*f+16*f33-f34;
		float scale = 0.0001;
		cout << "extremal point type: " << endl;
		cout << "secondDiff1:" << setw(7) << floor(secondDiff1 / scale + 0.5) * scale << endl;
		cout << setw(7) << floor(f10 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f11 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f13 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f14 / scale + 0.5) * scale << endl;
		cout << "secondDiff2:" << setw(7) << floor(secondDiff2 / scale + 0.5) * scale << endl;
		cout << setw(7) << floor(f20 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f21 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f23 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f24 / scale + 0.5) * scale << endl;
		cout << "secondDiff3:" << setw(7) << floor(secondDiff3 / scale + 0.5) * scale << endl;
		cout << setw(7) << floor(f30 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f31 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f33 / scale + 0.5) * scale << ", " <<
			setw(7) << floor(f34 / scale + 0.5) * scale << endl;
	}

	storeDispCell = true;
	return dispCell;
}

int General_Data::getCellX()
{
	return dispCellX + 1;
}

int General_Data::getCellY()
{
	return dispCellY + 1;
}

int General_Data::getCellZ()
{
	return dispCellZ + 1;
}

void General_Data::setCellX(int value)
{
	dispCellX = value - 1;
	storeDispCell = false;
}

void General_Data::setCellY(int value)
{
	dispCellY = value - 1;
	storeDispCell = false;
}

void General_Data::setCellZ(int value)
{
	dispCellZ = value - 1;
	storeDispCell = false;
}

bool General_Data::dispAdaptiveSampling(float p1[3], float p2[3], vector<float> *sampleP,
		vector<float> *sampleG, vector<float> *sampleV1, vector<float> *sampleV3)
{
	float epsilon = 1e-12;
	int depth = 0;
	float g1[3], g2[3], v1_1[3], v3_1[3], v1_2[3], v3_2[3], showPos1[3], showPos2[3];
	SelfAdjointEigenSolver<Matrix3f> eigensolver1, eigensolver2;
	if ( allCubic )
	{
		getGradientCubic(p1, g1);
		getGradientCubic(p2, g2);
		if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
		getEigensolverCubic(p1, &eigensolver1);
		getEigensolverCubic(p2, &eigensolver2);
	}
	else
	{
		getGradient(p1, g1);
		getGradient(p2, g2);
		if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
		getEigensolver(p1, &eigensolver1);
		getEigensolver(p2, &eigensolver2);
	}
	getV1(&eigensolver1, v1_1);
	getV3(&eigensolver1, v3_1);
	getV1(&eigensolver2, v1_2);
	getV3(&eigensolver2, v3_2);
	getShowPos(p1, showPos1);
	sampleP->push_back(showPos1[0]);
	sampleP->push_back(showPos1[1]);
	sampleP->push_back(showPos1[2]);
	sampleV3->push_back(v3_1[0]);
	sampleV3->push_back(v3_1[1]);
	sampleV3->push_back(v3_1[2]);
	dispAdaptiveRecursion(p1, p2, g1, g2, v1_1, v1_2, v3_1, v3_2, sampleP, sampleG, sampleV1, sampleV3, depth);
	getShowPos(p2, showPos2);
	sampleP->push_back(showPos2[0]);
	sampleP->push_back(showPos2[1]);
	sampleP->push_back(showPos2[2]);
	sampleV3->push_back(v3_2[0]);
	sampleV3->push_back(v3_2[1]);
	sampleV3->push_back(v3_2[2]);
	return true;
}

void General_Data::dispAdaptiveRecursion(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3], float v3_1[3],
	float v3_2[3], vector<float> *PL, vector<float> *GL, vector<float> *V1L, vector<float> *V3L, int depth)
{
	depth ++;
	float cosG, cosV1, cosV3;
	float normG = getNorm(g1)*getNorm(g2);
	float normV1 = getNorm(v1_1)*getNorm(v1_2);
	float normV3 = getNorm(v3_1)*getNorm(v3_2);
	if (normG == 0) cosG = 1;
	else cosG = getDotP(g1, g2)/normG;
	if (normV1 == 0) cosV1 = 1;
	else cosV1 = getDotP(v1_1, v1_2)/normV1;
	if (normV3 == 0) cosV3 = 1;
	else cosV3 = getDotP(v3_1, v3_2)/normV3;
	if ( ( cosG > thresh && fabs(cosV1) > thresh && fabs(cosV3) > thresh ) || depth > 10 ) return;

	float p[3];
	p[0] = (p1[0]+p2[0])/2;
	p[1] = (p1[1]+p2[1])/2;
	p[2] = (p1[2]+p2[2])/2;
	vector<float> PL1, PL2, GL1, GL2, V1L1, V1L2, V3L1, V3L2;

	float g[3], v1[3], v3[3];
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	if ( allCubic )
	{
		getGradientCubic(p, g);
		getEigensolverCubic(p, &eigensolver);
	}
	else
	{
		getGradient(p, g);
		getEigensolver(p, &eigensolver);
	}
	getV1(&eigensolver, v1);
	getV3(&eigensolver, v3);

	dispAdaptiveRecursion(p1, p, g1, g, v1_1, v1, v3_1, v3, &PL1, &GL1, &V1L1, &V3L1, depth);
	dispAdaptiveRecursion(p, p2, g, g2, v1, v1_2, v3, v3_2, &PL2, &GL2, &V1L2, &V3L2, depth);

	PL->insert(PL->end(),PL1.begin(),PL1.end());
	float showPos[3];
	getShowPos(p, showPos);
	PL->push_back(showPos[0]);
	PL->push_back(showPos[1]);
	PL->push_back(showPos[2]);
	PL->insert(PL->end(),PL2.begin(),PL2.end());
	GL->insert(GL->end(),GL1.begin(),GL1.end());
	GL->push_back(g[0]);
	GL->push_back(g[1]);
	GL->push_back(g[2]);
	GL->insert(GL->end(),GL2.begin(),GL2.end());
	V1L->insert(V1L->end(),V1L1.begin(),V1L1.end());
	V1L->push_back(v1[0]);
	V1L->push_back(v1[1]);
	V1L->push_back(v1[2]);
	V1L->insert(V1L->end(),V1L2.begin(),V1L2.end());
	V3L->insert(V3L->end(),V3L1.begin(),V3L1.end());
	V3L->push_back(v3[0]);
	V3L->push_back(v3[1]);
	V3L->push_back(v3[2]);
	V3L->insert(V3L->end(),V3L2.begin(),V3L2.end());
}

void General_Data::dispOrientV3(vector<float> *sampleV3)
{
	int edgeSampleNum = (int)sampleV3->size()/3;
	float oldEnd[3], oldV[3], newV[3];
	oldEnd[0] = (*sampleV3)[3*(edgeSampleNum-1)];
	oldEnd[1] = (*sampleV3)[3*(edgeSampleNum-1)+1];
	oldEnd[2] = (*sampleV3)[3*(edgeSampleNum-1)+2];
	oldV[0] = (*sampleV3)[0];
	oldV[1] = (*sampleV3)[1];
	oldV[2] = (*sampleV3)[2];
	for ( int i = 1 ; i < edgeSampleNum ; i ++ )
	{
		newV[0] = (*sampleV3)[3*i];
		newV[1] = (*sampleV3)[3*i+1];
		newV[2] = (*sampleV3)[3*i+2];
		float fSign = (float)sign(getDotP(newV, oldV));
		if (fSign == 0) fSign = 1;
		(*sampleV3)[3*i] = fSign*(*sampleV3)[3*i];
		(*sampleV3)[3*i+1] = fSign*(*sampleV3)[3*i+1];
		(*sampleV3)[3*i+2] = fSign*(*sampleV3)[3*i+2];
		oldV[0] = fSign*newV[0];
		oldV[1] = fSign*newV[1];
		oldV[2] = fSign*newV[2];
	}
}

void General_Data::dispPropagateXY(vector<float> *sampleV3, vector<float> *sampleX, vector<float> *sampleY)
{
	int edgeSampleNum = (int)sampleV3->size()/3;
	float epsilon = 1e-12;
	float ox[3], oz[3], v[3], y[3];
	ox[0] = 1;
	ox[1] = 0;
	ox[2] = 0;
	oz[0] = 0;
	oz[1] = 0;
	oz[2] = 1;
	for ( int i = 0 ; i < edgeSampleNum ; i ++ )
	{
		v[0] = (*sampleV3)[3*i];
		v[1] = (*sampleV3)[3*i+1];
		v[2] = (*sampleV3)[3*i+2];
		float axle[3];
		getCrossP(oz, v, axle);
		float axleNorm = getNorm(axle);
		if (axleNorm > epsilon)
		{
			axle[0] = axle[0] / axleNorm;
			axle[1] = axle[1] / axleNorm;
			axle[2] = axle[2] / axleNorm;
			float angleCos = getDotP(oz, v);
			if (angleCos > 1) angleCos = 1;
			if (angleCos < -1) angleCos = -1;
			float angle = acos(angleCos);
			float a = getDotP(ox, axle)*(1-angleCos);
			float crossP[3];
			getCrossP(axle, ox, crossP);
			float b = sin(angle);
			ox[0] = ox[0]*angleCos + axle[0]*a + crossP[0]*b;
			ox[1] = ox[1]*angleCos + axle[1]*a + crossP[1]*b;
			ox[2] = ox[2]*angleCos + axle[2]*a + crossP[2]*b;
			oz[0] = v[0];
			oz[1] = v[1];
			oz[2] = v[2];
		}
		sampleX->push_back(ox[0]);
		sampleX->push_back(ox[1]);
		sampleX->push_back(ox[2]);
		getCrossP(v, ox, y);
		sampleY->push_back(y[0]);
		sampleY->push_back(y[1]);
		sampleY->push_back(y[2]);
	}
}

void General_Data::dispCombineAB(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, float *sumAB)
{
	int f[4];
	float a[4], b[4], w[4];
	f[0] = edge1->f;
	f[1] = edge2->f;
	f[2] = edge3->f;
	f[3] = edge4->f;
	a[0] = edge1->a;
	a[1] = edge2->a;
	a[2] = -edge3->a;
	a[3] = -edge4->a;
	b[0] = edge1->b;
	b[1] = edge2->b;
	b[2] = -edge3->b;
	b[3] = -edge4->b;
	w[0] = edge1->w;
	w[1] = edge2->w;
	w[2] = -edge3->w;
	w[3] = -edge4->w;
	int d[4] = {0, 0, 1, 1};
	int F = f[0] + f[1] + f[2] + f[3];
	if ( F % 2 != 0 ) return;
	F = 0;
	float A = 0;
	float W = 0;
	for (int i = 0; i < 4; i ++)
	{
		if ( ( F + d[i] * f[i] ) % 2 != 0 )
		{
			a[i] = 2 * b[i] - a[i];
			w[i] = -w[i];
		}
		F += f[i];
		A += a[i];
		W += w[i];
	}
	A = smallAbsA(A);
	*sumAB = A + W;
}

void General_Data::saveDisplay(bool resultDebug)
{
	if (resultDebug)
	{
		showGridScalar = view3D->showGridScalar;
		showGridV1 = view3D->showGridV1;
		showGridV2 = view3D->showGridV2;
		showGridV3 = view3D->showGridV3;
		showGridGradient = view3D->showGridGradient;
		showEdgePoint = view3D->showEdgePoint;
		showFacePoint = view3D->showFacePoint;
		showCellPoint = view3D->showCellPoint;
		showCell = view3D->showCell;
		showSample = view3D->showSample;
		showFaceSampleGradient = view3D->showFaceSampleGradient;
	}
	else
	{
		showGridScalar = showGridV1 = showGridV2 = showGridV3 = showGridGradient = 
			showEdgePoint = showFacePoint = showCellPoint = showCell = showSample = showFaceSampleGradient = false;
	}
	pointRatio = view3D->pointRatio;
	curveRatio = view3D->curveRatio;
	surfaceRatio = view3D->surfaceRatio;
	localIntensityThreshMinG = view3D->localIntensityThreshMinG;
	localIntensityThreshMaxG = view3D->localIntensityThreshMaxG;
	eigenvalueThresh = view3D->eigenvalueThresh;
	isovalue = view3D->isovalue;
	saddleCurve = view3D->saddleCurve;
	saddlePoint = view3D->saddlePoint;
	hideCurve = view3D->hideCurve;
	hideSurface = view3D->hideSurface;
	hidePoint = view3D->hidePoint;
	minCurve = view3D->minCurve;
	minSurface = view3D->minSurface;
	minPoint = view3D->minPoint;
	maxCurve = view3D->maxCurve;
	maxSurface = view3D->maxSurface;
	maxPoint = view3D->maxPoint;
	showIsosurface = view3D->showIsosurface;
	isosurfaceFace = view3D->isosurfaceFace;
	cylinderCurve = view3D->cylinderCurve;
	shadingFace = view3D->shadingFace;
	shadingWireframe = view3D->shadingWireframe;
	hideSaliencyRatio = view3D->hideSaliencyRatio;
	hideLocalIntensity = view3D->hideLocalIntensity;
	hideEigenvalue = view3D->hideEigenvalue;
}

void General_Data::loadDisplay(bool resultDebug)
{
	if (resultDebug)
	{
		view3D->showGridScalar = showGridScalar;
		view3D->showGridV1 = showGridV1;
		view3D->showGridV2 = showGridV2;
		view3D->showGridV3 = showGridV3;
		view3D->showGridGradient = showGridGradient;
		view3D->showEdgePoint = showEdgePoint;
		view3D->showFacePoint = showFacePoint;
		view3D->showCellPoint = showCellPoint;
		view3D->showCell = showCell;
		view3D->showSample = showSample;
		view3D->showFaceSampleGradient = showFaceSampleGradient;
	}
	else
	{
		view3D->showGridScalar = view3D->showGridV1 = view3D->showGridV2 = view3D->showGridV3 = 
			view3D->showGridGradient = view3D->showEdgePoint = view3D->showFacePoint = 
			view3D->showCellPoint = view3D->showCell = view3D->showSample = view3D->showFaceSampleGradient = false;
	}
	view3D->pointRatio = pointRatio;
	view3D->curveRatio = curveRatio;
	view3D->surfaceRatio = surfaceRatio;
	view3D->localIntensityThreshMinG = localIntensityThreshMinG;
	view3D->localIntensityThreshMaxG = localIntensityThreshMaxG;
	view3D->eigenvalueThresh = eigenvalueThresh;
	view3D->isovalue = isovalue;
	view3D->saddleCurve = saddleCurve;
	view3D->saddlePoint = saddlePoint;
	view3D->hideCurve = hideCurve;
	view3D->hideSurface = hideSurface;
	view3D->hidePoint = hidePoint;
	view3D->minCurve = minCurve;
	view3D->minSurface = minSurface;
	view3D->minPoint = minPoint;
	view3D->maxCurve = maxCurve;
	view3D->maxSurface = maxSurface;
	view3D->maxPoint = maxPoint;
	view3D->showIsosurface = showIsosurface;
	view3D->isosurfaceFace = isosurfaceFace;
	view3D->cylinderCurve = cylinderCurve;
	view3D->shadingFace = shadingFace;
	view3D->shadingWireframe = shadingWireframe;
	view3D->hideSaliencyRatio = hideSaliencyRatio;
	view3D->hideLocalIntensity = hideLocalIntensity;
	view3D->hideEigenvalue = hideEigenvalue;
}

void General_Data::generateCubicVolume()
{
	switch (dataType)
	{
	case 1:
		generateCubicMRC();
		break;
	case 2:
		generateCubicDTI();
		break;
	}
}

void General_Data::generateCubicMRC()
{
	FILE* fin = fopen( dataName, "rb" ) ;
	FILE* oldFin = fin ;

	// Parse header
	int osizex, osizey, osizez, mode, nsize, lmarginx, lmarginy, lmarginz, rmarginx, rmarginy, rmarginz;
	fread( &osizex, sizeof( int ), 1, fin ) ;
	fread( &osizey, sizeof( int ), 1, fin ) ;
	fread( &osizez, sizeof( int ), 1, fin ) ;
	fread( &mode, sizeof( int ), 1, fin ) ;
	nsize = osizex;
	if ( nsize < osizey ) nsize = osizey;
	if ( nsize < osizez ) nsize = osizez;
	lmarginx = floor( (float)( nsize - osizex ) / 2 );
	lmarginy = floor( (float)( nsize - osizey ) / 2 );
	lmarginz = floor( (float)( nsize - osizez ) / 2 );
	rmarginx = osizex + lmarginx - 1;
	rmarginy = osizey + lmarginy - 1;
	rmarginz = osizez + lmarginz - 1;

	// Read volume
	fin = oldFin;
	fseek( fin, 1024, SEEK_SET ) ;
	char chard ;
	short shortd ;
	float floatd ;
	float d ;
	float* data = new float[ nsize * nsize * nsize ];
	for ( int i = 0 ; i < nsize ; i ++ )
		for ( int j = 0 ; j < nsize ; j ++ )
			for ( int k = 0 ; k < nsize ; k ++ )
			{
				if ( i<lmarginz || i>rmarginz || j<lmarginy || j>rmarginy || k<lmarginx || k>rmarginx )
				{
					if ( reverse ) data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = -maxScalar;
					else data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = minScalar;
					continue;
				}
				switch ( mode )
				{
				case 0:
					fread( &chard, sizeof( char ), 1, fin ) ;
					d = (float) chard ;
					break ;
				case 1:
					fread( &shortd, sizeof( short ), 1, fin ) ;
					d = (float) shortd ;
					break ;
				case 2:
					fread( &floatd, sizeof( float ), 1, fin ) ;
					d = floatd ;
					break ;
				}
				data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = d;
			}
	fclose( fin ) ;

	FILE* out = fopen( "../../Data/out.mrc", "wb" );

	fwrite( &nsize, sizeof( int ), 1, out );
	fwrite( &nsize, sizeof( int ), 1, out );
	fwrite( &nsize, sizeof( int ), 1, out );

	mode = 2 ;
	fwrite( &mode, sizeof ( int ), 1, out );

	int start[3] = {0,0,0};
	int interval[3] = { nsize - 1, nsize - 1, nsize - 1 };
	fwrite( start, sizeof( int ), 3, out );
	fwrite( interval, sizeof( int ), 3, out );

	float cella[3] = {1,1,1};
	float cellb[3] = {90,90,90};
	fwrite( cella, sizeof( float ), 3, out );
	fwrite( cellb, sizeof( float ), 3, out );

	int map[3] = {1,2,3} ;
	fwrite( map, sizeof( int ), 3, out );

	float Imin = 10000, Imax = -10000, Imean = 0;
	/*for (int i = 0 ; i < nsize * nsize * nsize ; i++ )
	{
		Imean += data[i];
		if ( data[i] < Imin )
		{
			Imin = data[i];
		}
		if ( data[i] > Imax )
		{
			Imax = data[i];
		}
	}
	Imean /= nsize * nsize * nsize;*/
	float D[3] = {Imin, Imax, Imean};
	fwrite( D, sizeof( float ), 3, out);

	int zero = 0;
	for (int i = 22 ; i < 256 ; i++ )
	{
		fwrite( &zero, sizeof( int ), 1, out);
	}

	for ( int i = 0 ; i < nsize ; i ++ )
		for ( int j = 0 ; j < nsize ; j ++ )
			for ( int k = 0 ; k < nsize ; k ++ )
			{
				float I = data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)];
				fwrite( &I, sizeof( float ), 1, out );
			}

	fclose( out );
	delete[] data;
}

void General_Data::generateCubicDTI()
{
	FILE* fin = fopen( dataName, "rb" ) ;

	int osizex, osizey, osizez, mode, nsize, lmarginx, lmarginy, lmarginz, rmarginx, rmarginy, rmarginz;
	osizex = width;
	osizey = height;
	osizez = slices;
	nsize = osizex;
	if ( nsize < osizey ) nsize = osizey;
	if ( nsize < osizez ) nsize = osizez;
	lmarginx = floor( (float)( nsize - osizex ) / 2 );
	lmarginy = floor( (float)( nsize - osizey ) / 2 );
	lmarginz = floor( (float)( nsize - osizez ) / 2 );
	rmarginx = osizex + lmarginx - 1;
	rmarginy = osizey + lmarginy - 1;
	rmarginz = osizez + lmarginz - 1;

	// Read FA
	float d ;
	float* data = new float[ nsize * nsize * nsize ];
	for ( int i = 0 ; i < nsize ; i ++ )
		for ( int j = 0 ; j < nsize ; j ++ )
			for ( int k = 0 ; k < nsize ; k ++ )
			{
				if ( i<lmarginz || i>rmarginz || j<lmarginy || j>rmarginy || k<lmarginx || k>rmarginx )
				{
					if ( reverse ) data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = -maxScalar;
					else data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = minScalar;
					continue;
				}
				fread( &d, sizeof( float ), 1, fin ) ;
				data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)] = d;
			}
	fclose( fin ) ;

	FILE* out = fopen( "../../Data/out.mrc", "wb" );

	fwrite( &nsize, sizeof( int ), 1, out );
	fwrite( &nsize, sizeof( int ), 1, out );
	fwrite( &nsize, sizeof( int ), 1, out );

	mode = 2 ;
	fwrite( &mode, sizeof ( int ), 1, out );

	int start[3] = {0,0,0};
	int interval[3] = { nsize - 1, nsize - 1, nsize - 1 };
	fwrite( start, sizeof( int ), 3, out );
	fwrite( interval, sizeof( int ), 3, out );

	float cella[3] = {1,1,1};
	float cellb[3] = {90,90,90};
	fwrite( cella, sizeof( float ), 3, out );
	fwrite( cellb, sizeof( float ), 3, out );

	int map[3] = {1,2,3} ;
	fwrite( map, sizeof( int ), 3, out );

	float Imin = 10000, Imax = -10000, Imean = 0;
	/*for (int i = 0 ; i < nsize * nsize * nsize ; i++ )
	{
		Imean += data[i];
		if ( data[i] < Imin )
		{
			Imin = data[i];
		}
		if ( data[i] > Imax )
		{
			Imax = data[i];
		}
	}
	Imean /= nsize * nsize * nsize;*/
	float D[3] = {Imin, Imax, Imean};
	fwrite( D, sizeof( float ), 3, out);

	int zero = 0;
	for (int i = 22 ; i < 256 ; i++ )
	{
		fwrite( &zero, sizeof( int ), 1, out);
	}

	for ( int i = 0 ; i < nsize ; i ++ )
		for ( int j = 0 ; j < nsize ; j ++ )
			for ( int k = 0 ; k < nsize ; k ++ )
			{
				float I = data[getIndex(1, 0, k, j, i, nsize, nsize, nsize)];
				fwrite( &I, sizeof( float ), 1, out );
			}

	fclose( out );
	delete[] data;
}

void General_Data::saveOFF(ofstream &ofs)
{
	ofs << "OFF" << endl;
	ofs << "# " << endl;
	int vertNum = vertices->size();
	//int geoNum = quads->size() + segments->size() + points->size();
	ofs << vertNum << " " << quads->size() << " 0" << endl;
	ofs << "# " << segments->size() << " " << points->size() << endl;
	ofs << "# " << maxIntensity << " " << minIntensity << " " << *maxEigenvalue << " " << *minEigenvalue << endl;
	for (int i = 0; i < vertices->size(); i++)
	{
		Vertex vertex = (*vertices)[i];
		ofs << vertex.position[0] << " " << vertex.position[1] << " " << vertex.position[2] << endl;
	}
	for (int i = 0; i < quads->size(); i++)
	{
		Quad quad = (*quads)[i];
		ofs << 4 << " " << quad.vertIdxs[0] << " " << quad.vertIdxs[1] << " " << quad.vertIdxs[2] << " " << quad.vertIdxs[3] << " "
			<< quad.relativeSaliencies[0] << " " << quad.relativeSaliencies[1] << " " << quad.relativeSaliencies[2] << " "
			<< quad.localIntensity << " " << quad.firstEigenvalue << " " << quad.type << endl;
	}
	for (int i = 0; i < segments->size(); i++)
	{
		Segment segment = (*segments)[i];
		ofs << 2 << " " << segment.vertIdxs[0] << " " << segment.vertIdxs[1] << " "
		//ofs << 4 << " " << segment.vertIdxs[0] << " " << segment.vertIdxs[1] << " " << segment.vertIdxs[0] << " " << segment.vertIdxs[1] << " "
			<< segment.relativeSaliencies[0] << " " << segment.relativeSaliencies[1] << " " << segment.relativeSaliencies[2] << " "
			<< segment.localIntensity << " " << segment.firstEigenvalue << " " << segment.type << endl;
	}
	for (int i = 0; i < points->size(); i++)
	{
		Point point = (*points)[i];
		ofs << 1 << " " << point.vertIdx << " "
		//ofs << 4 << " " << point.vertIdx << " " << point.vertIdx << " " << point.vertIdx << " " << point.vertIdx << " "
			<< point.relativeSaliencies[0] << " " << point.relativeSaliencies[1] << " " << point.relativeSaliencies[2] << " "
			<< point.localIntensity << " " << point.firstEigenvalue << " " << point.type << endl;
	}
}

void General_Data::loadOFF(ifstream &ifs)
{
	char line[256], word[64];
	int vertNum, quadNum, segNum, pointNum;
	vertices = new vector<Vertex>;
	quads = new vector<Quad>;
	segments = new vector<Segment>;
	points = new vector<Point>;
	ifs.getline(line, 256);
	ifs.getline(line, 256);
	ifs >> vertNum >> quadNum >> word;
	ifs >> word >> segNum >> pointNum;
	ifs >> word >> maxIntensity >> minIntensity >> *maxEigenvalue >> *minEigenvalue;

	for (int i = 0; i < vertNum; i++)
	{
		Vertex vertex;
		ifs >> vertex.position[0] >> vertex.position[1] >> vertex.position[2];
		vertices->push_back(vertex);
	}
	for (int i = 0; i < quadNum; i++)
	{
		Quad quad;
		ifs >> word >> quad.vertIdxs[0] >> quad.vertIdxs[1] >> quad.vertIdxs[2] >> quad.vertIdxs[3]
			>> quad.relativeSaliencies[0] >> quad.relativeSaliencies[1] >> quad.relativeSaliencies[2]
			>> quad.localIntensity >> quad.firstEigenvalue >> quad.type;
		quads->push_back(quad);
	}
	for (int i = 0; i < segNum; i++)
	{
		Segment segment;
		ifs >> word >> segment.vertIdxs[0] >> segment.vertIdxs[1]
			>> segment.relativeSaliencies[0] >> segment.relativeSaliencies[1] >> segment.relativeSaliencies[2]
			>> segment.localIntensity >> segment.firstEigenvalue >> segment.type;
		segments->push_back(segment);
	}
	for (int i = 0; i < pointNum; i++)
	{
		Point point;
		ifs >> word >> point.vertIdx
			>> point.relativeSaliencies[0] >> point.relativeSaliencies[1] >> point.relativeSaliencies[2]
			>> point.localIntensity >> point.firstEigenvalue >> point.type;
		points->push_back(point);
	}
}

void General_Data::saveVPT(ofstream &ofs)
{
	ofs << view3D->currentMat[0] << " " << view3D->currentMat[1] << " " << view3D->currentMat[2] << " " << view3D->currentMat[3] << " ";
	ofs << view3D->currentMat[4] << " " << view3D->currentMat[5] << " " << view3D->currentMat[6] << " " << view3D->currentMat[7] << " ";
	ofs << view3D->currentMat[8] << " " << view3D->currentMat[9] << " " << view3D->currentMat[10] << " " << view3D->currentMat[11] << " ";
	ofs << view3D->currentMat[12] << " " << view3D->currentMat[13] << " " << view3D->currentMat[14] << " " << view3D->currentMat[15] << endl;
	int view[4];
	glGetIntegerv(GL_VIEWPORT, view);
	double model[16], proj[16], cx, cy, cz;
	glGetDoublev(GL_MODELVIEW_MATRIX, model);
	glGetDoublev(GL_PROJECTION_MATRIX, proj);
	gluUnProject( (view[2]-view[0])/2 , (view[3]-view[1])/2, 0, model, proj, view, &cx, &cy, &cz);
	ofs << cx << " " << cy << " " << cz << " ";
	ofs << cx-model[2] << " " << cy-model[6] << " " << cz-model[10] << " ";
	ofs << model[1] << " " << model[5] << " " << model[9] << " " << 20 << " " << endl;
}

void General_Data::loadVPT(ifstream &ifs)
{
	ifs >> view3D->currentMat[0] >> view3D->currentMat[1] >> view3D->currentMat[2] >> view3D->currentMat[3];
	ifs >> view3D->currentMat[4] >> view3D->currentMat[5] >> view3D->currentMat[6] >> view3D->currentMat[7];
	ifs >> view3D->currentMat[8] >> view3D->currentMat[9] >> view3D->currentMat[10] >> view3D->currentMat[11];
	ifs >> view3D->currentMat[12] >> view3D->currentMat[13] >> view3D->currentMat[14] >> view3D->currentMat[15];
	/*ifs >> view3D->lookAt[0] >> view3D->lookAt[1] >> view3D->lookAt[2] >> view3D->lookAt[3] >> view3D->lookAt[4]
	>> view3D->lookAt[5] >> view3D->lookAt[6] >> view3D->lookAt[7] >> view3D->lookAt[8];
	view3D->currentMat[0] = view3D->currentMat[5] = view3D->currentMat[10] = view3D->currentMat[15] = 1;
	view3D->currentMat[1] = view3D->currentMat[2] = view3D->currentMat[3] = 
		view3D->currentMat[4] = view3D->currentMat[6] = view3D->currentMat[7] = 
		view3D->currentMat[8] = view3D->currentMat[9] = view3D->currentMat[11] = 
		view3D->currentMat[12] = view3D->currentMat[13] = view3D->currentMat[14] = 0;*/
}

void General_Data::edgePhase(float *coef)
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Edge Phase Begin:" << endl;
	change = 0.017f * gridSize;
	maxEdgeIndex = (gridx-1) * gridy * gridz + gridx * (gridy-1) * gridz + gridx * gridy * (gridz-1);
	edges = new Edge[ maxEdgeIndex ] ;
	grid_2DGradMag = new float[ maxGridIndex ] ;
	storeGrid_2DGradMag = new bool[ maxGridIndex ] ;
	for ( int index = 0 ; index < maxGridIndex ; index ++ )
	{
		storeGrid_2DGradMag[index] = false;
	}
	totalEdgePoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	CGNSParallelEdge parallel_edge(gridx,gridy,gridz,totalEdgePoints,edges,allCubic,globalVec,
		dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,coef);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_edge,auto_partitioner());

	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::facePhase(float *coef)
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Face Phase Begin:" << endl;
	change = 0.017f * gridSize;
	maxFaceIndex = gridx * (gridy-1) * (gridz-1) + (gridx-1) * gridy * (gridz-1) + (gridx-1) * (gridy-1) * gridz;
	faces = new Face[ maxFaceIndex ] ;
	totalFacePoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	CGNSParallelFace parallel_face(gridx,gridy,gridz,totalFacePoints,edges,faces,dataType,
		gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,coef);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_face,auto_partitioner());

	delete[] grid_2DGradMag;
	delete[] storeGrid_2DGradMag;
	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}

void General_Data::cellPhase(float *coef)
{
	clock_t timeStart, timeEnd;
	timeStart = clock();
	cout << "Cell Phase Begin:" << endl;
	points = new vector<Point>;
	vertices = new vector<Vertex>;
	conpoints = new concurrent_vector<Point>;
	convertices = new concurrent_vector<Vertex>;
	change = 0.017f * gridSize;
	dispCellX = (int)floor( gridx / 2.0 ) - 1;
	dispCellY = (int)floor( gridy / 2.0 ) - 1;
	dispCellZ = (int)floor( gridz / 2.0 ) - 1;
	maxCellIndex = (gridx-1) * (gridy-1) * (gridz-1);
	cells = new Cell[ maxCellIndex ] ;
	totalCellPoints = new int(0);

	tbb::task_scheduler_init init(task_scheduler_init::automatic);  
	CGNSParallelCell parallel_cell(gridx,gridy,gridz,totalCellPoints,edges,faces,cells,conpoints,dataType,gridPoints,change,
		gridSize,maxEigenvalue, minEigenvalue,convertices,halfDataSize,rightBackTop,thickness,coef);
	parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_cell,auto_partitioner());

	concurrent_vector<Point>::iterator iter = conpoints->begin();
	while( iter != conpoints->end())
	{
		Point p;
		p.vertIdx = (*iter).vertIdx;
	    p.relativeSaliencies[0] = (*iter).relativeSaliencies[0];
		p.relativeSaliencies[1] = (*iter).relativeSaliencies[1];
		p.relativeSaliencies[2] = (*iter).relativeSaliencies[2];
	    p.localIntensity = (*iter).localIntensity;
		p.firstEigenvalue = (*iter).firstEigenvalue;
		p.type=(*iter).type;		
		points->push_back(p);
		iter++;
	}

	concurrent_vector<Vertex>::iterator iterv = convertices->begin();
	while( iterv != convertices->end())
	{
		Vertex p;
		p.position[0]=(*iterv).position[0];
		p.position[1]=(*iterv).position[1];
		p.position[2]=(*iterv).position[2];
		vertices->push_back(p);
		iterv++;
	}
	timeEnd = clock();
	cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
	cout << "Done!" << endl;
	cout << "****************************************************************" << endl;
}
