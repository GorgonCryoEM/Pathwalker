/**
 * $Id$
 */

/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */
#ifdef _WIN32
#pragma warning(disable:4819)
#include <malloc.h>
#endif	//_WIN32

#include <cstring>
#include <ctime>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "emdata.h"
#include "util.h"
#include "fundamentals.h"
#include "lapackblas.h"
#include "lbfgsb.h"
using namespace EMAN;
#include "steepest.h"
#include "emassert.h"
#include "randnum.h"
#include "mcqd.h"
#include <algorithm>

#include <cmath>
//#include <omp.h>
using namespace std;
using std::complex;

/* Subroutine */ 

int Util::coveig(int n, float *covmat, float *eigval, float *eigvec)
{
	// n size of the covariance/correlation matrix
	// covmat --- covariance/correlation matrix (n by n)
	// eigval --- returns eigenvalues
	// eigvec --- returns eigenvectors

	//ENTERFUNC;

	int i;

	// make a copy of covmat so that it will not be overwritten
	for ( i = 0 ; i < n * n ; i++ )   eigvec[i] = covmat[i];

	char NEEDV = 'V';
	char UPLO = 'U';
	int lwork = -1;
	int info = 0;
	float *work, wsize;

	//  query to get optimal workspace
	ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, &wsize, &lwork, &info);
	lwork = (int)wsize;

	work = (float *)calloc(lwork, sizeof(float));
	//  calculate eigs
	ssyev_(&NEEDV, &UPLO, &n, eigvec, &n, eigval, work, &lwork, &info);
	free(work);
	//EXITFUNC;
	return info;
}
