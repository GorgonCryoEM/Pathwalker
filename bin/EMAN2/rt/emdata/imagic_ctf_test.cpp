/*
 * Author: Grant Tang, 11/04/2008 (gtang@bcm.edu)
 * Copyright (c) 2000-2007 Baylor College of Medicine
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
 * */
#include <iostream>
#include "emdata.h"
#include "ctf.h"
using namespace std;
using namespace EMAN;

int main()
{
	cout << "Start Imagec Ctf test..." << endl;

	EMData * img = new EMData();
	img->read_image("e1_test.hed", 0);
	Ctf * ctf_ = img->get_attr("ctf");

	cout << "ctf_->to_string(): " << ctf_->to_string() << endl;

	delete img;
	delete ctf_;
	return 0;
}
