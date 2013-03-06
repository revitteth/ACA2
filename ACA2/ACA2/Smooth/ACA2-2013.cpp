//============================================================================
// Name        : ACA2-2013.cpp
// Author      : George Rokos
// Description : 2nd Assessed Coursework for ACA 2013
//============================================================================

#include <cstdlib>

#include <time.h>
#include <iostream>
#include "Mesh.hpp"
#include "Smooth.hpp"
#include "Smooth_CL.hpp"
#include "Smooth_tbb.hpp"
#include "Smooth_timed.hpp"
#include <tbb/tbb.h>
#include "Timer.hpp"
#include "Information.hpp"

int main(int argc, char **argv){
	if(argc!=2){
	std::cerr << "Usage: " << argv[0] << " mesh_file" << std::endl;
	}

	//platformInfo();

	Mesh *mesh = new Mesh(argv[1]);
	Mesh *mesh_cl = new Mesh(argv[1]);
	Mesh *mesh_tbb = new Mesh(argv[1]);

	Quality q = mesh->get_mesh_quality();

	reportSmoothHeaders(q);

	#ifdef SMOOTH_HPP_
		Timer* t1 = new Timer(tbb::tick_count::now());
		//smooth(mesh, 200);
		t1->Stop(tbb::tick_count::now());
		reportSmooth(mesh, t1, "default");
		delete mesh;
	#endif /* SMOOTH_HPP_ */

	#ifdef SMOOTH_CL_HPP_
		Timer* t_cl = new Timer(tbb::tick_count::now());
		//smooth_cl(mesh_cl, 200);
		t_cl->Stop(tbb::tick_count::now());
		reportSmooth(mesh_cl, t_cl, "openCL1");
		delete mesh_cl;
	#endif /* SMOOTH_HPP_ */

	#ifdef SMOOTH_tbb_HPP_
		Timer* t_tbb = new Timer(tbb::tick_count::now());
		//smooth_tbb(mesh_tbb, 200);
		t_tbb->Stop(tbb::tick_count::now());
		reportSmooth(mesh_tbb, t_tbb, "tbb1");
		delete mesh_tbb;
	#endif /* SMOOTH_HPP_ */

	//For individual loop timing
	//smooth_timer_start(mesh_tbb, 200);



	std::cin.get();

	return EXIT_SUCCESS;
}
