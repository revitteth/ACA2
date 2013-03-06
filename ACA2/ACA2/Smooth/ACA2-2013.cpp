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
  //Quality q_cl = mesh_cl->get_mesh_quality(); unecessary as same as input parameter

  Timer* t1 = new Timer(tbb::tick_count::now());
  //smooth(mesh, 200);
  t1->Stop(tbb::tick_count::now());

  Timer* t_cl = new Timer(tbb::tick_count::now());
  //smooth_cl(mesh_cl, 200);
  t_cl->Stop(tbb::tick_count::now());

  Timer* t_tbb = new Timer(tbb::tick_count::now());
  //smooth_tbb(mesh_tbb, 200);
  t_tbb->Stop(tbb::tick_count::now());

  //For individual loop timing
  smooth_timer_start(mesh_tbb, 200);

  reportSmoothHeaders(q);
  reportSmooth(mesh, t1, "default");
  reportSmooth(mesh_cl, t_cl, "openCL1");
  reportSmooth(mesh_tbb, t_tbb, "tbb1");

  delete mesh;
  delete mesh_cl;
  delete mesh_tbb;

  std::cin.get();

  return EXIT_SUCCESS;
}
