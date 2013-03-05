//============================================================================
// Name        : ACA2-2013.cpp
// Author      : George Rokos
// Description : 2nd Assessed Coursework for ACA 2013
//============================================================================

#include <cstdlib>
#include <time.h>
#include <iostream>
#include <tbb/tbb.h>
#include <CL/cl.h>

#include "Mesh.hpp"
#include "Smooth.hpp"
#include "Smooth_CL.hpp"
#include "Timer.hpp"
#include "Information.hpp"

int main(int argc, char **argv){
  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " mesh_file" << std::endl;
  }

  //platformInfo();

  Mesh *mesh = new Mesh(argv[1]);
  Mesh *mesh_cl = new Mesh(argv[1]);

  Quality q = mesh->get_mesh_quality();

  Timer* t1 = new Timer(tbb::tick_count::now());
  smooth(mesh, 200);
  t1->Stop(tbb::tick_count::now());

  Timer* t_cl = new Timer(tbb::tick_count::now());
  smooth_cl(mesh_cl, 200);
  t_cl->Stop(tbb::tick_count::now());


  reportSmoothHeaders(q);
  reportSmooth(mesh, t1, "default");
  reportSmooth(mesh_cl, t_cl, "openCL1");

  delete mesh;
  delete mesh_cl;

  std::cin.get();

  return EXIT_SUCCESS;
}
