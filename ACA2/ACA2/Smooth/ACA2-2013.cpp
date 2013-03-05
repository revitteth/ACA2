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
#include <tbb/tbb.h>
#include <CL/cl.h>

int main(int argc, char **argv){
  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " mesh_file" << std::endl;
  }

  Mesh *mesh = new Mesh(argv[1]);
  Mesh *mesh_cl = new Mesh(argv[1]);

  Quality q = mesh->get_mesh_quality();
  Quality q_cl = mesh_cl->get_mesh_quality();

  std::cout << "Initial quality:\n"
            << "Quality mean:  " << q.mean << std::endl
            << "Quality min:   " << q.min << std::endl;
	std::cout << "Initial quality (CL):\n"
            << "Quality mean:  " << q_cl.mean << std::endl
            << "Quality min:   " << q_cl.min << std::endl;

  tbb::tick_count start_time = tbb::tick_count::now();
  smooth(mesh, 200);
  tbb::tick_count end_time = tbb::tick_count::now();

  double time_smooth = (end_time - start_time).seconds();

  tbb::tick_count start_time_cl = tbb::tick_count::now();
  smooth_cl(mesh_cl, 200);
  tbb::tick_count end_time_cl = tbb::tick_count::now();

  double time_smooth_cl = (end_time_cl - start_time_cl).seconds();


  q = mesh->get_mesh_quality();

  std::cout<<"After smoothing:\n"
           << "Quality mean:  " << q.mean << std::endl
           << "Quality min:   " << q.min << std::endl;
    std::cout<<"After smoothing (CL):\n"
           << "Quality mean:  " << q_cl.mean << std::endl
           << "Quality min:   " << q_cl.min << std::endl;

  if((q.mean>0.90)&&(q.min>0.55))
    std::cout << "Test passed"<< std::endl;
  else
    std::cout << "Test failed"<< std::endl;

  if((q_cl.mean>0.90)&&(q_cl.min>0.55))
    std::cout << "Test passed (CL)"<< std::endl;
  else
    std::cout << "Test failed (CL)"<< std::endl;

  std::cout<<"BENCHMARK: \t\t" << time_smooth << "s" << std::endl;
  std::cout<<"BENCHMARK (CL): \t" << time_smooth_cl << "s" << std::endl;

  delete mesh;

  std::cin.get();

  return EXIT_SUCCESS;
}
