//============================================================================
// Name        : Smooth.hpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel prototype
//============================================================================

#ifndef SMOOTH_HPP_timed
#define SMOOTH_HPP_timed

#include "Mesh.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

#include <vector>
#include <utility>


#include "SVD2x2.hpp"
#include "Smooth.hpp"
#include "Timer.hpp"


#include <iostream>
#include <iomanip>

void smooth_timed(Mesh* mesh, size_t niter, std::vector<double> &times, std::vector<double> &passes){
  // For the specified number of iterations, loop over all mesh vertices.
	  Timer* time_overall = new Timer(tbb::tick_count::now());
     Timer* time_loop_1 = new Timer(tbb::tick_count::now());
     Timer* time_loop_2 = new Timer(tbb::tick_count::now());
		Timer* general = new Timer(tbb::tick_count::now());
		  		  	Timer* for_loop_1 = new Timer(tbb::tick_count::now());


 for(size_t iter=0; iter<niter; ++iter){

																		time_loop_1->Start(tbb::tick_count::now());


    for(size_t vid=0; vid<mesh->NNodes; ++vid){

																			 time_loop_2->Start(tbb::tick_count::now());

	// If this is a corner node, it cannot be moved.
																				 general->Start(tbb::tick_count::now());
      if(mesh->isCornerNode(vid))
        continue;
																						general->Stop(tbb::tick_count::now());
																						times[3]+= general->getTime();
																						passes[3]++;
																					general->Start(tbb::tick_count::now());
	  // Find the quality of the worst element adjacent to vid
      double worst_q=1.0;
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin();
          it!=mesh->NEList[vid].end(); ++it){
        worst_q = min(worst_q, mesh->element_quality(*it));
      }
																						general->Stop(tbb::tick_count::now());
																						times[4]+= general->getTime();
																						passes[4]++;
																					general->Start(tbb::tick_count::now());

      /* Find the barycentre (centre of mass) of the cavity. A cavity is
       * defined as the set containing vid and all its adjacent vertices and
       * triangles. Since we work on metric space, all lengths have to measured
       * using the metric. The metric tensor is a 2x2 symmetric matrix which
       * encodes the ideal length and orientation of an edge containing vid. As
       * an edge is defined by two vertices, we calculate the edge length using
       * the value of the metric in the middle of the edge, i.e. the average of
       * the two metric tensors of the vertices defining the edge.
       */
      const double * m0 = &mesh->metric[3*vid];

      double x0 = mesh->coords[2*vid];
      double y0 = mesh->coords[2*vid+1];

      double A[4] = {0.0, 0.0, 0.0, 0.0};
      double q[2] = {0.0, 0.0};
																						general->Stop(tbb::tick_count::now());
																						times[5]+= general->getTime();
																						passes[5]++;
																					general->Start(tbb::tick_count::now());

      // Iterate over all edges and assemble matrices A and q.
      for(std::vector<size_t>::const_iterator it=mesh->NNList[vid].begin();

          it!=mesh->NNList[vid].end(); ++it){
																		//	  for_loop_1->Start(tbb::tick_count::now());
        size_t il = *it;

        const double *m1 = &mesh->metric[3*il];

        // Find the metric in the middle of the edge.
        double ml00 = 0.5*(m0[0] + m1[0]);
        double ml01 = 0.5*(m0[1] + m1[1]);
        double ml11 = 0.5*(m0[2] + m1[2]);

        double x = mesh->coords[2*il] - x0;
        double y = mesh->coords[2*il+1] - y0;

        // Calculate and accumulate the contribution of
        // this vertex to the barycentre of the cavity.
        q[0] += (ml00*x + ml01*y);
        q[1] += (ml01*x + ml11*y);

        A[0] += ml00;
        A[1] += ml01;
        A[3] += ml11;
																			//	for_loop_1->Stop(tbb::tick_count::now());
																			//	times[7]+= for_loop_1->getTime();
																			//	passes[7]++;

      }
																						general->Stop(tbb::tick_count::now());
																						times[6]+= general->getTime();
																						passes[6]++;
																					general->Start(tbb::tick_count::now());

      // The metric tensor is symmetric, i.e. ml01=ml10, so A[2]=A[1].
      A[2]=A[1];

      // Displacement vector for vid
      double p[2];

      /* The displacement p for vid is found by solving the linear system:
       * ┌─       ─┐   ┌    ┐   ┌    ┐
       * │A[0] A[1]│   │p[0]│   │q[0]│
       * │         │ x │    │ = │    │
       * │A[2] A[3]│   │p[1]│   │q[0]│
       * └─       ─┘   └    ┘   └    ┘
       */
      svd_solve_2x2(A, p, q);

      /* If this is a surface vertex, restrict the displacement
       * to the surface. The new displacement is the projection
       * of the old displacement on the surface.
       */
      if(mesh->isSurfaceNode(vid)){
        p[0] -= p[0]*fabs(mesh->normals[2*vid]);
        p[1] -= p[1]*fabs(mesh->normals[2*vid+1]);
      }

																						general->Stop(tbb::tick_count::now());
																						times[8]+= general->getTime();
																						passes[8]++;
																					general->Start(tbb::tick_count::now());
      // Update the coordinates
      mesh->coords[2*vid] += p[0];
      mesh->coords[2*vid+1] += p[1];
																						general->Stop(tbb::tick_count::now());
																						times[9]+= general->getTime();
																						passes[9]++;
																					general->Start(tbb::tick_count::now());
      /************************************************************************
       * At this point we must also interpolate the metric tensors from all   *
       * neighbouring vertices in order to calculate the new value of vid's   *
       * metric tensor at the new location. This is a quite complex procedure *
       * and has been omitted for simplicity of the exercise. A vertex will   *
       * always use its original metric tensor, no matter whether it has been *
       * relocated or not.                                                    *
       ************************************************************************/

      /* Find the quality of the worst element after smoothing. If an element
       * of the cavity was inverted, i.e. if vid was relocated outside the
       * interior convex hull of the cavity, then the calculated area of that
       * element will be negative and mesh->element_quality() will return a
       * negative number. In such a case, the smoothing operation has to be
       * rejected.
       */
      double new_worst_q=1.0;
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin();
          it!=mesh->NEList[vid].end(); ++it){
        new_worst_q = min(new_worst_q, mesh->element_quality(*it));
      }

      /* If quality is worse than before, either because of element inversion
       * or just because relocating vid to the barycentre of the cavity does
       * not improve quality, revert the changes.
       */
      if(new_worst_q < worst_q){
        mesh->coords[2*vid] -= p[0];
        mesh->coords[2*vid+1] -= p[1];
      }

																					general->Stop(tbb::tick_count::now());
																					times[10]+= general->getTime();
																					passes[10]++;

	  
																					time_loop_2->Stop(tbb::tick_count::now());
																					times[2]+= time_loop_2->getTime();
																					passes[2]++;



    }

																					time_loop_1->Stop(tbb::tick_count::now());
																					times[1]+= time_loop_1->getTime();
																					passes[1]++;



  }

																					time_overall->Stop(tbb::tick_count::now());
																					times[0] = time_overall->getTime();
																					passes[0]++;
		 

}


void smooth_timer_start(Mesh* mesh, size_t niter) {
	int no = 11;
	//	std::vector<Timer> timings(no);
		std::vector<double> times(no, 0);
		std::vector<double> passes(no, 0);

	  //Timer* time_overall = new Timer(tbb::tick_count::now());
	  smooth_timed(mesh, niter, times, passes);
	  //time_overall->Stop(tbb::tick_count::now());
	  //std::cout << "time_overall" << time_overall->getTime();

	  for (int i = 0; i<times.size(); i++) {
		  std::cout << " " << std::setw(12)  << times[i];
		  std::cout << " " << std::setw(12)  << passes[i];
		  std::cout << " " << std::setw(20)  << times[i]/passes[i];
		  std::cout << " \n ";
	  }


}

#endif /* SMOOTH_HPP_ */
