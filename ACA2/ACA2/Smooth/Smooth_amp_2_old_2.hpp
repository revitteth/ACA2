//============================================================================
// Name        : Smooth.hpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel prototype
//============================================================================

#ifndef SMOOTH_amp_2_HPP_
#define SMOOTH_amp_2_HPP_

#include "Mesh.hpp"
#include <algorithm>
#include <cmath>
#include <amp.h>

#include "SVD2x2.hpp"
#include "Smooth.hpp"

 #include <string>
#include <iostream>
#include <iomanip>


double JW_element_quality(Mesh* mesh, size_t eid) restrict(amp) {
    const size_t *n = &mesh->ENList[3*eid];


  // Pointers to the coordinates of each vertex
  const double *c0 = &mesh->coords[2*n[0]];
  const double *c1 = &mesh->coords[2*n[1]];
  const double *c2 = &mesh->coords[2*n[2]];

  // Pointers to the metric tensor at each vertex
  const double *m0 = &mesh->metric[3*n[0]];
  const double *m1 = &mesh->metric[3*n[1]];
  const double *m2 = &mesh->metric[3*n[2]];

  // Metric tensor averaged over the element
  double m00 = (m0[0] + m1[0] + m2[0])/3;
  double m01 = (m0[1] + m1[1] + m2[1])/3;
  double m11 = (m0[2] + m1[2] + m2[2])/3;

  // l is the length of the perimeter, measured in metric space
  double l =
    sqrt((c0[1] - c1[1])*((c0[1] - c1[1])*m11 + (c0[0] - c1[0])*m01) +
         (c0[0] - c1[0])*((c0[1] - c1[1])*m01 + (c0[0] - c1[0])*m00))+
    sqrt((c0[1] - c2[1])*((c0[1] - c2[1])*m11 + (c0[0] - c2[0])*m01) +
         (c0[0] - c2[0])*((c0[1] - c2[1])*m01 + (c0[0] - c2[0])*m00))+
    sqrt((c2[1] - c1[1])*((c2[1] - c1[1])*m11 + (c2[0] - c1[0])*m01) +
         (c2[0] - c1[0])*((c2[1] - c1[1])*m01 + (c2[0] - c1[0])*m00));

  // Area in physical space
  // Pointers to the coordinates of each vertex
  /*const double *c0 = &coords[2*n[0]];
  const double *c1 = &coords[2*n[1]];
  const double *c2 = &coords[2*n[2]];*/

  // DOUBLE ELEMENT AREA CONDENSED
  double a =  mesh->get_orientation() * 0.5 *
            ( (c0[1] - c2[1]) * (c0[0] - c1[0]) -
              (c0[1] - c1[1]) * (c0[0] - c2[0]) );

  // Area in metric space
  double a_m = a*sqrt(m00*m11 - m01*m01);

  // Function
  double f = min(l/3.0, 3.0/l);
  double F = pow(f * (2.0 - f), 3.0);

  // This is the 2D Lipnikov functional.
  double quality = 12.0 * sqrt(3.0) * a_m * F / (l*l);

  return quality;
}


void smooth_amp_2(Mesh* mesh, size_t niter){
  // For the specified number of iterations, loop over all mesh vertices.
  for(size_t iter=0; iter<niter; ++iter){
	concurrency::parallel_for((size_t)0, (size_t)mesh->NNodes, [&](size_t vid)
    { 
    //for(size_t vid=0; vid<mesh->NNodes; ++vid){
      // If this is a corner node, it cannot be moved.
      if(mesh->isCornerNode(vid))
        //continue;
        return; 

      // Find the quality of the worst element adjacent to vid
      double worst_q=1.0;
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin(); it!=mesh->NEList[vid].end(); ++it)
	  {
        //worst_q = min(worst_q, mesh->element_quality(*it));
		  //const size_t *n = &mesh->ENList[3*(*it)];
		  //std::cout << JW_element_quality(mesh, (*it)) << "\n";
		 worst_q = min(worst_q, JW_element_quality(mesh, (*it)));
      }

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

      // Iterate over all edges and assemble matrices A and q.
      for(std::vector<size_t>::const_iterator it=mesh->NNList[vid].begin(); it!=mesh->NNList[vid].end(); ++it)
	  {
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
      }

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

      // Update the coordinates
      mesh->coords[2*vid] += p[0];
      mesh->coords[2*vid+1] += p[1];

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
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin(); it!=mesh->NEList[vid].end(); ++it)
	  {
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
    });
  }
}

#endif /* SMOOTH_HPP_ */
