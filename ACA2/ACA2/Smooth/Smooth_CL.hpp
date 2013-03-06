//============================================================================
// Name        : Smooth.hpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel prototype
//============================================================================

#ifndef SMOOTH_CL_HPP_
#define SMOOTH_CL_HPP_

#define __CL_ENABLE_EXCEPTIONS
#define __NO_STD_VECTOR // Use cl::vector instead of STL version

#include <CL/cl.hpp>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

#include "Mesh.hpp"
#include <algorithm>
#include <cmath>
#include "SVD2x2.hpp"

//using namespace cl;

void setupOpenCl()
{
	// Create the two input vectors
    const int LIST_SIZE = 1000;
    int *A = new int[LIST_SIZE]; 
    int *B = new int[LIST_SIZE];
    for(int i = 0; i < LIST_SIZE; i++) {
        A[i] = i;
        B[i] = LIST_SIZE - i;
    }
 
   try { 
        // Get available platforms
        cl::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
 
        // Select the default platform and create a context using this platform and the GPU
        cl_context_properties cps[3] = { 
            CL_CONTEXT_PLATFORM, 
            (cl_context_properties)(platforms[0])(), 
            0 
        };
        cl::Context context(CL_DEVICE_TYPE_GPU, cps);
 
        // Get a list of devices on this platform
        cl::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
 
        // Create a command queue and use the first device
        cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
 
        // Read source file
        std::ifstream sourceFile("..\\ACA2\\Smooth\\vector_add_kernel.cl");
        std::string sourceCode(
            std::istreambuf_iterator<char>(sourceFile),
            (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));
 
        // Make program of the source code in the context
        cl::Program program = cl::Program(context, source);
 
        // Build program for these specific devices
        program.build(devices);
 
        // Make kernel
        cl::Kernel kernel(program, "vector_add");
 
        // Create memory buffers
        cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_ONLY, LIST_SIZE * sizeof(int));
        cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_ONLY, LIST_SIZE * sizeof(int));
        cl::Buffer bufferC = cl::Buffer(context, CL_MEM_WRITE_ONLY, LIST_SIZE * sizeof(int));
 
        // Copy lists A and B to the memory buffers
        queue.enqueueWriteBuffer(bufferA, CL_TRUE, 0, LIST_SIZE * sizeof(int), A);
        queue.enqueueWriteBuffer(bufferB, CL_TRUE, 0, LIST_SIZE * sizeof(int), B);
 
        // Set arguments to kernel
        kernel.setArg(0, bufferA);
        kernel.setArg(1, bufferB);
        kernel.setArg(2, bufferC);
 
        // Run the kernel on specific ND range
        cl::NDRange global(LIST_SIZE);
        cl::NDRange local(1);
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
 
        // Read buffer C into a local list
        int *C = new int[LIST_SIZE];
        queue.enqueueReadBuffer(bufferC, CL_TRUE, 0, LIST_SIZE * sizeof(int), C);
 
        for(int i = 0; i < LIST_SIZE; i ++)
             std::cout << A[i] << " + " << B[i] << " = " << C[i] << std::endl; 
    } catch(cl::Error error) {
       std::cout << error.what() << "(" << error.err() << ")" << std::endl;
    }
}


double clWorstQuality(
	size_t* vid, 
	std::vector<size_t>* ENList, 
	std::vector< std::set<size_t> >* NEList, 
	std::vector<double>* metric, 
	std::vector<double>* coords, 
	int orientation)
{ 
   try { 
        // Get available platforms
        cl::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
 
        // Select the default platform and create a context using this platform and the GPU
        cl_context_properties cps[3] = { 
            CL_CONTEXT_PLATFORM, 
            (cl_context_properties)(platforms[0])(), 
            0 
        };
        cl::Context context(CL_DEVICE_TYPE_GPU, cps);
 
        // Get a list of devices on this platform
        cl::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
 
        // Create a command queue and use the first device
        cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
 
        // Read source file
        std::ifstream sourceFile("..\\ACA2\\Smooth\\worst_q.cl");
        std::string sourceCode(
            std::istreambuf_iterator<char>(sourceFile),
            (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));
 
        // Make program of the source code in the context
        cl::Program program = cl::Program(context, source);
 
        // Build program for these specific devices
        program.build(devices);
 
        // Make kernel
        cl::Kernel kernel(program, "worst_q");
 
        // Create memory buffers
        cl::Buffer buf_vid = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(vid));
		cl::Buffer buf_ENList = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(ENList));
		cl::Buffer buf_NEList = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(NEList));
		cl::Buffer buf_metric = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(metric));
		cl::Buffer buf_coords = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(coords));
		cl::Buffer buf_orientation = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(orientation));
		cl::Buffer buf_worst_q = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(double));

 
        // Copy to memory buffers
        queue.enqueueWriteBuffer(buf_vid, CL_TRUE, 0, sizeof(vid), vid);
		queue.enqueueWriteBuffer(buf_ENList, CL_TRUE, 0, sizeof(ENList), ENList);
		queue.enqueueWriteBuffer(buf_NEList, CL_TRUE, 0, sizeof(NEList), NEList);
		queue.enqueueWriteBuffer(buf_metric, CL_TRUE, 0, sizeof(metric), metric);
		queue.enqueueWriteBuffer(buf_coords, CL_TRUE, 0, sizeof(coords), coords);
		queue.enqueueWriteBuffer(buf_orientation, CL_TRUE, 0, sizeof(orientation), &orientation);
 
        // Set arguments to kernel
        kernel.setArg(0, buf_vid);
        kernel.setArg(1, buf_ENList);
		kernel.setArg(2, buf_NEList);
		kernel.setArg(3, buf_metric);
		kernel.setArg(4, buf_coords);
		kernel.setArg(5, buf_orientation);
		kernel.setArg(6, buf_worst_q);
 
        // Run the kernel on specific ND range
        cl::NDRange global(sizeof(NEList[(*vid)]));
        cl::NDRange local(5); // how many to split the list into ? 
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
 
        // Read buffer C into a local list
        double *worst_q = new double;
        queue.enqueueReadBuffer(buf_worst_q, CL_TRUE, 0, sizeof(double), worst_q);
		return (*worst_q);

    } catch(cl::Error error) {
       std::cout << error.what() << "(" << error.err() << ")" << std::endl;
    }
}

void smooth_cl(Mesh* mesh, size_t niter){
	// For the specified number of iterations, loop over all mesh vertices.
	for(size_t iter=0; iter<niter; ++iter)
	{
		for(size_t vid = 0; vid < mesh->NNodes; ++vid)
		{
			// If this is a corner node, it cannot be moved.
			if(mesh->isCornerNode(vid))
				continue;

			// Find the quality of the worst element adjacent to vid
			double worst_q=1.0;
			worst_q = clWorstQuality(&vid, &mesh->ENList, &mesh->NEList, &mesh->metric, &mesh->coords, mesh->get_orientation());
			//for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin(); it!=mesh->NEList[vid].end(); ++it)
			//{
			//	worst_q = min(worst_q, mesh->element_quality(*it));
			//}

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
			if(mesh->isSurfaceNode(vid))
			{
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
			if(new_worst_q < worst_q)
			{
				mesh->coords[2*vid] -= p[0];
				mesh->coords[2*vid+1] -= p[1];
			}
		}
	}
}

#endif
