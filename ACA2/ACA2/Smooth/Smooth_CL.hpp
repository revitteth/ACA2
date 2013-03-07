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


void smooth_cl(Mesh* mesh, size_t niter){

		size_t nelist_size = mesh->NEList.size();
		size_t enlist_size = mesh->ENList.size();
		size_t nnlist_size = mesh->NNList.size();
		size_t coords_size = mesh->coords.size();
		size_t metric_size = mesh->metric.size();
		size_t normals_size = mesh->normals.size();
		size_t nnodes = mesh->NNodes;
		int orientation = mesh->get_orientation();
		size_t nelist_max_set_size = 0;

		// construct flat ne list with set size (stride) = nelist_max_set_size. Set all other values to max_size_t or something.

		for(unsigned i = 0; i < mesh->NEList.size(); i++)
		{
			if(nelist_max_set_size < mesh->NEList[i].size())
				nelist_max_set_size = mesh->NEList[i].size();
		}
		std::cout << nelist_max_set_size << std::endl;
		
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
        std::ifstream sourceFile("..\\ACA2\\Smooth\\nodes_solve_kernel.cl");
        std::string sourceCode(
            std::istreambuf_iterator<char>(sourceFile),
            (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

        // Make program of the source code in the context
        cl::Program program = cl::Program(context, source);
 
        // Build program for these specific devices
        program.build(devices);
 
        // Make kernel
        cl::Kernel kernel(program, "nodes_solve");

        // Create memory buffers (pointer and size in order to iterate easily)
		cl::Buffer buf_nelist =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, nelist_size * sizeof(size_t));					//nelist pointer
		cl::Buffer buf_nelist_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//nelist size
		cl::Buffer buf_enlist =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, enlist_size * sizeof(size_t));					//enlist pointer
		cl::Buffer buf_enlist_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//enlist size
		cl::Buffer buf_coords =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, coords_size * sizeof(size_t));
		cl::Buffer buf_coords_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//coords size
		cl::Buffer buf_metric =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, metric_size * sizeof(size_t));					//coords pointer
		cl::Buffer buf_metric_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//coords size
		cl::Buffer buf_normals =			cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, normals_size * sizeof(size_t));					//normals pointer
		cl::Buffer buf_normals_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//normals size
		cl::Buffer buf_orientation =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(int));									//orientation (value)

		// Copy to the memory buffers
		queue.enqueueWriteBuffer(buf_nelist, CL_TRUE, 0, nelist_size * sizeof(size_t), &mesh->NEList[0]);
		queue.enqueueWriteBuffer(buf_nelist_size, CL_TRUE, 0, sizeof(size_t), &nelist_size);
		queue.enqueueWriteBuffer(buf_enlist, CL_TRUE, 0, enlist_size * sizeof(size_t), &mesh->ENList);
		queue.enqueueWriteBuffer(buf_enlist_size, CL_TRUE, 0, sizeof(size_t), &enlist_size);
		queue.enqueueWriteBuffer(buf_coords, CL_TRUE, 0, coords_size * sizeof(size_t), &mesh->coords[0]);
		queue.enqueueWriteBuffer(buf_coords_size, CL_TRUE, 0, sizeof(size_t), &coords_size);
		queue.enqueueWriteBuffer(buf_metric, CL_TRUE, 0, metric_size * sizeof(size_t), &mesh->metric[0]);
		queue.enqueueWriteBuffer(buf_metric_size, CL_TRUE, 0, sizeof(size_t), &metric_size);
		queue.enqueueWriteBuffer(buf_normals, CL_TRUE, 0, normals_size * sizeof(size_t), &mesh->normals[0]);
		queue.enqueueWriteBuffer(buf_normals_size, CL_TRUE, 0, sizeof(size_t), &normals_size);
		queue.enqueueWriteBuffer(buf_orientation, CL_TRUE, 0, sizeof(int), &orientation);
 
		// Set arguments to kernel
		kernel.setArg(0,  buf_nelist);
		kernel.setArg(1,  buf_nelist_size);
		kernel.setArg(2,  buf_enlist);
		kernel.setArg(3,  buf_enlist_size);
		kernel.setArg(4,  buf_coords);
		kernel.setArg(5,  buf_coords_size);
		kernel.setArg(6,  buf_metric);
		kernel.setArg(7,  buf_metric_size);
		kernel.setArg(8,  buf_normals);
		kernel.setArg(9,  buf_normals_size);
		kernel.setArg(10, buf_orientation);
        
		// For the specified number of iterations, loop over all mesh vertices.
		for(size_t iter=0; iter<niter; ++iter)
		{

			// Run the kernel on specific ND range
			cl::NDRange global(nnodes);
			cl::NDRange local(1);
			queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
 
			// Read buffer C into a local list
			double *temp = new double[mesh->coords.size()];
			queue.enqueueReadBuffer(buf_coords, CL_TRUE, 0, coords_size * sizeof(size_t), temp);
			std::vector<double> tmp (temp, temp + mesh->coords.size());
			mesh->coords = tmp;
		}
 
    } 
	catch(cl::Error error) 
	{
		std::cout << error.what() << "(" << error.err() << ")" << std::endl;
    }
}

#endif
