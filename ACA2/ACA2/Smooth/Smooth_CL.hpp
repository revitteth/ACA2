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

		// NELIST - nelist_cl, nelist_size, nelist_set_size
		size_t nelist_size = mesh->NEList.size();
		size_t nelist_set_size = 0;
		for(unsigned i = 0; i < mesh->NEList.size(); i++)
		{
			if(nelist_set_size < mesh->NEList[i].size())
				nelist_set_size = mesh->NEList[i].size();
		}
		size_t* nelist_cl = new size_t[nelist_size * nelist_set_size];
		std::vector<size_t> ne_vec;
		for(int i = 0; i < nelist_size; i++)
		{
			ne_vec.clear();
			std::copy(mesh->NEList[i].begin(), mesh->NEList[i].end(), std::back_inserter(ne_vec));
			for(int j = 0; j < nelist_set_size; j++)
			{
				if(j <= mesh->NEList[i].size())
					nelist_cl[i+j] = ne_vec[j];
				else
					nelist_cl[i+j] = 0;
			}
		}

		// nnlist - nnlist_cl, nnlist_size, nnlist_set_size
		size_t nnlist_size = mesh->NNList.size();
		size_t nnlist_set_size = 0;
		for(unsigned i = 0; i < mesh->NNList.size(); i++)
		{
			if(nnlist_set_size < mesh->NNList[i].size())
				nnlist_set_size = mesh->NNList[i].size();
		}
		size_t* nnlist_cl = new size_t[nnlist_size * nnlist_set_size];
		std::vector<size_t> nn_vec;
		for(int i = 0; i < nnlist_size; i++)
		{
			nn_vec.clear();
			std::copy(mesh->NNList[i].begin(), mesh->NNList[i].end(), std::back_inserter(nn_vec));
			for(int j = 0; j < nnlist_set_size; j++)
			{
				if(j <= mesh->NNList[i].size())
					nnlist_cl[i+j] = nn_vec[j];
				else
					nnlist_cl[i+j] = 0;
			}
		}

		// ENLIST - enlist_cl, enlist_size
		size_t enlist_size = mesh->ENList.size();
		size_t* enlist_cl = &mesh->ENList[0];

		// COORDS - coords_cl, coords_size
		size_t coords_size = mesh->coords.size();
		float* coords_cl = new float[coords_size];
		for(std::vector<double>::size_type i = 0; i != coords_size; i++)
		{
			coords_cl[i] = (float)mesh->coords[i];
		}

		// METRIC - metric_cl, metric_size
		size_t metric_size = mesh->metric.size();
		float* metric_cl = new float[metric_size];
		for(std::vector<double>::size_type i = 0; i != metric_size; i++)
		{
			metric_cl[i] = (float)mesh->coords[i];
		}

		// NORMALS - normals_cl, normals_size
		size_t normals_size = mesh->normals.size();
		float* normals_cl = new float[normals_size];
		for(std::vector<double>::size_type i = 0; i != normals_size; i++)
		{
			normals_cl[i] = (float)mesh->normals[i];
		}

		// ORIENTATION - orientation
		int orientation = mesh->get_orientation();

		
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
		cl::Buffer buf_nelist =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, nelist_size * nelist_set_size * sizeof(size_t));	//nelist pointer
		cl::Buffer buf_nelist_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//nelist size	
		cl::Buffer buf_nelist_set_size =	cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//nelistset size
		cl::Buffer buf_nnlist =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, nnlist_size * nnlist_set_size * sizeof(size_t));	//nnlist pointer
		cl::Buffer buf_nnlist_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//nnlist size
		cl::Buffer buf_nnlist_set_size =	cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//nelistset size
		cl::Buffer buf_enlist =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, enlist_size * sizeof(size_t));					//enlist pointer
		cl::Buffer buf_enlist_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//enlist size
		cl::Buffer buf_coords =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, coords_size * sizeof(float));
		cl::Buffer buf_coords_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//coords size
		cl::Buffer buf_metric =				cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, metric_size * sizeof(float));					//coords pointer
		cl::Buffer buf_metric_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//coords size
		cl::Buffer buf_normals =			cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, normals_size * sizeof(float));					//normals pointer
		cl::Buffer buf_normals_size =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(size_t));									//normals size
		cl::Buffer buf_orientation =		cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(int));									//orientation (value)

		// Copy to the memory buffers
		queue.enqueueWriteBuffer(buf_nelist, CL_TRUE, 0, nelist_size * nelist_set_size * sizeof(size_t), &nelist_cl[0]);
		queue.enqueueWriteBuffer(buf_nelist_size, CL_TRUE, 0, sizeof(size_t), &nelist_size);
		queue.enqueueWriteBuffer(buf_nelist_set_size, CL_TRUE, 0, sizeof(size_t), &nelist_set_size);
		queue.enqueueWriteBuffer(buf_nnlist, CL_TRUE, 0, nnlist_size * nnlist_set_size * sizeof(size_t), &nnlist_cl[0]);
		queue.enqueueWriteBuffer(buf_nnlist_size, CL_TRUE, 0, sizeof(size_t), &nnlist_size);
		queue.enqueueWriteBuffer(buf_nnlist_set_size, CL_TRUE, 0, sizeof(size_t), &nnlist_set_size);
		queue.enqueueWriteBuffer(buf_enlist, CL_TRUE, 0, enlist_size * sizeof(size_t), &enlist_cl[0]);
		queue.enqueueWriteBuffer(buf_enlist_size, CL_TRUE, 0, sizeof(size_t), &enlist_size);
		queue.enqueueWriteBuffer(buf_coords, CL_TRUE, 0, coords_size * sizeof(float), &coords_cl[0]);
		queue.enqueueWriteBuffer(buf_coords_size, CL_TRUE, 0, sizeof(size_t), &coords_size);
		queue.enqueueWriteBuffer(buf_metric, CL_TRUE, 0, metric_size * sizeof(float), &metric_cl[0]);
		queue.enqueueWriteBuffer(buf_metric_size, CL_TRUE, 0, sizeof(size_t), &metric_size);
		queue.enqueueWriteBuffer(buf_normals, CL_TRUE, 0, normals_size * sizeof(float), &normals_cl[0]);
		queue.enqueueWriteBuffer(buf_normals_size, CL_TRUE, 0, sizeof(size_t), &normals_size);
		queue.enqueueWriteBuffer(buf_orientation, CL_TRUE, 0, sizeof(int), &orientation);
 
		// Set arguments to kernel
		kernel.setArg(0,  buf_nelist);
		kernel.setArg(1,  buf_nelist_size);
		kernel.setArg(2,  buf_nelist_set_size);
		kernel.setArg(3,  buf_nnlist);
		kernel.setArg(4,  buf_nnlist_size);
		kernel.setArg(5,  buf_nnlist_set_size);
		kernel.setArg(6,  buf_enlist);
		kernel.setArg(7,  buf_enlist_size);
		kernel.setArg(8,  buf_coords);
		kernel.setArg(9,  buf_coords_size);
		kernel.setArg(10, buf_metric);
		kernel.setArg(11, buf_metric_size);
		kernel.setArg(12, buf_normals);
		kernel.setArg(13, buf_normals_size);
		kernel.setArg(14, buf_orientation);
        
		// For the specified number of iterations, loop over all mesh vertices.
		for(size_t iter=0; iter<niter; ++iter)
		{
			// Run the kernel on specific ND range
			cl::NDRange global(mesh->NNodes);
			cl::NDRange local(1);
			queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
 
			// Read buffer C into a local list
			float *temp = new float[mesh->coords.size()];
			queue.enqueueReadBuffer(buf_coords, CL_TRUE, 0, coords_size * sizeof(size_t), temp);
			std::vector<float> tmp (temp, temp + mesh->coords.size());
			//mesh->coords = tmp;
		}
 
    } 
	catch(cl::Error error) 
	{
		std::cout << error.what() << "(" << error.err() << ")" << std::endl;
    }
}

#endif
