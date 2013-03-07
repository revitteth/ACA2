 #pragma OPENCL EXTENSION cl_khr_fp64: enable

__kernel void nodes_solve(
	__global size_t* nelist, 
	__global size_t* nelist_size,
	__global size_t* coords)
{

    const int vid = get_global_id(0);

	coords[vid] = nelist_size;
}