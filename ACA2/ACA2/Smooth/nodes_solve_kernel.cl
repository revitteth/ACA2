 #pragma OPENCL EXTENSION cl_khr_fp64: enable

__kernel void nodes_solve(
	__global size_t* nelist, 
	__global size_t* nelist_size,
	__global size_t* enlist,
	__global size_t* enlist_size,
	__global size_t* coords,
	__global size_t* coords_size,
	__global size_t* metric,
	__global size_t* metric_size,
	__global size_t* normals,
	__global size_t* normals_size,
	__global int* orientation
	)
{

    const int vid = get_global_id(0);

	coords[vid] = nelist_size;
}