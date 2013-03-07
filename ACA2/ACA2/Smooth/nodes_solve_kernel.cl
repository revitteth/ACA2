int isCornerNode(size_t vid, size_t *normals)
{
	if(normals[2*vid] == 1 | normals[2*vid] == -1)
	{
		if(normals[2*vid+1] == 1 | normals[2*vid+1] == -1)
		{
		return 1;
		}
	}
	else
	{
		return 0;
	}
}

float element_quality(
	size_t eid, 
	size_t *enlist,
	float *coords,
	float *metric)
{
	//const size_t *n = &enlist[3*eid];

	//const float *c0 = &coords[2*n[0]];
	//const float *c1 = &coords[2*n[1]];
	//const float *c2 = &coords[2*n[2]];
}

__kernel void nodes_solve(
	__global size_t* nelist, 
	__global size_t* nelist_size,
	__global size_t* enlist,
	__global size_t* enlist_size,
	__global float* coords,
	__global size_t* coords_size,
	__global float* metric,
	__global size_t* metric_size,
	__global size_t* normals,
	__global size_t* normals_size,
	__global int* orientation
	)
{

    const int vid = get_global_id(0);

	coords[vid] = *nelist_size;

	if(isCornerNode(vid, *normals) == 1)
		return;

	float worst_q = 1.0;

	for(size_t i = 0; i < *nelist_size; i++)
	{
		worst_q = min(worst_q, 2.0); //implement element_quality!!!
		element_quality(i, *enlist, *coords, *metric);
	}
		



}