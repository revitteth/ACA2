
__kernel void nodes_solve(
	__global size_t* nelist, 
	__global size_t* nelist_size,
	__global size_t* nelist_set_size,
	__global size_t* nnlist,
	__global size_t* nnlist_size,
	__global size_t* nnlist_set_size,
	__global size_t* enlist,
	__global size_t* enlist_size,
	__global float* coords,
	__global size_t* coords_size,
	__global float* metric,
	__global size_t* metric_size,
	__global float* normals,
	__global size_t* normals_size,
	__global int* orientation
	)
{

    const int vid = get_global_id(0);

	//coords[vid] = *nelist_size;

	// if it is a corner node
	if (normals[2*vid] == 1 | normals[2*vid] == -1)
	{
		if (normals[2*vid+1] == 1 | normals[2*vid+1] == -1)
		{
			return;
		}
	}

	float worst_q = 1.0;

	for (unsigned i = vid*(*nelist_set_size); i < (vid+1)*(*nelist_set_size); i++)
	{
		size_t* n = enlist[3*i];

		float m00 = (metric[3*n[0]] + metric[3*n[1]] + metric[3*n[2]]/3);
		float m01 = (metric[3*n[0]+1] + metric[3*n[1]+1] + metric[3*n[2]+1]/3);
		float m11 = (metric[3*n[0]+2] + metric[3*n[1]+2] + metric[3*n[2]+1]/3);

		float l =
			sqrt((coords[2*n[0]+1] - coords[2*n[1]+1])*((coords[2*n[0]+1] - coords[2*n[1]+1])*m11 + (coords[2*n[0]] - coords[2*n[1]])*m01) +
				 (coords[2*n[0]] - coords[2*n[1]])*((coords[2*n[0]+1] - coords[2*n[1]+1])*m01 + (coords[2*n[0]] - coords[2*n[1]])*m00))+
			sqrt((coords[2*n[0]+1] - coords[2*n[2]+1])*((coords[2*n[0]+1] - coords[2*n[2]+1])*m11 + (coords[2*n[0]] - coords[2*n[2]])*m01) +
				 (coords[2*n[0]] - coords[2*n[2]])*((coords[2*n[0]+1] - coords[2*n[2]+1])*m01 + (coords[2*n[0]] - coords[2*n[2]])*m00))+
			sqrt((coords[2*n[2]+1] - coords[2*n[1]+1])*((coords[2*n[2]+1] - coords[2*n[1]+1])*m11 + (coords[2*n[2]] - coords[2*n[1]])*m01) +
				 (coords[2*n[2]] - coords[2*n[1]])*((coords[2*n[2]+1] - coords[2*n[1]+1])*m01 + (coords[2*n[2]] - coords[2*n[1]])*m00));

		float a =
			*orientation * 0.5 *
            ( (coords[2*n[0]+1] - coords[2*n[2]+1]) * (coords[2*n[0]] - coords[2*n[1]]) -
              (coords[2*n[0]+1] - coords[2*n[1]+1]) * (coords[2*n[0]] - coords[2*n[2]]) );

		float a_m = a*sqrt(m00*m11 - m01*m01);

		float f = min(l/3.0, 3.0/l);
		float inter = f * (2.0-f);
		float F = inter * inter * inter;

		float quality = 12.0 * sqrt(3.0) * a_m * F / (l*l);

		worst_q = min(worst_q, quality);
	}






}