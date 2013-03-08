
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



	float m0 = metric[3*vid];

	float x0 = coords[2*vid];
	float y0 = coords[2*vid+1];

	float A[4] = {0.0, 0.0, 0.0, 0.0};
	float q[2] = {0.0, 0.0};

	for (unsigned i = vid*(*nelist_set_size); i < (vid+1)*(*nelist_set_size); i++)
	{
		float m1 = metric[3*i];

		float ml00 = 0.5*(metric[3*vid] + metric[3*i]);
		float ml01 = 0.5*(metric[3*vid+1] + metric[3*i+1]);
		float ml11 = 0.5*(metric[3*vid+2] + metric[3*i+2]);

		float x = coords[2*i] - x0;
		float y = coords[2*i+1] - y0;

		q[0] += (ml00*x + ml01*y);
		q[1] += (ml01*x + ml11*y);

		A[0] += ml00;
		A[1] += ml01;
		A[3] += ml11;
	}

	A[2]=A[1];

	float p[2];

	//svd_solve_2x2(A, p, q);
	float AAT[4];
	float eigenvalues[2];
	float U[4];
	float V[4];

	AAT[0] = A[0]*A[0] + A[2]*A[2];
	AAT[1] = A[0]*A[1] + A[2]*A[3];
	AAT[2] = AAT[1];
	AAT[3] = A[1]*A[1] + A[3]*A[3];

	//calc_eigenvals(AAT, eigenvals)
	float b, discriminant;

	b = A[0] + A[3];

	discriminant = sqrt(b*b - 4*(A[0]*A[3] - A[1]*A[2]));

	eigenvalues[0] = (b + discriminant) * 0.5;
	eigenvalues[1] = (b - discriminant) * 0.5;

	//calc_eigenvects(AAT, eigenvects, v)
	float eigenvectors[4];
	float D[4];
	float proportion;

	D[1] = AAT[1];
	D[2] = AAT[2];

	D[0] = AAT[0] - eigenvalues[0];
	D[3] = AAT[3] - eigenvalues[0];

	if (D[1] != 0.0)
	{
		proportion = -(D[0]/D[1]);
		eigenvectors[0] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[2] = proportion * eigenvectors[0];
	}
	else if (D[3] != 0.0)
	{
		proportion = -(D[2]/D[3]);
		eigenvectors[0] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[2] = proportion * eigenvectors[0];
	}
	else if (D[0] == 0.0)
	{
		eigenvectors[0] = 1.0;
		eigenvectors[2] = 0.0;
	}
	else if (D[2] == 0.0)
	{
		eigenvectors[2] = 1.0;
		eigenvectors[0] = 0.0;
	}

	D[0] = AAT[0] - eigenvalues[1];
	D[3] = AAT[3] - eigenvalues[1];

	if (D[1] != 0.0)
	{
		proportion = -(D[0]/D[1]);
		eigenvectors[1] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[3] = proportion * eigenvectors[1];
	}
	else if (D[3] != 0.0)
	{
		proportion = -(D[2]/D[3]);
		eigenvectors[1] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[3] = proportion * eigenvectors[1];
	}
	else if (D[0] == 0.0)
	{
		eigenvectors[1] = 1.0;
		eigenvectors[3] = 0.0;
	}
	else if (D[2] == 0.0)
	{
		eigenvectors[3] = 1.0;
		eigenvectors[1] = 0.0;
	}

	//back to svd
	if(eigenvalues[0] < 1E-12)
		eigenvalues[0] = 0.0;
	else
		eigenvalues[0] = 1.0 / sqrt(eigenvalues[0]);
	if(eigenvalues[1] < 1E-12)
		eigenvalues[1] = 0.0;
	else
		eigenvalues[1] = 1.0 / sqrt(eigenvalues[1]);

	U[0] = eigenvalues[0] * (A[0]*V[0] + A[1]*V[2]);
	U[1] = eigenvalues[0] * (A[2]*V[0] + A[3]*V[2]);
	U[2] = eigenvalues[1] * (A[0]*V[1] + A[1]*V[3]);
	U[3] = eigenvalues[1] * (A[2]*V[1] + A[3]*V[3]);

	AAT[0] = V[0]*eigenvalues[0];
	AAT[1] = V[1]*eigenvalues[1];
	AAT[2] = V[2]*eigenvalues[0];
	AAT[3] = V[3]*eigenvalues[1];

	V[0] = AAT[0]*U[0] + AAT[1]*U[2];
	V[1] = AAT[0]*U[1] + AAT[1]*U[3];
	V[2] = AAT[2]*U[0] + AAT[3]*U[2];
	V[3] = AAT[2]*U[1] + AAT[3]*U[3];

	p[0] = V[0]*q[0] + V[1]*q[1];
	p[1] = V[2]*q[0] + V[3]*q[1];
	
	//end svd

	int ne = 0;
	int nn = 0;

	size_t maxwell = 4294967295;

	for (unsigned i = 0; i < *nelist_set_size; i++)
	{
		if(nelist[vid*(*nelist_set_size)+i] != maxwell)
			ne++;
	}

	for (unsigned i = 0; i < *nnlist_set_size; i++)
	{
		if(nnlist[vid*(*nnlist_set_size)+i] != maxwell)
			nn++;
	}

	if(ne < nn){
		p[0] -= p[0]*fabs(normals[2*vid]);
		p[1] -= p[1]*fabs(normals[2*vid+1]);
	}

	float sdf = coords[vid];

	coords[2*vid] += p[0];
	coords[2*vid+1] += p[1];

	float new_worst_q = 1.0;

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

		new_worst_q = min(new_worst_q, quality);
	}

	if(new_worst_q < worst_q)
	{
		coords[2*vid] -= p[0];
		coords[2*vid+1] -= p[1];
		return;
	}









}