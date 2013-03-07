#pragma OPENCL EXTENSION cl_khr_fp64 : enable

bool isCornerNode(size_t vid)
{
	return fabs(normals[2*vid])==1.0 && fabs(normals[2*vid+1]==1.0);
}

bool isSurfaceNode(size_t vid)
{
  return NEList_size < NNList_size; // this could be really wrong!!!!!!!
  // maybe should be NE_List at offset.size < NNList at offset.size
}

double element_area(size_t eid)
{
	const size_t *n = &ENlist[3*eid];

	// Pointers to the coordinates of each vertex
	const double *c0 = &coords[2*n[0]];
	const double *c1 = &coords[2*n[1]];
	const double *c2 = &coords[2*n[2]];

	return (orientation * 0.5 * ((c0[1] - c2[1]) * (c0[0] - c1[0]) - (c0[1] - c1[1]) * (c0[0] - c2[0])));
}

double element_quality(size_t eid)
{
	const size_t *n = &ENList[3*eid];

	// Pointers to the coordinates of each vertex
	const double *c0 = &coords[2*n[0]];
	const double *c1 = &coords[2*n[1]];
	const double *c2 = &coords[2*n[2]];

	// Pointers to the metric tensor at each vertex
	const double *m0 = &metric[3*n[0]];
	const double *m1 = &metric[3*n[1]];
	const double *m2 = &metric[3*n[2]];

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
	double a = element_area(eid);

	// Area in metric space
	double a_m = a*sqrt(m00*m11 - m01*m01);

	// Function
	double f = min(l/3.0, 3.0/l);
	double F = pow(f * (2.0 - f), 3.0);

	// This is the 2D Lipnikov functional.
	double quality = 12.0 * sqrt(3.0) * a_m * F / (l*l);

	return quality;
}

void calc_eigenvalues(const double A[4], double eigenvalues[2])
{
	double b, discriminant;

	b = A[0]+A[3];

	discriminant = sqrt(b*b - 4*(A[0]*A[3] - A[1]*A[2]));

	eigenvalues[0] = (b + discriminant) * 0.5;
	eigenvalues[1] = (b - discriminant) * 0.5;
}


void calc_eigenvectors(const double A[4], const double eigenvalues[2], double eigenvectors[4])
{
	double D[4];
	double proportion;

	D[1] = A[1];
	D[2] = A[2];

	D[0] = A[0] - eigenvalues[0];
	D[3] = A[3] - eigenvalues[0];

	if(D[1] != 0.0)
	{
		proportion = -(D[0]/D[1]);
		eigenvectors[0] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[2] = proportion * eigenvectors[0];
	}
	else if(D[3] != 0.0)
	{
		proportion = -(D[2]/D[3]);
		eigenvectors[0] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[2] = proportion * eigenvectors[0];
	}
	else if(D[0] == 0.0)
	{
		eigenvectors[0] = 1.0;
		eigenvectors[2] = 0.0;
	}
	else if(D[2] == 0.0)
	{
		eigenvectors[2] = 1.0;
		eigenvectors[0] = 0.0;
	}

	D[0] = A[0] - eigenvalues[1];
	D[3] = A[3] - eigenvalues[1];

	if(D[1] != 0.0)
	{
		proportion = -(D[0]/D[1]);
		eigenvectors[1] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[3] = proportion * eigenvectors[1];
	}
	else if(D[3] != 0.0)
	{
		proportion = -(D[2]/D[3]);
		eigenvectors[1] = sqrt(1.0 / (1 + proportion*proportion));
		eigenvectors[3] = proportion * eigenvectors[1];
	}
	else if(D[0] == 0.0)
	{
		eigenvectors[1] = 1.0;
		eigenvectors[3] = 0.0;
	}
	else if(D[2] == 0.0)
	{
		eigenvectors[3] = 1.0;
		eigenvectors[1] = 0.0;
	}
}

void svd_solve_2x2(const double A[4], double p[2], const double q[2]){

	double AAT[4]; // This will be used to store either A*Atransp or Atransp*A
	double eigenvalues[2];
	double U[4];
	double V[4];

	// Caclulate Atransp*A
	AAT[0] = A[0]*A[0] + A[2]*A[2];
	AAT[1] = A[0]*A[1] + A[2]*A[3];
	AAT[2] = AAT[1];
	AAT[3] = A[1]*A[1] + A[3]*A[3];

	// Calculate the eigenvalues of AT*A
	calc_eigenvalues(AAT, eigenvalues);

	// Calculate the right singular vector V:
	calc_eigenvectors(AAT, eigenvalues, V);

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
}



__kernel void nodes_solve(
		__global size_t *NEList,
		__global size_t *NEList_size
		__global size_t *NEList_Offset,
		__global size_t *NEList_Offset_size,
		__global size_t *ENList,
		__global size_t *ENList_size
		__global double *coords,
		__global double *coords_size,
		__global double *metric,
		__global double *metric_size,
		__global int *orientation,
		__global double *normals,
		__global double *normals_size,
		__global size_t *NNList_size
		) 
{

    // Get the index of the current element to be processed
    int vid = get_global_id(0);

	// calculate number of adjacent nodes count = (end-start)/sizepernode
	int adjacentNodeCount = abs(NEList_Offset[vid+1]-NEList_Offset[vid])/sizeof(size_t);

	if(isCornerNode(vid))
		return; // was continue, but illegal not in loop (how to end thread?)

	double worst_q = 1.0;

	//psss pointer to nelist for each iteration as a size_t ptr NEList_vid
	for(size_t i = 0; i < NEList_size; i++)
	{
		worst_q = min(worst_q, element_quality(i));
	}

	const double *m0 = metric[3*vid];
	double x0 = coords[2*vid];
	double y0 = coords[2*vid+1];

	double A[4] = {0.0, 0.0, 0.0, 0.0};
	double q[2] = {0.0, 0.0};

	// Iterate over all edges and assemble matrices A and q
	for(size_t j = 0; j < NEList_size; j++)
	{
		size_t il = j;

		const double *m1 = &metric[3*il];

		// Find the metric in the middle of the edge.
		double ml00 = 0.5*(m0[0] + m1[0]);
		double ml01 = 0.5*(m0[1] + m1[1]);
		double ml11 = 0.5*(m0[2] + m1[2]);

		double x = coords[2*il] - x0;
		double y = coords[2*il+1] - y0;

		// Calculate and accumulate the contribution of
		// this vertex to the barycentre of the cavity.
		q[0] += (ml00*x + ml01*y);
		q[1] += (ml01*x + ml11*y);

		A[0] += ml00;
		A[1] += ml01;
		A[3] += ml11;
	}

	A[2] = A[1];

	double p[2];

	svd_solve_2x2(A, p, q);


	if(isSurfaceNode(vid))
	{
		p[0] -= p[0]*fabs(normals[2*vid]);
		p[1] -= p[1]*fabs(normals[2*vid+1]);
	}


	coords[2*vid] += p[0];
	coords[2*vid+1] += p[1];

	double new_worst_q = 1.0;
	for(size_t i = 0; i < NEList_size; i++)
	{
		new_worst_q = min(new_worst_q, element_quality(i));
	}

	if(new_worst_q < worst_q)
	{
		coords[2*vid] -= p[0];
		coords[2*vid+1] -= p[1];
	}


	// now the coords which have been updated should be made available in the program.
}