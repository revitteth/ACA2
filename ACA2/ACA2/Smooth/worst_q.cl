__kernel void vector_add(
	__global size_t *vid,
	__global vector<size_t> *ENList, 
	__global vector< set<size_t> > *NEList,
	__global vector<double> *metric,
	__global vector<double> *coords,
	__global int *orientation,
	__global double *worst_q) 
{
 
	set<size_t>::const_iterator it = get_global_id(vid);

	worst_q = min(worst_q, element_quality(
		it, vid, ENList, NEList, metric, coords, orientation
		));
}

double element_quality(
	size_t eid, 
	size_t vid, 
	std::vector<size_t> ENList, 
	std::vector< std::set<size_t> > NEList, 
	std::vector<double> metric, 
	std::vector<double> coords, 
	int orientation) const
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
	double f = std::min(l/3.0, 3.0/l);
	double F = pow(f * (2.0 - f), 3.0);

	// This is the 2D Lipnikov functional.
	double quality = 12.0 * sqrt(3.0) * a_m * F / (l*l);

	return quality;
}