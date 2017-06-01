#ifndef EDGEMEMO_H
#define EDGEMEMO_H

#include "vertexmemo.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>

#include <armadillo>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon;

class EdgeMemo
{
public:
	EdgeMemo(VertexMemo *v1, VertexMemo *v2);
	VertexMemo *source, *target;
	double compute_length();
	double get_current_length();
	void compute_force(const double& K_s, const double& K_d);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();

private:
	double initial_length;
	double current_length;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif