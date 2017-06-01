#ifndef VERTEXMEMO_H
#define VERTEXMEMO_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>

#include <armadillo>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon;

class VertexMemo
{
public:
	VertexMemo(Kernel::Point_2 pnt, size_t location, size_t id);
	void set_vertex(Kernel::Point_2 pnt);
	void set_normal(Kernel::Vector_3 vertex_normal);
	void set_velocity(const arma::vec& vel);
	void set_ref_point(Kernel::Point_3 ref);
	void set_sdf(double dia);
	void set_area_factor(double value);
	double compute_length();
	void compute_sdf_force(const double& K_sdf, const double& threshold_dia);
	double get_sdf();
	double get_area_factor();
	Kernel::Point_2 get_vertex();
	arma::vec get_normal();
	arma::vec get_velocity();
	arma::vec get_sdf_force();
	size_t pgId; // 0 for outer_boundary, and 1,2,3,... for holes 0,1,2,...
	size_t loc; // location in vector of all vertices
	size_t index; // global matrix index for movable vertice
	bool isMovable;
	void compute_force(const double& K_sdf, const double& K_s, const double& K_d, const double& threshold_dia);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();

private:
	Kernel::Point_2 v;
	arma::vec normal; // vertex normal
	arma::vec velocity;
	Kernel::Point_2 ref_point;
	double initial_length;
	double current_length;
	double sdf;
	double area_factor; // sum of projected area of all faces around the vertex divided by average projected area for all vertices
	arma::vec sdf_force;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif