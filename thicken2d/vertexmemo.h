/*
MAD Lab, University at Buffalo
Copyright (C) 2018  Prakhar Jaiswal <prakharj@buffalo.edu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


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
	VertexMemo(Kernel::Point_2 pnt, size_t location);
	void set_vertex(Kernel::Point_2 pnt);
	void set_vertex(double x, double y);
	void set_normal(Kernel::Vector_3 vertex_normal);
	void set_normal(double x, double y);
	void set_velocity(const arma::vec& vel);
	void set_velocity(double x, double y);
	void set_ref_point(Kernel::Point_3 ref);
	void set_ref_point(double x, double y);
	void set_sdf(double dia);
	double compute_length();
	//void compute_sdf_force(const double& K_sdf, const double& threshold_dia);
	double get_sdf();
	Kernel::Point_2 get_vertex();
	arma::vec get_normal();
	arma::vec get_velocity();
	arma::vec get_sdf_force();
	size_t loc; // location in vector of all vertices
	size_t index; // global matrix index for movable vertice
	bool isMovable;
	void compute_force(const double& K_sdf, const double& K_s, const double& K_d);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();
	void setSDFForceMag(double f);
	void setMaxSDFMag(double f);
	void setCornerForceMag(double f);

private:
	Kernel::Point_2 v;
	arma::vec normal; // vertex normal
	arma::vec velocity;
	Kernel::Point_2 ref_point;
	double initial_length;
	double current_length;
	double sdf;
	double sdf_force;
	double corner_force;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif