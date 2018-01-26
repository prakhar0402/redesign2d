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