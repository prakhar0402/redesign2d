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


#include "vertexmemo.h"

VertexMemo::VertexMemo(Kernel::Point_2 pnt, size_t location)
{
	v = pnt;
	loc = location;
	normal.set_size(2);
	normal.fill(0.0);
	velocity.set_size(2);
	velocity.fill(0.0);
	sdf = 0.0;
	sdf_force = 0.0;
	corner_force = 0.0;
	ref_point = Kernel::Point_2(0.0, 0.0);
	initial_length = compute_length();
	current_length = 0.0;
	force.set_size(2);
	Jpos.set_size(2, 2);
	Jvel.set_size(2, 2);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
	isMovable = false;
}

void VertexMemo::set_vertex(Kernel::Point_2 pnt)
{
	v = pnt;
}

void VertexMemo::set_vertex(double x, double y)
{
	v = Kernel::Point_2(x, y);
}

void VertexMemo::set_normal(Kernel::Vector_3 vertex_normal)
{
	normal << CGAL::to_double(vertex_normal.x()) << CGAL::to_double(vertex_normal.y());
}

void VertexMemo::set_normal(double x, double y)
{
	normal << x << y;
}

void VertexMemo::set_velocity(const arma::vec& vel)
{
	velocity = vel;
}

void VertexMemo::set_velocity(double x, double y)
{
	velocity << x << y;
}

void VertexMemo::set_sdf(double dia)
{
	sdf = dia;
}

void VertexMemo::set_ref_point(Kernel::Point_3 ref)
{
	ref_point = Kernel::Point_2(ref.x(), ref.y());
	initial_length = compute_length();
}

void VertexMemo::set_ref_point(double x, double y)
{
	ref_point = Kernel::Point_2(x, y);
	initial_length = compute_length();
}

double VertexMemo::compute_length()
{
	current_length = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(v, ref_point)));
	return current_length;
}

//void VertexMemo::compute_sdf_force(const double& K_sdf, const double& threshold_dia)
//{
//	sdf_force = 0.0;
//	if (threshold_dia > sdf)
//		sdf_force = K_sdf*(1.0 - sdf / threshold_dia);
//		//sdf_force = -K_sdf*area*std::log(sdf / threshold_dia);
//		//sdf_force = K_sdf*area*(threshold_dia - sdf);
//		//sdf_force = K_sdf*area*CGAL::square(threshold_dia - sdf);
//}

void VertexMemo::compute_force(const double& K_sdf, const double& K_s, const double& K_d)
{
	if (isMovable)
	{
		compute_length();
		//compute_sdf_force(K_sdf, threshold_dia);

		arma::vec pos1, pos2, vec12;
		pos1 << CGAL::to_double(v.x()) << CGAL::to_double(v.y());
		pos2 << CGAL::to_double(ref_point.x()) << CGAL::to_double(ref_point.y());
		vec12 = pos2 - pos1; // vector pointing from pos1 to pos2

		arma::mat mv12 = vec12 * arma::trans(vec12); // product of different elements of vec12

		double coef1 = 2.0*K_s*(1 - initial_length / current_length);
		double coef2 = 2.0*K_s*initial_length / pow(current_length, 3);

		force = coef1*vec12 + K_d*(-velocity) + (sdf_force + corner_force)*normal;

		Jpos = -coef2*mv12;
		Jpos.diag() -= coef1;

		Jvel.diag() -= K_d;
	}
}

double VertexMemo::get_sdf()
{
	return sdf;
}

Kernel::Point_2 VertexMemo::get_vertex()
{
	return v;
}

arma::vec VertexMemo::get_normal()
{
	return normal;
}

arma::vec VertexMemo::get_velocity()
{
	return velocity;
}

arma::vec VertexMemo::get_sdf_force()
{
	return sdf_force*normal;
}

arma::vec VertexMemo::get_force()
{
	return force;
}

arma::mat VertexMemo::get_Jacobian_pos()
{
	return Jpos;
}

arma::mat VertexMemo::get_Jacobian_vel()
{
	return Jvel;
}

void VertexMemo::setSDFForceMag(double f)
{
	sdf_force = f;
}

void VertexMemo::setMaxSDFMag(double f)
{
	sdf_force = std::max(sdf_force, f);
}

void VertexMemo::setCornerForceMag(double f)
{
	corner_force = f;
	if (corner_force != 0.0) isMovable = true;
}