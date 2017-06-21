#include "edgememo.h"

EdgeMemo::EdgeMemo(VertexMemo *v1, VertexMemo *v2)
{
	source = v1;
	target = v2;
	initial_length = compute_length();
	current_length = initial_length;
	force.set_size(2);
	Jpos.set_size(2, 2);
	Jvel.set_size(2, 2);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

double EdgeMemo::compute_length()
{
	current_length = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(source->get_vertex(), target->get_vertex())));
	return current_length;
}

double EdgeMemo::get_current_length()
{
	return current_length;
}

void EdgeMemo::compute_force(const double& K_s, const double& K_d)
{
	compute_length();

	arma::vec pos1, pos2, vec12;
	pos1 << CGAL::to_double(source->get_vertex().x()) << CGAL::to_double(source->get_vertex().y());
	pos2 << CGAL::to_double(target->get_vertex().x()) << CGAL::to_double(target->get_vertex().y());
	vec12 = pos2 - pos1; // vector pointing from pos1 to pos2

	arma::mat mv12 = vec12 * arma::trans(vec12); // product of different elements of vec12
	
	double coef1 = K_s*(1 - initial_length / current_length);
	double coef2 = K_s*initial_length / pow(current_length, 3);

	force = coef1*vec12 + K_d*(target->get_velocity() - source->get_velocity());

	Jpos = -coef2*mv12;
	Jpos.diag() -= coef1;

	Jvel.diag() -= K_d;
}

arma::vec EdgeMemo::get_force()
{
	return force;
}

arma::mat EdgeMemo::get_Jacobian_pos()
{
	return Jpos;
}

arma::mat EdgeMemo::get_Jacobian_vel()
{
	return Jvel;
}