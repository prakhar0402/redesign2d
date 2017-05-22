#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/connect_holes.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/enum.h>
#include <CGAL/boost/graph/graph_traits_Arrangement_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2							Point;
typedef CGAL::Polygon_2<Kernel>					Polygon;
typedef CGAL::Polygon_with_holes_2<Kernel>		Pwh;
typedef std::list<Pwh>							Pwh_list;

typedef Kernel::Vector_2 Vector;
typedef std::map<Polygon::Edge_const_iterator, Vector> Edge_vector_map;
typedef std::map<Polygon::Vertex_const_iterator, Vector> Vertex_vector_map;
typedef std::map<Polygon::Edge_const_iterator, double> Edge_double_map;


// FILE FORMAT:
// #number of polygons with holes <num_poly>
// #number of holes in polygon 1 <num_holes>
// #number of vertices in outer boundary of polygon 1 <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// #number of vertices in hole 1 of polygon 1 <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// #number of vertices in hole 2 of polygon 1 <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// ...
// ...
// #number of vertices in hole num_holes of polygon 1 <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// ...
// ...
// #number of holes in polygon num_poly <num_holes>
// #number of vertices in outer boundary of polygon num_poly <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// #number of vertices in hole 1 of polygon num_poly <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// #number of vertices in hole 2 of polygon num_poly <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>
// ...
// ...
// #number of vertices in hole num_holes of polygon num_poly <num_V>
// <x1> <y1>
// <x2> <y2>
// ...
// <x_num_V> <y_num_V>


// writing polygon to a string
std::string write(Polygon pg, bool isInPwh = false)
{
	std::ostringstream oss;
	if (!isInPwh)
		oss << "1\n0\n";
	oss << pg.size() << "\n";
	Polygon::Vertex_const_iterator it;
	for (it = pg.vertices_begin(); it != pg.vertices_end(); it++)
		oss << it->x() << " " << it->y() << "\n";
	return oss.str();
}

// writing polygon with holes to a string
std::string write(Pwh pwh, bool isInList = false)
{
	std::ostringstream oss;
	if (!isInList)
		oss << "1\n";
	oss << pwh.number_of_holes() << "\n";
	oss << write(pwh.outer_boundary(), true);
	Pwh::Hole_const_iterator hit;
	for (hit = pwh.holes_begin(); hit != pwh.holes_end(); hit++)
		oss << write(*hit, true);
	return oss.str();
}

// writing polygon with holes list to a string
std::string write(Pwh_list pwh_list)
{
	std::ostringstream oss;
	oss << pwh_list.size() << "\n";
	Pwh_list::const_iterator  it;
	for (it = pwh_list.begin(); it != pwh_list.end(); it++)
		oss << write(*it, true);
	return oss.str();
}

// reading polygon with holes from a file
bool read(std::string filename, Pwh_list& result)
{
	std::ifstream file(filename);

	if (!file.is_open())
	{
		std::cerr << "Failed to open the input file." << std::endl;
		return false;
	}

	size_t num_poly, num_holes, num_V;
	file >> num_poly;

	double x, y;
	for (size_t i = 0; i < num_poly; i++)
	{
		Pwh pwh;
		file >> num_holes;
		file >> num_V;
		for (size_t j = 0; j < num_V; j++)
		{
			file >> x;
			file >> y;
			Point p(x, y);
			pwh.outer_boundary().push_back(p);
		}
		for (size_t k = 0; k < num_holes; k++)
		{
			Polygon pg;
			file >> num_V;
			for (size_t j = 0; j < num_V; j++)
			{
				file >> x;
				file >> y;
				Point p(x, y);
				pg.push_back(p);
			}
			pwh.add_hole(pg);
		}
		result.push_back(pwh);
	}
	file.close();
	return true;
}

// get vertex normals of a simple polygon
void getNormals(const Polygon& pg, std::ofstream& ofile, bool isHole = false)
{
	Edge_double_map edm;
	boost::associative_property_map<Edge_double_map> Edge_length_map(edm);

	int n = 0;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ei++)
		boost::put(Edge_length_map, ei, CGAL::sqrt(ei->to_vector().squared_length()));

	Vertex_vector_map vvm;
	boost::associative_property_map<Vertex_vector_map> Vertex_normal_map(vvm);

	Vector normal;
	Polygon::Vertex_iterator vi = pg.vertices_begin();
	Polygon::Edge_const_iterator ei = pg.edges_begin();

	normal = -((pg.edges_end() - 1)->to_vector() + ei->to_vector()).perpendicular(pg.orientation());
	if (isHole)
		normal = -normal;
	normal = normal / CGAL::sqrt(normal.squared_length());
	boost::put(Vertex_normal_map, vi, normal);
	ofile << normal << std::endl;
	vi++; ei++;

	for (; vi != pg.vertices_end(); vi++, ei++)
	{
		normal = -((ei - 1)->to_vector() + ei->to_vector()).perpendicular(pg.orientation());
		if (isHole)
			normal = -normal;
		normal = normal / CGAL::sqrt(normal.squared_length());
		boost::put(Vertex_normal_map, vi, normal);
		ofile << normal << std::endl;
	}
}


// get vertex normals of a polygon with holes
void getNormals(const Pwh& pwh, std::ofstream& ofile)
{
	getNormals(pwh.outer_boundary(), ofile);
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		getNormals(*hi, ofile, true);
}


int main(int argc, char *argv[])
{
	std::string FILENAME = "data/example.dat";
	Pwh_list slice;
	read(FILENAME, slice);

	Pwh pwh = slice.front();

	std::ofstream ofile("data/normal.dat");
	getNormals(pwh, ofile);
	ofile.close();

	return 0;
}