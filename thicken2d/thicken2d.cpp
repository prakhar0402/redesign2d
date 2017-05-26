#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <random>
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
#include <CGAL/number_utils.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2							Point_2;
typedef CGAL::Polygon_2<Kernel>					Polygon;
typedef CGAL::Polygon_with_holes_2<Kernel>		Pwh;
typedef std::list<Pwh>							Pwh_list;

typedef Kernel::Vector_2 Vector_2;
typedef std::map<Polygon::Edge_const_iterator, Vector_2> Edge_vector_map;
typedef std::map<Polygon::Vertex_const_iterator, Vector_2> Vertex_vector_map;
typedef std::map<Polygon::Edge_const_iterator, double> Edge_double_map;

typedef Kernel::Segment_3 Segment;
typedef Kernel::Ray_3 Ray;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::list<Segment> Seg_list;
typedef Seg_list::iterator Iterator;
typedef CGAL::AABB_segment_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_segment_traits;
typedef CGAL::AABB_tree<AABB_segment_traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> Segment_intersection;

// Constants
size_t num_samples = 25;
double cone_half_angle = CGAL_PI / 3.0;
std::uniform_real_distribution<double> unif(-cone_half_angle, cone_half_angle);
std::default_random_engine re;


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
			Point_2 p(x, y);
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
				Point_2 p(x, y);
				pg.push_back(p);
			}
			pwh.add_hole(pg);
		}
		result.push_back(pwh);
	}
	file.close();
	return true;
}


// get vertices as vector of Point_3 for a simple polygon
void getVertexVector(const Polygon& pg, std::vector<Point_3>& vertex_vec)
{
	for (Polygon::Vertex_iterator vi = pg.vertices_begin(); vi != pg.vertices_end(); vi++)
		vertex_vec.push_back(Point_3(vi->x(), vi->y(), 0.0));
}

// get vertices as vector of Point_3 for a polygon with holes
void getVertexVector(const Pwh& pwh, std::vector<Point_3>& vertex_vec)
{
	getVertexVector(pwh.outer_boundary(), vertex_vec);
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		getVertexVector(*hi, vertex_vec);
}



// get vertex normals of a simple polygon
void getNormals(const Polygon& pg, std::vector<Vector_3>& normal_vec, bool isHole = false)
{
	Edge_double_map edm;
	boost::associative_property_map<Edge_double_map> Edge_length_map(edm);

	int n = 0;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ei++)
		boost::put(Edge_length_map, ei, CGAL::sqrt(CGAL::to_double(ei->to_vector().squared_length())));

	Vertex_vector_map vvm;
	boost::associative_property_map<Vertex_vector_map> Vertex_normal_map(vvm);

	Vector_2 normal;
	Polygon::Vertex_iterator vi = pg.vertices_begin();
	Polygon::Edge_const_iterator ei = pg.edges_begin();

	normal = -((pg.edges_end() - 1)->to_vector() + ei->to_vector()).perpendicular(pg.orientation());
	if (isHole)
		normal = -normal;
	normal = normal / CGAL::sqrt(CGAL::to_double(normal.squared_length()));
	normal_vec.push_back(Vector_3(normal.x(), normal.y(), 0.0));
	boost::put(Vertex_normal_map, vi, normal);
	vi++; ei++;

	for (; vi != pg.vertices_end(); vi++, ei++)
	{
		normal = -((ei - 1)->to_vector() + ei->to_vector()).perpendicular(pg.orientation());
		if (isHole)
			normal = -normal;
		normal = normal / CGAL::sqrt(CGAL::to_double(normal.squared_length()));
		normal_vec.push_back(Vector_3(normal.x(), normal.y(), 0.0));
		boost::put(Vertex_normal_map, vi, normal);
	}
}


// get vertex normals of a polygon with holes
void getNormals(const Pwh& pwh, std::vector<Vector_3>& normal_vec)
{
	getNormals(pwh.outer_boundary(), normal_vec);
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		getNormals(*hi, normal_vec, true);
}


// sample rays in a cone vertexed at source with specified axis direction
// TODO: fix the random angles for all cones to be same?
void sampleRays(const Point_3& source, const Vector_3& axis, std::vector<Ray>& rays)
{
	double angle, cos, sin;
	for (size_t i = 0; i < num_samples; i++)
	{
		angle = unif(re);
		cos = std::cos(angle);
		sin = std::sin(angle);
		rays.push_back(Ray(source, Vector_3(axis.x()*cos - axis.y()*sin, axis.x()*sin + axis.y()*cos, 0.0)));
	}
}

// convert polygon with holes into list of 3D edge segments
void getEdgeSegmentList(const Pwh& pwh, Seg_list& seg_list)
{
	for (Polygon::Edge_const_iterator ei = pwh.outer_boundary().edges_begin(); ei != pwh.outer_boundary().edges_end(); ei++)
		seg_list.push_back(Segment(Point_3(ei->source().x(), ei->source().y(), 0.0), Point_3(ei->target().x(), ei->target().y(), 0.0)));
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		for (Polygon::Edge_const_iterator ei = hi->edges_begin(); ei != hi->edges_end(); ei++)
			seg_list.push_back(Segment(Point_3(ei->source().x(), ei->source().y(), 0.0), Point_3(ei->target().x(), ei->target().y(), 0.0)));

	//Tree tree(seg_list.begin(), seg_list.end());
}

// get minimum distance of all intersections
// TODO: igonre outside intersections IMPORTANT
double closestIntersectionDistance(const Point_3& source, const std::list<Ray_intersection>& intersections)
{
	double dist = 0.0, min = 0.0;
	for (std::list<Ray_intersection>::const_iterator it = intersections.begin(); it != intersections.end(); it++)
	{
		const Point_3* p = boost::get<Point_3>(&(it->value().first));
		if (p)
		{
			dist = CGAL::to_double(CGAL::squared_distance(*p, source));
			if (min == 0.0) min = dist;
			if (dist != 0.0 && min > dist) // TODO: use a small number like 1.0e-9 instead of 0.0?
				min = dist;
		}
	}
	return CGAL::sqrt(min);
}

// remove outliers and get average
double niceAverage(std::vector<double> numbers)
{
	double sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
	double mean = sum / numbers.size();

	std::vector<double> diff(numbers.size());
	std::transform(numbers.begin(), numbers.end(), diff.begin(), [mean](double x) { return x - mean; });
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / numbers.size());

	size_t n = numbers.size() / 2;
	std::nth_element(numbers.begin(), numbers.begin() + n, numbers.end());
	double median = numbers.at(n);

	double lower = median - stdev;
	double upper = median + stdev;

	size_t count = 0;
	sum = 0.0;
	for (std::vector<double>::iterator it = numbers.begin(); it != numbers.end(); it++)
	{
		if (*it >= lower && *it <= upper)
		{
			count++;
			sum += *it;
		}
	}
	return (count == 0) ? 0.0 : sum / count;
}

// get SDF values
void getSDF(const Pwh& pwh, std::vector<double>& sdf)
{
	Seg_list seg_list;
	getEdgeSegmentList(pwh, seg_list);

	Tree tree(seg_list.begin(), seg_list.end());

	std::vector<Point_3> vertex_vec;
	std::vector<Vector_3> normal_vec;

	getVertexVector(pwh, vertex_vec);
	getNormals(pwh, normal_vec);

	std::vector<Ray> rays;
	std::list<Ray_intersection> intersections;

	std::vector<double> dists;
	for (size_t i = 0; i < vertex_vec.size(); i++)
	{
		rays.clear(); dists.clear();
		sampleRays(vertex_vec.at(i), -normal_vec.at(i), rays);

		for (std::vector<Ray>::iterator ri = rays.begin(); ri != rays.end(); ri++)
		{
			intersections.clear();
			tree.all_intersections(*ri, std::back_inserter(intersections));

			dists.push_back(closestIntersectionDistance(vertex_vec.at(i), intersections));
		}

		sdf.push_back(niceAverage(dists));
	}
}

// AABB intersection test
void intersect(Polygon& pg, Point_2& source, Vector_2& dir)
{
	Seg_list seg_list;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ei++)
	{
		seg_list.push_back(Segment(Point_3(ei->source().x(), ei->source().y(), 0.0), Point_3(ei->target().x(), ei->target().y(), 0.0)));
	}

	Tree tree(seg_list.begin(), seg_list.end());

	Ray ray(Point_3(source.x(), source.y(), 0.0), Vector_3(dir.x(), dir.y(), 0.0));

	std::list<Ray_intersection> intersections;
	tree.all_intersections(ray, std::back_inserter(intersections));

	for (std::list<Ray_intersection>::const_iterator it = intersections.begin(); it != intersections.end(); it++)
	{
		const Point_3* p = boost::get<Point_3>(&(it->value().first));
		if (p)
			std::cout << "intersection object is a point " << *p << std::endl;
	}
}




int main(int argc, char *argv[])
{
	re.seed(0);

	std::string FILENAME = "data/example.dat";
	Pwh_list slice;
	read(FILENAME, slice);

	Pwh pwh = slice.front();

	//std::ofstream ofile("data/normal.dat");
	//getNormals(pwh, ofile);
	//ofile.close();

	//intersect(pwh.outer_boundary(), Point_2(-3.75, -2.375), Vector_2(-1.0, 0.0));

	std::vector<double> sdf;

	getSDF(pwh, sdf);

	std::ofstream ofile("data/SDF.dat");
	std::ostream_iterator<double> output_iterator(ofile, "\n");
	std::copy(sdf.begin(), sdf.end(), output_iterator);
	ofile.close();

	return 0;
}