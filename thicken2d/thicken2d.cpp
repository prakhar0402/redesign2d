#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <random>
#include <algorithm>
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

#include "vertexmemo.h"
#include "edgememo.h"

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
std::string FILENAME = "data/example.dat";

size_t num_samples = 10;
double cone_half_angle = CGAL_PI / 12.0;
std::uniform_real_distribution<double> unif(-cone_half_angle, cone_half_angle);
std::default_random_engine re;

double T1 = 0.4;

double K_SDF = 64.0;
double K_S = 32.0;
double K_D = 16.0;

double VERTEX_MASS = 1.0;
double MASS_INV = 1.0 / VERTEX_MASS;

double TIME_STEP = 1.0;
double TOTAL_TIME = 3.0;
size_t STEPS = (size_t)std::ceil(TOTAL_TIME / TIME_STEP);

// Global Variables
std::vector<Point_3> vertex_vec;
std::vector<Vector_3> normal_vec;
std::vector<std::pair<double, double>> angles;
std::vector<double> sdf_vec;
std::vector<VertexMemo> vmemo_vec;
std::vector<EdgeMemo> ememo_vec;



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

// write normals to a file
void writeNormals(std::string filename)
{
	std::ofstream ofile(filename);
	for (std::vector<Vector_3>::iterator it = normal_vec.begin(); it != normal_vec.end(); it++)
		ofile << *it << "\n";
	ofile.close();
}

// write angles to a file
void writeAngles(std::string filename)
{
	std::ofstream ofile(filename);
	for (std::vector<std::pair<double, double>>::iterator it = angles.begin(); it != angles.end(); it++)
		ofile << it->first << " " << it->second << "\n";
	ofile.close();
}

// get vertices as vector of Point_3 for a simple polygon
void getVertexVector(const Polygon& pg)
{
	for (Polygon::Vertex_iterator vi = pg.vertices_begin(); vi != pg.vertices_end(); vi++)
		vertex_vec.push_back(Point_3(vi->x(), vi->y(), 0.0));
}

// get vertices as vector of Point_3 for a polygon with holes
void getVertexVector(const Pwh& pwh)
{
	vertex_vec.clear();
	getVertexVector(pwh.outer_boundary());
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		getVertexVector(*hi);
}

// get counter-clockwise angle between two vectors
double ccwAngle(const Vector_2& v1, const Vector_2& v2)
{
	double dot = CGAL::to_double(v1.x()*v2.x() + v1.y()*v2.y());
	double det = CGAL::to_double(v1.x()*v2.y() - v1.y()*v2.x());
	double angle = std::atan2(det, dot);
	return (angle >= 0) ? angle : 2.0 * CGAL_PI + angle;
}

// computes the perimeter of a polygon
double getPerimeter(const Polygon& pg)
{
	double dist = 0.0;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ei++)
		dist += CGAL::sqrt(CGAL::to_double(ei->squared_length()));
	return dist;
}

// TODO: check and correct orientation of holes


// resamples a polygon with equidistant points
Polygon resample(const Polygon& pg, double delta)
{
	// initialize result holder
	Polygon rpg;

	double dist; int N = 0;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ei++)
	{
		Vector_2 edge_vec = ei->to_vector();
		dist = CGAL::sqrt(CGAL::to_double(edge_vec.squared_length()));
		N = (int)std::ceil(dist / delta);
		edge_vec = edge_vec / N;
		for (int i = 0; i < N; i++)
		{
			Point_2 pnt(ei->source().x() + edge_vec.x()*i, ei->source().y() + edge_vec.y()*i);
			rpg.push_back(pnt);
		}
	}
	return rpg;


	
	//// initialize some variables
	//Polygon::Vertex_const_iterator vi = pg.vertices_begin();
	//Point_2 current(vi->x(), vi->y());
	//vi++;
	//Point_2 next(vi->x(), vi->y());
	//rpg.push_back(current);
	//double dist = 0.0, edgelength = 0.0;
	//bool loopEnd = false;

	//// loop until reaching end of polygon
	//while (!loopEnd || vi != pg.vertices_begin()+1)
	//{
	//	edgelength = CGAL::sqrt(CGAL::to_double((next - current).squared_length())); // get the current edge length
	//	if (edgelength + dist < delta)
	//	{
	//		// next vertex too close, increment distance and jump to next edge
	//		dist += edgelength;
	//		vi++;
	//		current = next;
	//		if (vi == pg.vertices_end())
	//		{
	//			loopEnd = true;
	//			vi = pg.vertices_begin();
	//		}
	//		next = Point_2(vi->x(), vi->y());

	//		std::cout << "-:-" << std::endl;
	//	}
	//	else
	//	{
	//		// resample a new point at a distance (delta - distance covered) from current on the edge <current, next>
	//		Point_2 pnt(current.x() + (delta - dist)*(next.x() - current.x()) / edgelength, current.y() + (delta - dist)*(next.y() - current.y()) / edgelength);

	//		// insert the new point in the result and at the next position and increment size of loop
	//		rpg.push_back(pnt);
	//		current = pnt;
	//		dist = 0.0; // reinitialize distance
	//	}
	//}
	//return rpg;
}

// resamples a polygon with equidistant points
Pwh resample(const Pwh& pwh, double delta)
{
	Pwh rpwh;
	rpwh.outer_boundary() = resample(pwh.outer_boundary(), delta);
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		rpwh.add_hole(resample(*hi, delta));
	return rpwh;
}

// get vertex normals of a simple polygon
void getNormalsAndAngles(const Polygon& pg, bool isHole = false)
{
	Vector_2 normal, edge1, edge2;
	Polygon::Vertex_iterator vi = pg.vertices_begin();
	Polygon::Edge_const_iterator ei = pg.edges_begin();

	edge1 = (pg.edges_end() - 1)->to_vector();
	edge2 = ei->to_vector();
	normal = (edge1 + edge2).perpendicular(pg.orientation());
	if (!isHole)
		normal = -normal;
	normal = normal / CGAL::sqrt(CGAL::to_double(normal.squared_length()));
	normal_vec.push_back(Vector_3(normal.x(), normal.y(), 0.0));
	angles.push_back(std::pair<double, double>(ccwAngle(-normal, -edge1), ccwAngle(edge2, -normal)));
	vi++; ei++;

	for (; vi != pg.vertices_end(); vi++, ei++)
	{
		edge1 = (ei - 1)->to_vector();
		edge2 = ei->to_vector();
		normal = (edge1 + edge2).perpendicular(pg.orientation());
		if (!isHole)
			normal = -normal;
		normal = normal / CGAL::sqrt(CGAL::to_double(normal.squared_length()));
		normal_vec.push_back(Vector_3(normal.x(), normal.y(), 0.0));
		angles.push_back(std::pair<double, double>(ccwAngle(-normal, -edge1), ccwAngle(edge2, -normal)));
	}
}

// get vertex normals of a polygon with holes
void getNormalsAndAngles(const Pwh& pwh)
{
	normal_vec.clear();
	angles.clear();
	getNormalsAndAngles(pwh.outer_boundary());
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
		getNormalsAndAngles(*hi, true);
}

// sample rays in a cone vertexed at source with specified axis direction
// TODO: fix the random angles for all cones to be same?
void sampleRays(const Point_3& source, const Vector_3& axis, std::pair<double, double> edge_angle, std::vector<Ray>& rays)
{
	double angle, cos, sin;
	for (size_t i = 0; i < num_samples; i++)
	{
		angle = unif(re);
		if (angle >= -edge_angle.first && angle <= edge_angle.second)
		{
			cos = std::cos(angle);
			sin = std::sin(angle);
			rays.push_back(Ray(source, Vector_3(axis.x()*cos - axis.y()*sin, axis.x()*sin + axis.y()*cos, 0.0)));
		}
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

// TODO: fix smoothing and then smooth sdf values
// smooth
std::vector<double> smooth(std::vector<double> numbers, size_t oneSideWindow)
{
	size_t windowSz = 2 * oneSideWindow + 1;
	std::vector<double> sNumbers;
	numbers.insert(numbers.end(), numbers.begin(), numbers.begin() + oneSideWindow);
	numbers.insert(numbers.begin(), numbers.end() - windowSz + 1, numbers.end() - oneSideWindow);
	double sum = std::accumulate(numbers.begin(), numbers.begin() + windowSz - 1, 0.0);
	for (std::vector<double>::iterator it = numbers.begin() + windowSz - 1; it != numbers.end(); it++)
	{
		sum += *it;
		sNumbers.push_back(sum / windowSz);
		sum -= *(it - windowSz + 1);
	}
	return sNumbers;
}

// get SDF values
void computeSDF(const Pwh& pwh, bool SDFexists = false)
{
	getVertexVector(pwh);
	getNormalsAndAngles(pwh);

	writeNormals("data/normals.dat");
	writeAngles("data/angles.dat");

	sdf_vec.clear();

	// use pre-existing sdf values if already computed before
	std::string sdf_file = FILENAME.substr(0, FILENAME.size() - 4) + ".sdf";
	std::ifstream ifile(sdf_file);
	if (SDFexists && ifile.good()) // if the file exists and can be read
	{
		// read sdf values from the file
		double val;
		while (ifile >> val)
			sdf_vec.push_back(val);
	}
	else // if file doesn't exist or can't be read, compute face sdf values and write to file
	{
		ifile.close();

		// compute SDF values
		Seg_list seg_list;
		getEdgeSegmentList(pwh, seg_list);

		Tree tree(seg_list.begin(), seg_list.end());

		std::vector<Ray> rays;
		std::list<Ray_intersection> intersections;

		std::vector<double> dists;
		for (size_t i = 0; i < vertex_vec.size(); i++)
		{
			rays.clear(); dists.clear();
			sampleRays(vertex_vec.at(i), -normal_vec.at(i), angles.at(i), rays);

			for (std::vector<Ray>::iterator ri = rays.begin(); ri != rays.end(); ri++)
			{
				intersections.clear();
				tree.all_intersections(*ri, std::back_inserter(intersections));

				dists.push_back(closestIntersectionDistance(vertex_vec.at(i), intersections));
			}

			sdf_vec.push_back(niceAverage(dists));
		}

		if (SDFexists) // TODO: remove this if statement
		{
			// write sdf values to a file
			std::ofstream ofile(sdf_file);
			std::ostream_iterator<double> output_iterator(ofile, "\n");
			std::copy(sdf_vec.begin(), sdf_vec.end(), output_iterator);
			ofile.close();
		}
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






void setMovable(
	const Polygon& pg,
	size_t startIdx,
	size_t curIdx,
	size_t n_ring = 0
)
{
	size_t idx = curIdx - n_ring;
	for (size_t i = 0; i <= 2*n_ring; i++, idx++)
	{
		if ((idx - startIdx) >= pg.size())
			idx = startIdx + (idx - startIdx) % pg.size();
		vmemo_vec[idx].isMovable = true;
	}
	//Vertex_memo_map[vi].isMovable = true;
	//if (n_ring <= 0) return;

	//Polyhedron::Halfedge_const_handle he1;
	//Polyhedron::Halfedge_const_handle he2;
	//he1 = vi->halfedge();
	//he2 = he1;
	//do
	//{
	//	setMovable(mesh, Vertex_memo_map, he2->opposite()->vertex(), n_ring - 1);
	//	he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
	//} while (he2 != he1);
}


void populate_memos(
	const Polygon& pg,
	size_t startIdx,
	size_t pgId
)
{
	size_t loc = startIdx;
	// loop over each vertex to comoute vertex sdf
	for (Polygon::Vertex_const_iterator vi = pg.vertices_begin(); vi != pg.vertices_end(); ++vi, ++loc)
	{
		VertexMemo vmemo(*vi, loc, pgId);

		vmemo.set_normal(normal_vec[loc]);
		vmemo.set_sdf(sdf_vec[loc]);
		vmemo.set_ref_point(vertex_vec[loc] - sdf_vec[loc] * normal_vec[loc] / 2.0);
		vmemo.compute_force(K_SDF, K_S, K_D, T1);

		vmemo_vec.push_back(vmemo);
	}

	size_t loc1 = startIdx, loc2 = startIdx + 1;
	for (Polygon::Edge_const_iterator ei = pg.edges_begin(); ei != pg.edges_end(); ++ei, ++loc1, ++loc2)
	{
		if ((loc2 - startIdx) == pg.size())
			loc2 = startIdx;

		if (sdf_vec[loc1] <= T1)
			setMovable(pg, startIdx, loc1, 3);

		EdgeMemo ememo(&(vmemo_vec[loc1]), &(vmemo_vec[loc2]));
		ememo.compute_force(K_S, K_D);

		ememo_vec.push_back(ememo);
	}
}

// compute element forces and Jacobians on vertex elements
size_t populate_memos(
	const Pwh& pwh,
	bool recomputeSDF = false
)
{
	vmemo_vec.clear();
	ememo_vec.clear();

	if (recomputeSDF || sdf_vec.size() == 0)
		computeSDF(pwh, !recomputeSDF);

	size_t startIdx = 0, pgId = 0;
	populate_memos(pwh.outer_boundary(), startIdx, pgId);
	startIdx += pwh.outer_boundary().size();
	for (Pwh::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); hi++)
	{
		populate_memos(*hi, startIdx, ++pgId);
		startIdx += hi->size();
	}

	size_t movable = 0;
	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->isMovable)
			it->index = ++movable;
		else
			it->index = 0;
	}

	return movable;
}

// return location indicies of 2x2 submatrix starting at (row, col) in 2 row location format
arma::umat getLocationIndices(
	int row,
	int col
)
{
	arma::umat loc;
	loc << row << row + 1 << row << row + 1 << arma::endr
		<< col << col << col + 1 << col + 1 << arma::endr;
	return loc;
}

// assembly
void global_assembly(
	arma::vec& FORCE,
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
	arma::vec& VELOCITY,
	const size_t& N
)
{
	FORCE.set_size(N * 2);
	VELOCITY.set_size(N * 2);

	FORCE.fill(0.0);
	VELOCITY.fill(0.0);

	arma::umat locations;
	arma::vec JPOSvalues, JVELvalues;

	int s_idx, e_idx, count = 0;
	// compute the number of entries in sparse Jacobian matrices and fill FORCE and VELOCITY
	for (std::vector<EdgeMemo>::iterator it = ememo_vec.begin(); it != ememo_vec.end(); it++)
	{
		s_idx = ((int)it->source->index - 1) * 2;
		e_idx = ((int)it->target->index - 1) * 2;

		if (s_idx >= 0)
		{
			FORCE(arma::span(s_idx, s_idx + 1)) += it->get_force();
			count += 4;
		}
		if (e_idx >= 0) // if the opposite vertex is also movable
		{
			FORCE(arma::span(e_idx, e_idx + 1)) -= it->get_force();
			count += 4;
		}
		if (s_idx >= 0 && e_idx >= 0)
			count += 8;
	}

	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->index > 0)
		{
			s_idx = (it->index - 1) * 2;
			FORCE(arma::span(s_idx, s_idx + 1)) += it->get_force();
			VELOCITY(arma::span(s_idx, s_idx + 1)) = it->get_velocity(); //TODO: this is unnecessary during updation
			count += 4;
		}
	}

	locations.set_size(2, count);
	JPOSvalues.set_size(count);
	JVELvalues.set_size(count);

	count = 0;
	// fill sparse matrices JPOS and JVEL
	for (std::vector<EdgeMemo>::iterator it = ememo_vec.begin(); it != ememo_vec.end(); it++)
	{
		s_idx = ((int)it->source->index - 1) * 2;
		e_idx = ((int)it->target->index - 1) * 2;

		if (s_idx >= 0)
		{
			JPOSvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 3)) = getLocationIndices(s_idx, s_idx);
			count += 4;
		}
		if (e_idx >= 0) // if the opposite vertex is also movable
		{
			JPOSvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 3)) = getLocationIndices(e_idx, e_idx);
			count += 4;
		}
		if (s_idx >= 0 && e_idx >= 0)
		{
			JPOSvalues(arma::span(count, count + 3)) = -arma::vectorise(it->get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 3)) = -arma::vectorise(it->get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 3)) = getLocationIndices(s_idx, e_idx);
			count += 4;

			JPOSvalues(arma::span(count, count + 3)) = -arma::vectorise(it->get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 3)) = -arma::vectorise(it->get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 3)) = getLocationIndices(e_idx, s_idx);
			count += 4;
		}
	}

	// loop over each vertex
	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->index > 0)
		{
			s_idx = (it->index - 1) * 2;

			JPOSvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 3)) = arma::vectorise(it->get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 3)) = getLocationIndices(s_idx, s_idx);
			count += 4;
		}
	}

	JPOS = arma::sp_mat(true, locations, JPOSvalues, N * 2, N * 2);
	JVEL = arma::sp_mat(true, locations, JVELvalues, N * 2, N * 2);
}


// compute update in velocity and position for one time step
void march_one_time_step(
	arma::vec& delta_VEL,
	arma::vec& delta_POS,
	const size_t& N,
	const arma::vec& FORCE,
	const arma::sp_mat& JPOS,
	const arma::sp_mat& JVEL,
	const arma::vec& VELOCITY
)
{
	arma::sp_mat A;
	arma::vec b;

	A = arma::speye<arma::sp_mat>(N * 2, N * 2) - TIME_STEP*MASS_INV*JVEL - TIME_STEP*TIME_STEP*MASS_INV*JPOS;
	b = TIME_STEP*MASS_INV*(FORCE + TIME_STEP*JPOS*VELOCITY);

	delta_VEL = arma::spsolve(A, b); // using SuperLU solver, TODO: use settings to compute faster, such as symmetric
	//delta_VEL = arma::solve(arma::mat(A), b);
	//delta_VEL = arma::inv_sympd(arma::mat(A))*b;
	delta_POS = TIME_STEP*(VELOCITY + delta_VEL);
}

// create a new Pwh from vmemo_vec
Pwh createPwhFromVmemo(Pwh &pwh)
{
	pwh.clear();
	pwh.outer_boundary();
	Polygon pg;
	size_t id = 0;
	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->pgId != id)
		{
			if (id == 0)
				pwh = Pwh(Polygon(pg));
			else
				pwh.add_hole(Polygon(pg));
			id++;
			pg.clear();
		}
		pg.push_back(it->get_vertex());
	}
	if (id == 0)
		pwh = Pwh(Polygon(pg));
	else
		pwh.add_hole(Polygon(pg));
	return pwh;
}

// create a new Pwh from vmemo_vec
Pwh createPwhFromVmemo()
{
	Pwh pwh;
	Polygon pg;
	size_t id = 0;
	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->pgId != id)
		{
			if (id == 0)
				pwh = Pwh(Polygon(pg));
			else
				pwh.add_hole(Polygon(pg));
			id++;
			pg.clear();
		}
		pg.push_back(it->get_vertex());
	}
	if (id == 0)
		pwh = Pwh(Polygon(pg));
	else
		pwh.add_hole(Polygon(pg));
	return pwh;
}



// update the positions, velocities, and memos
void update(
	const arma::vec& delta_POS,
	const arma::vec& delta_VEL,
	const size_t& N,
	Pwh& pwh,
	arma::vec& FORCE,
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
	arma::vec& VELOCITY
)
{
	size_t idx;
	Vector_2 del;
	for (std::vector<VertexMemo>::iterator it = vmemo_vec.begin(); it != vmemo_vec.end(); it++)
	{
		if (it->index > 0)
		{
			idx = (it->index - 1) * 2;
			del = Vector_2(delta_POS(idx), delta_POS(idx + 1));
			it->set_vertex(it->get_vertex() + del);
		}
	}

	createPwhFromVmemo(pwh);

	// TODO: decide whether to update normals and sdf_force every time step or not
	// updating normals
	getVertexVector(pwh);
	getNormalsAndAngles(pwh);

	for (size_t i = 0; i < vmemo_vec.size(); i++)
	{
		vmemo_vec[i].set_normal(normal_vec[i]);
		vmemo_vec[i].compute_force(K_SDF, K_S, K_D, T1);
		ememo_vec[i].compute_force(K_S, K_D);
	}

	global_assembly(FORCE, JPOS, JVEL, VELOCITY, N);
}

// computes the maxima of change of location of vertices
double max_change(
	arma::vec delta_POS
)
{
	size_t N = delta_POS.n_rows;
	arma::mat delta_XY = delta_POS;
	arma::vec dist;
	delta_XY = arma::square(delta_XY); // square each element
	delta_XY.reshape(2, N / 2); // reshape into 2x(N/2) matrix with squared [x y] coordinate change
	arma::inplace_trans(delta_XY); // transpose
	dist = arma::sum(delta_XY, 1); // sum elements in each row
	dist = arma::sqrt(dist); // square root of each element
	return dist.max();
}

// march one step and then update
double march_and_update(
	const size_t& N,
	Pwh& pwh,
	arma::vec& FORCE,
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
	arma::vec& VELOCITY
)
{
	arma::vec delta_VEL;
	arma::vec delta_POS;
	march_one_time_step(delta_VEL, delta_POS, N, FORCE, JPOS, JVEL, VELOCITY);
	update(delta_POS, delta_VEL, N, pwh, FORCE, JPOS, JVEL, VELOCITY);
	return max_change(delta_POS);
}





// TODO: set parameters for resolution and resampling

int main(int argc, char *argv[])
{
	re.seed(0);

	Pwh_list slice;
	read(FILENAME, slice);

	std::cout << "Starting resampling..." << std::endl;
	Pwh pwh = resample(slice.front(), T1 / 4.0);
	std::cout << "Resampling done!" << std::endl;

	// populate memos - compute and store normals, sdf, forces, and Jacobians
	std::cout << "Populating memos..." << std::endl;
	size_t nMove = populate_memos(pwh, true); // nMove is the number of movable vertices
	std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
	std::cout << std::endl;

	arma::vec max_change_vec(STEPS);
	if (STEPS > 0)
	{
		arma::vec FORCE;
		arma::sp_mat JPOS;
		arma::sp_mat JVEL;
		arma::vec VELOCITY;

		// assemble
		std::cout << "Assembling..." << std::endl;
		global_assembly(FORCE, JPOS, JVEL, VELOCITY, nMove);
		std::cout << "Assembling completed!" << std::endl;
		std::cout << std::endl;

		// numerical time integration
		max_change_vec.fill(0.0);
		std::cout << "Performing numerical integration..." << std::endl;
		for (size_t i = 0; i < STEPS; i++)
		{
			std::cout << "Time step " << i << std::endl;
			max_change_vec[i] = march_and_update(nMove, pwh, FORCE, JPOS, JVEL, VELOCITY);
			std::cout << "Max change = " << max_change_vec[i] << std::endl;
		}
		std::cout << "Integration completed!" << std::endl;
		std::cout << std::endl;

		// populate memos - compute and store normals, sdf, forces, and Jacobians
		std::cout << "Populating memos..." << std::endl;
		//nMove = populate_memos(pwh, false); // nMove is the number of movable vertices
		std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
		std::cout << std::endl;
	}

	//intersect(pwh.outer_boundary(), Point_2(-3.75, -2.375), Vector_2(-1.0, 0.0));


	std::ofstream outfile("data/exampleOut.dat");
	outfile << write(pwh);
	outfile.close();

	return 0;
}