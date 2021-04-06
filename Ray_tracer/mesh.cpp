#include "mesh.h"
#include "plane.h"
#include <fstream>
#include <string>
#include <limits>

// Consider a triangle to intersect a ray if the ray intersects the plane of the
// triangle with barycentric weights in [-weight_tolerance, 1+weight_tolerance]
static const double weight_tolerance = 1e-4;

// Read in a mesh from an obj file.  Populates the bounding box and registers
// one part per triangle (by setting number_parts).
void Mesh::Read_Obj(const char* file)
{
    std::ifstream fin(file);
    if(!fin)
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    ivec3 e;
    vec3 v;
    box.Make_Empty();
    while(fin)
    {
        getline(fin,line);

        if(sscanf(line.c_str(), "v %lg %lg %lg", &v[0], &v[1], &v[2]) == 3)
        {
            vertices.push_back(v);
            box.Include_Point(v);
        }

        if(sscanf(line.c_str(), "f %d %d %d", &e[0], &e[1], &e[2]) == 3)
        {
            for(int i=0;i<3;i++) e[i]--;
            triangles.push_back(e);
        }
    }
    number_parts=triangles.size();
}

// Check for an intersection against the ray.  See the base class for details.
Hit Mesh::Intersection(const Ray& ray, int part) const
{
    // TODO;
	/* std::cout << "vertices size: " << vertices.size() << std::endl;
	for (unsigned i = 0; i < vertices.size(); ++i) {
		std::cout << vertices.at(i) << ", ";
	}
	std::cout << std::endl;
	std::cout << "triangles size: " << triangles.size() << std::endl;
	for (unsigned i = 0; i < triangles.size(); ++i) {
		std::cout << triangles.at(i) << ", ";
	}
	std::cout << std::endl; */

	double a = 0.0;
	bool intersect = Intersect_Triangle(ray, part, a);
	if (intersect) {
		return {this, a, part};
	
	}
    return {NULL, a, part};
}

// Compute the normal direction for the triangle with index part.
vec3 Mesh::Normal(const vec3& point, int part) const
{
    assert(part>=0);
    // TODO;
	ivec3 triangle = triangles.at(part);
	vec3 first = vertices.at(triangle[1]) - vertices.at(triangle[0]);
	vec3 second = vertices.at(triangle[2]) - vertices.at(triangle[0]);
	vec3 normal = cross(first, second).normalized();
    // return vec3();
	return normal;
}

// This is a helper routine whose purpose is to simplify the implementation
// of the Intersection routine.  It should test for an intersection between
// the ray and the triangle with index tri.  If an intersection exists,
// record the distance and return true.  Otherwise, return false.
// This intersection should be computed by determining the intersection of
// the ray and the plane of the triangle.  From this, determine (1) where
// along the ray the intersection point occurs (dist) and (2) the barycentric
// coordinates within the triangle where the intersection occurs.  The
// triangle intersects the ray if dist>small_t and the barycentric weights are
// larger than -weight_tolerance.  The use of small_t avoid the self-shadowing
// bug, and the use of weight_tolerance prevents rays from passing in between
// two triangles.
bool Mesh::Intersect_Triangle(const Ray& ray, int tri, double& dist) const
{
    // TODO;
	ivec3 triangle = triangles.at(tri);
	vec3 first = vertices.at(triangle[1]) - vertices.at(triangle[0]);
	vec3 second = vertices.at(triangle[2]) - vertices.at(triangle[0]);
	vec3 normal = cross(first, second).normalized();
	Plane plane(vertices.at(0), normal);
	Hit h = plane.Intersection(ray, tri);
	if (h.object == 0) {
		return false;
	}
	
	dist = h.dist;
	vec3 intersection_point = ray.Point(dist);

	double areaOfFullTriangle = 0.5 * (cross(first, second)).magnitude();
	vec3 side = intersection_point - vertices.at(triangle[0]);

	double area2 = 0.5 * dot(cross(first, side), normal);
	double area1 = 0.5 * dot(cross(side, second), normal);
	double two = area2 / areaOfFullTriangle;
	double one = area1 / areaOfFullTriangle;
	double zero = 1 - two - one;



	/* double area2_0 = 0.5 * cross(first, side).magnitude();
	double area1_0 = 0.5 * cross(second, side).magnitude();
	double two_0 = area2 / areaOfFullTriangle;
	double one_0 = area1 / areaOfFullTriangle;
	double zero_0 = 1 - two - one;

		std::cout << "intersection_point: " << intersection_point << std::endl;
		std::cout << "side: " << side << std::endl;
		std::cout << "normal: " << normal << std::endl;
		std::cout << "area2: " << area2 << std::endl;
		std::cout << "area1: " << area1 << std::endl;
		std::cout << "one: " << one << " two: " << two << " zero: " << zero << std::endl;
		std::cout << "area2_0: " << area2_0 << std::endl;
		std::cout << "area1_0: " << area1_0 << std::endl;
		std::cout << "one_0: " << one_0 << " two_0: " << two_0 << " zero_0: " << zero_0 << std::endl; */


	if (two < -weight_tolerance || one < -weight_tolerance || zero < -weight_tolerance) {
		return false;
	}
	
	
	return true;
}

// Compute the bounding box.  Return the bounding box of only the triangle whose
// index is part.
Box Mesh::Bounding_Box(int part) const
{
    Box b;
    TODO;
    return b;
}
