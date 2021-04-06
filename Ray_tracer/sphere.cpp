#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    // TODO;
	double discriminant = dot(ray.direction, ray.endpoint - center) * dot(ray.direction, ray.endpoint - center) - (ray.endpoint - center).magnitude_squared() + radius * radius;
	if (discriminant < 0) {
		return {NULL, 0, 0};
	}
	else if (discriminant == 0) {
		double t = -dot(ray.direction, ray.endpoint - center);
		return {this, t, part};
	}
	else {
		double t1, t2;
		t1 = -dot(ray.direction, ray.endpoint - center) + sqrt(discriminant);
		t2 = -dot(ray.direction, ray.endpoint - center) - sqrt(discriminant);
		if (t1 > small_t && t1 < t2) {
			return {this, t1, part};
		}
		else if (t2 > small_t && t2 < t1) {
			return {this, t2, part};
		}
		else if (t1 < 0 && t2 > small_t) {
			return {this, t2, part};
		}
		else if (t2 < 0 && t1 > small_t) {
			return {this, t1, part};
		}
		else {
			return {NULL, 0, 0};
		}
	}
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    // TODO; // compute the normal direction
	normal = (point - center).normalized();
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    return box;
}
