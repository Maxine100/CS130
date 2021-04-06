#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    TODO; // determine the color
	Ray reflected_ray;
	reflected_ray.endpoint = intersection_point;
	reflected_ray.direction = (ray.direction - normal * 2 * dot(ray.direction, normal)).normalized();

	color = shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);
	color = (1 - reflectivity) * color + reflectivity * shader->world.Cast_Ray(reflected_ray, recursion_depth);

    return color;
}
