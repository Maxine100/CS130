#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    // TODO; //determine the color
	vec3 Ia, Id, Is;
	Ia = color_ambient * world.ambient_color * world.ambient_intensity;

	for (unsigned i = 0; i < world.lights.size(); ++i) {

		vec3 light_ray = (world.lights.at(i)->position - intersection_point);

		vec3 L = world.lights.at(i)->Emitted_Light(light_ray);

		double ndotlight_ray = dot(normal, light_ray.normalized());
		if (ndotlight_ray < 0) {
			ndotlight_ray = 0;
		}

		vec3 reflected = (-light_ray + normal * 2 * dot(light_ray, normal)).normalized();
		double vdotreflected = dot((world.camera.position - intersection_point).normalized(), reflected);
		if (vdotreflected < 0) {
			vdotreflected = 0;
		}



		

		// Compute shadow ray.
		// If !shadow ray hits an object add light's diffuse and specular components.

		if (world.enable_shadows) {
			Ray shadow;
			shadow.endpoint = intersection_point;
			shadow.direction = (world.lights.at(i)->position - intersection_point);
			double shadow_extent = shadow.direction.magnitude();
			shadow.direction = shadow.direction.normalized();

			Hit hit = world.Closest_Intersection(shadow);
			
			vec3 intersectionOfShadowAndObject;
			intersectionOfShadowAndObject = shadow.direction * hit.dist;
			if (hit.object != 0 && shadow_extent > intersectionOfShadowAndObject.magnitude()) {
			
			}
			else {
				Id += color_diffuse * L * ndotlight_ray;
				Is += color_specular * L * pow(vdotreflected, specular_power);
			}


		}
		else {
			Id += color_diffuse * L * ndotlight_ray;
		 	Is += color_specular * L * pow(vdotreflected, specular_power);
		}
		
	color = Ia + Id + Is;
	}

    return color;
}
