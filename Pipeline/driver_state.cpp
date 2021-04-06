#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    // state.image_color=0;
    state.image_depth=0;
	state.image_color = new pixel[width * height];
	state.image_depth = new float[width * height];
    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
	for (int i = 0; i < width * height; ++i) {
		state.image_color[i] = make_pixel(0, 0, 0);
		state.image_depth[i] = 100.0;
	}
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;

	int startOfFirstVertex, startOfSecondVertex, startOfThirdVertex;
	data_geometry* pointer[3];
	data_geometry* copy = new data_geometry[3];
	data_vertex in;
	const data_geometry* constCopy[3];

	switch (type) {
		case render_type::invalid: {

			break;
		}
		case render_type::indexed: {
			// std::cout << "state.num_vertices: " << state.num_vertices << "; state.floats_per_vertex: " << state.floats_per_vertex << std::endl;
			

			for (int i = 0; i < state.num_triangles; ++i) {
				// std::cout << "Triangle #" << i + 1 << std::endl;

				startOfFirstVertex = state.index_data[3 * i];
				startOfSecondVertex = state.index_data[3 * i + 1];
				startOfThirdVertex = state.index_data[3 * i + 2];

				/* std::cout << "	Vertices: " << startOfFirstVertex << " " << startOfSecondVertex << " " << startOfThirdVertex << std::endl;
				
				std::cout << "		" << state.vertex_data[3 * startOfFirstVertex] << " " << state.vertex_data[3 * startOfFirstVertex + 1] << " " << state.vertex_data[3 * startOfFirstVertex + 2] << std::endl;
				std::cout << "		" << state.vertex_data[3 * startOfSecondVertex] << " " << state.vertex_data[3 * startOfSecondVertex + 1] << " " << state.vertex_data[3 * startOfSecondVertex + 2] << std::endl;
				std::cout << "		" << state.vertex_data[3 * startOfThirdVertex] << " " << state.vertex_data[3 * startOfThirdVertex + 1] << " " << state.vertex_data[3 * startOfThirdVertex + 2] << std::endl; */

				data_geometry first, second, third;
				first.data = &state.vertex_data[3 * startOfFirstVertex];
				second.data = &state.vertex_data[3 * startOfSecondVertex];
				third.data = &state.vertex_data[3 * startOfThirdVertex];

				pointer[0] = &first;
				pointer[1] = &second;
				pointer[2] = &third;

				/* std::cout << "		";
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < state.floats_per_vertex; ++k) {
						std::cout << "pointer[" << j << "]->data[" << k << "]: " << pointer[j]->data[k] << "; ";
					}
				} */

				for (int k = 0; k < 3; ++k) {
					copy[k].data = pointer[k]->data;
					in.data = pointer[k]->data;
					state.vertex_shader(in, copy[k], state.uniform_data);
				}

				for (int k = 0; k < 3; ++k) {
					constCopy[k] = &copy[k];
				}

				clip_triangle(state, constCopy, 0);
			}

			break;
		}
		case render_type::triangle: {
			/* std::cout << "state.num_vertices: " << state.num_vertices << "; state.floats_per_vertex: " << state.floats_per_vertex << std::endl << std::endl;

			for (int i = 0; i < state.num_vertices * state.floats_per_vertex; ++i) {
				// std::cout << "state.vertex_data[" << i << "]: " << state.vertex_data[i] << "; ";
			}
			std::cout << std::endl << std::endl; */

			for (int i = 0; i < state.num_vertices * state.floats_per_vertex; i += 3 * state.floats_per_vertex) {

				// std::cout << "i: " << i << std::endl;

				startOfFirstVertex = i;
				startOfSecondVertex = i + 1 * state.floats_per_vertex;
				startOfThirdVertex = i + 2 * state.floats_per_vertex;
				// std::cout << "startOfFirstVertex: " << startOfFirstVertex << "; startOfSecondVertex: " << startOfSecondVertex << "; startOfThirdVertex: " << startOfThirdVertex << std::endl << std::endl;

				data_geometry first, second, third;
				first.data = &state.vertex_data[startOfFirstVertex];
				second.data = &state.vertex_data[startOfSecondVertex];
				third.data = &state.vertex_data[startOfThirdVertex];

				/* for (int j = 0; j < state.floats_per_vertex; ++j) {
					std::cout << "first.data[" << j << "]: " << first.data[j] << "; ";
				}
				std::cout << std::endl;
				for (int j = 0; j < state.floats_per_vertex; ++j) {
					std::cout << "second.data[" << j << "]: " << second.data[j] << "; ";
				}
				std::cout << std::endl;
				for (int j = 0; j < state.floats_per_vertex; ++j) {
					std::cout << "third.data[" << j << "]: " << third.data[j] << "; ";
				}
				std::cout << std::endl << std::endl; */
				pointer[0] = &first;
				pointer[1] = &second;
				pointer[2] = &third;

				/* for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < state.floats_per_vertex; ++k) {
						std::cout << "pointer[" << j << "]->data[" << k << "]: " << pointer[j]->data[k] << "; ";
					}
				}
				std::cout << std::endl << std::endl; */
				
				for (int k = 0; k < 3; ++k) {
					copy[k].data = pointer[k]->data;
					in.data = pointer[k]->data;
					/* std::cout << "in: ";
					for (int j = 0; j < state.floats_per_vertex; ++j) {
						std::cout << in.data[j] << ", ";
					}
					std::cout << std::endl; */
					state.vertex_shader(in, copy[k], state.uniform_data);
				}
				/* for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < state.floats_per_vertex; ++k) {
						std::cout << "copy[" << j << "].data[" << k << "]: " << copy[j].data[k] << "; ";
					}
				}
				std::cout << std::endl << std::endl; */

				for (int k = 0; k < 3; ++k) {
					constCopy[k] = &copy[k];
				}

				/* for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < state.floats_per_vertex; ++k) {
						std::cout << "constCopy[" << j << "]->data[" << k << "]: " << constCopy[j]->data[k] << "; ";
					}
				}
				std::cout << std::endl << std::endl; */

				clip_triangle(state, constCopy, 0);
				// rasterize_triangle(state, constCopy);
			}
			break;
		}
		case render_type::fan: {
			/* std::cout << "state.num_vertices: " << state.num_vertices << "; state.floats_per_vertex: " << state.floats_per_vertex << std::endl << std::endl;
			for (int i = 0; i < state.num_vertices * state.floats_per_vertex; ++i) {
				std::cout << "state.vertex_data[" << i << "]: " << state.vertex_data[i] << "; ";
			}
			std::cout << std::endl << std::endl; */

			data_geometry center;
			center.data = &state.vertex_data[0];

			for (int i = state.floats_per_vertex; i < state.floats_per_vertex * state.num_vertices - state.floats_per_vertex; i += state.floats_per_vertex) {
				// std::cout << "i: " << i << std::endl;

				startOfFirstVertex = i;
				startOfSecondVertex = i + state.floats_per_vertex;

				// std::cout << "startOfFirstVertex: " << startOfFirstVertex << "; startOfSecondVertex: " << startOfSecondVertex << std::endl;

				data_geometry first, second;
				first.data = &state.vertex_data[startOfFirstVertex];
				second.data = &state.vertex_data[startOfSecondVertex];

				pointer[0] = &center;
				pointer[1] = &first;
				pointer[2] = &second;

				for (int k = 0; k < 3; ++k) {
					copy[k].data = pointer[k]->data;
					in.data = pointer[k]->data;
					state.vertex_shader(in, copy[k], state.uniform_data);
				}
				for (int k = 0; k < 3; ++k) {
					constCopy[k] = &copy[k];
				}
				clip_triangle(state, constCopy, 0);
			}
			break;
		}
		case render_type::strip: {
			/* std::cout << "state.num_vertices: " << state.num_vertices << "; state.floats_per_vertex: " << state.floats_per_vertex << std::endl << std::endl;
			for (int i = 0; i < state.num_vertices * state.floats_per_vertex; ++i) {
				std::cout << "state.vertex_data[" << i << "]: " << state.vertex_data[i] << "; ";
			}
			std::cout << std::endl << std::endl; */

			for (int i = 0; i < state.floats_per_vertex * state.num_vertices - 2 * state.floats_per_vertex; i += state.floats_per_vertex) {
				// std::cout << "i: " << i << std::endl;
				startOfFirstVertex = i;
				startOfSecondVertex = i + state.floats_per_vertex;
				startOfThirdVertex = i + 2 * state.floats_per_vertex;
				// std::cout << "startOfFirstVertex: " << startOfFirstVertex << "; startOfSecondVertex: " << startOfSecondVertex << "; startOfThirdVertex: " << startOfThirdVertex << std::endl;

				data_geometry first, second, third;
				first.data = &state.vertex_data[startOfFirstVertex];
				second.data = &state.vertex_data[startOfSecondVertex];
				third.data = &state.vertex_data[startOfThirdVertex];

				pointer[0] = &first;
				pointer[1] = &second;
				pointer[2] = &third;

				for (int k = 0; k < 3; ++k) {
					copy[k].data = pointer[k]->data;
					in.data = pointer[k]->data;
					state.vertex_shader(in, copy[k], state.uniform_data);
				}
				for (int k = 0; k < 3; ++k) {
					constCopy[k] = &copy[k];
				}
				clip_triangle(state, constCopy, 0);
			}
			break;
		}
	}
}



/* void interpolate() {








} */



// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
	int amountOfIns = 0;
	bool firstIn = false;
	bool secondIn = false;
	bool thirdIn = false;
	double delta;
	const data_geometry* goesIn[3];
	data_geometry* newOne = new data_geometry;
	data_geometry* newTwo = new data_geometry;
	data_geometry* newVertices[3];
	for (int i = 0; i < 3; ++i) {
		float* newArray = new float[state.floats_per_vertex];
		newVertices[i] = (data_geometry*)(in[i]);
		for (int j = 0; j < state.floats_per_vertex; ++j) {
			newArray[j] = in[i]->data[j];
		}
		newVertices[i]->data = newArray;
	}

	/* for (int i = 0; i < 3; ++i) {
		std::cout << "newVertices[" << i << "]:" << std::endl;
		std::cout << "	gl_Position: ";
		for (int j = 0; j < 4; ++j) {
			std::cout << newVertices[i]->gl_Position[j] << " ";
		}
		std::cout << std::endl;
		std::cout << "	data: ";
		for (int j = 0; j < state.floats_per_vertex; ++j) {
			std::cout << newVertices[i]->data[j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl; */

	int index = face / 2;
	int sign = 1;
	if (face % 2 == 1) {
		sign = -1;
	}

	// std::cout << "face: " << face << "; index: " << index << "; sign: " << sign << std::endl;
	

    if(face==6)
    {
	// std::cout << "face = 6" << std::endl;
        rasterize_triangle(state, in);
        return;
    }
	

	if (sign * newVertices[0]->gl_Position[index] <= newVertices[0]->gl_Position[3]) {
		// std::cout << "		first is inside." << std::endl;
		firstIn = true;
		++amountOfIns;
	}
	if (sign * newVertices[1]->gl_Position[index] <= newVertices[1]->gl_Position[3]) {
		// std::cout << "		second is inside." << std::endl;
		secondIn = true;
		++amountOfIns;
	}
	if (sign * newVertices[2]->gl_Position[index] <= newVertices[2]->gl_Position[3]) {
		// std::cout << "		third is inside." << std::endl;
		thirdIn = true;
		++amountOfIns;
	}
	
	if (amountOfIns == 0) {
		// std::cout << "		None are inside." << std::endl;
		return;
	}
	else if (amountOfIns == 1) {
		// std::cout << "		One is inside." << std::endl;
		if (firstIn) {
			// std::cout << "			first is inside." << std::endl;
			delta = (sign * newVertices[1]->gl_Position[3] - newVertices[1]->gl_Position[index]) / (newVertices[0]->gl_Position[index] - sign * newVertices[0]->gl_Position[3] - newVertices[1]->gl_Position[index] + sign * newVertices[1]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[1]->gl_Position[i] = delta * newVertices[0]->gl_Position[i] + (1 - delta) * newVertices[1]->gl_Position[i];
				// newVertices[1]->data[i] = delta * newVertices[0]->data[i] + (1 - delta) * newVertices[1]->data[i];
			}
			delta = (sign * newVertices[2]->gl_Position[3] - newVertices[2]->gl_Position[index]) / (newVertices[0]->gl_Position[index] - sign * newVertices[0]->gl_Position[3] - newVertices[2]->gl_Position[index] + sign * newVertices[2]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[2]->gl_Position[i] = delta * newVertices[0]->gl_Position[i] + (1 - delta) * newVertices[2]->gl_Position[i];
				// newVertices[2]->data[i] = delta * newVertices[0]->data[i] + (1 - delta) * newVertices[2]->data[i];
			}
		}
		if (secondIn) {
			// std::cout << "			second is inside." << std::endl;
			delta = (sign * newVertices[0]->gl_Position[3] - newVertices[0]->gl_Position[index]) / (newVertices[1]->gl_Position[index] - sign * newVertices[1]->gl_Position[3] - newVertices[0]->gl_Position[index] + sign* newVertices[0]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[0]->gl_Position[i] = delta * newVertices[1]->gl_Position[i] + (1 - delta) * newVertices[0]->gl_Position[i];
				// newVertices[0]->data[i] = delta * newVertices[1]->data[1] + (1 - delta) * newVertices[0]->data[i];
			}
			delta = (sign * newVertices[2]->gl_Position[3] - newVertices[2]->gl_Position[index]) / (newVertices[1]->gl_Position[index] - sign * newVertices[1]->gl_Position[3] - newVertices[2]->gl_Position[index] + sign * newVertices[2]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[2]->gl_Position[i] = delta * newVertices[1]->gl_Position[i] + (1 - delta) * newVertices[2]->gl_Position[i];
				// newVertices[2]->data[i] = delta * newVertices[1]->data[i] + (1 - delta) * newVertices[2]->data[i];
			}
		}
		if (thirdIn) {
			// std::cout << "			third is inside." << std::endl;
			delta = (sign * newVertices[0]->gl_Position[3] - newVertices[0]->gl_Position[index]) / (newVertices[2]->gl_Position[index] - sign * newVertices[2]->gl_Position[3] - newVertices[0]->gl_Position[index] + sign* newVertices[0]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[0]->gl_Position[i] = delta * newVertices[2]->gl_Position[i] + (1 - delta) * newVertices[0]->gl_Position[i];
				// newVertices[0]->data[i] = delta * newVertices[2]->data[i] + (1 - delta) * newVertices[0]->data[i];
			}
			delta = (sign * newVertices[1]->gl_Position[3] - newVertices[1]->gl_Position[index]) / (newVertices[2]->gl_Position[index] - sign * newVertices[2]->gl_Position[3] - newVertices[1]->gl_Position[index] + sign* newVertices[1]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newVertices[1]->gl_Position[i] = delta * newVertices[2]->gl_Position[i] + (1 - delta) * newVertices[1]->gl_Position[i];
				// newVertices[1]->data[i] = delta * newVertices[2]->data[i] + (1 - delta) * newVertices[0]->data[i];
			}
		}
		for (int i = 0; i < 3; ++i) {
			goesIn[i] = newVertices[i];
		}
		clip_triangle(state, goesIn, face + 1);
	}




	else if (amountOfIns == 2) {
		// std::cout << "		Two are inside." << std::endl;
		if (firstIn && secondIn) {
			// std::cout << "			first and second are inside." << std::endl;
			newOne->gl_Position = newVertices[2]->gl_Position;
			float* newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[2]->data[i];
			}
			newOne->data = newarray;
			newTwo->gl_Position = newVertices[2]->gl_Position;
			newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[2]->data[i];
			}
			newTwo->data = newarray;
		
			delta = (sign * newVertices[2]->gl_Position[3] - newVertices[2]->gl_Position[index]) / (newVertices[0]->gl_Position[index] - sign * newVertices[0]->gl_Position[3] - newVertices[2]->gl_Position[index] + sign * newVertices[2]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newOne->gl_Position[i] = delta * newVertices[0]->gl_Position[i] + (1 - delta) * newOne->gl_Position[i];
				// newOne->data[i] = delta * newVertices[0]->data[i] + (1 - delta) * newOne->data[i];
			}
			delta = (sign * newVertices[2]->gl_Position[3] - newVertices[2]->gl_Position[index]) / (newVertices[1]->gl_Position[index] - sign * newVertices[1]->gl_Position[3] - newVertices[2]->gl_Position[index] + sign * newVertices[2]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newTwo->gl_Position[i] = delta * newVertices[1]->gl_Position[i] + (1 - delta) * newTwo->gl_Position[i];
				// newTwo->data[i] = delta * newVertices[1]->data[i] + (1 - delta) * newTwo->data[i];
			}
			goesIn[0] = newVertices[0];
			goesIn[1] = newVertices[1];
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
			goesIn[0] = newVertices[0];
			goesIn[1] = newOne;
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
		}
		if (firstIn && thirdIn) {
			// std::cout << "			first and third are inside." << std::endl;
			newOne->gl_Position = newVertices[1]->gl_Position;
			float* newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[1]->data[i];
			}
			newOne->data = newarray;
			newTwo->gl_Position = newVertices[1]->gl_Position;
			newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[1]->data[i];
			}
			newTwo->data = newarray;
	
			delta = (sign * newVertices[1]->gl_Position[3] - newVertices[1]->gl_Position[index]) / (newVertices[0]->gl_Position[index] - sign * newVertices[0]->gl_Position[3] - newVertices[1]->gl_Position[index] + sign * newVertices[1]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newOne->gl_Position[i] = delta * newVertices[0]->gl_Position[i] + (1 - delta) * newOne->gl_Position[i];
				// newOne->data[i] = delta * newVertices[0]->data[i] + (1 - delta) * newOne->data[i];
			}
			delta = (sign * newVertices[1]->gl_Position[3] - newVertices[1]->gl_Position[index]) / (newVertices[2]->gl_Position[index] - sign * newVertices[2]->gl_Position[3] - newVertices[1]->gl_Position[index] + sign* newVertices[1]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newTwo->gl_Position[i] = delta * newVertices[2]->gl_Position[i] + (1 - delta) * newTwo->gl_Position[i];
				// newTwo->data[i] = delta * newVertices[2]->data[i] + (1 - delta) * newTwo->data[i];
			}
			goesIn[0] = newVertices[0];
			goesIn[1] = newVertices[2];
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
			goesIn[0] = newVertices[0];
			goesIn[1] = newOne;
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
		}
		if (secondIn && thirdIn) {
			// std::cout << "			second and third are inside." << std::endl;
			newOne->gl_Position = newVertices[0]->gl_Position;
			float* newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[0]->data[i];
			}
			newOne->data = newarray;
			newTwo->gl_Position = newVertices[0]->gl_Position;
			newarray = new float[state.floats_per_vertex];
			for (int i = 0; i < state.floats_per_vertex; ++i) {
				newarray[i] = newVertices[0]->data[i];
			}
			newTwo->data = newarray;
			
			delta = (sign * newVertices[0]->gl_Position[3] - newVertices[0]->gl_Position[index]) / (newVertices[1]->gl_Position[index] - sign * newVertices[1]->gl_Position[3] - newVertices[0]->gl_Position[index] + sign* newVertices[0]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newOne->gl_Position[i] = delta * newVertices[1]->gl_Position[i] + (1 - delta) * newOne->gl_Position[i];
				// newOne->data[i] = delta * newVertices[1]->data[i] + (1 - delta) * newOne->data[i];
			}
			delta = (sign * newVertices[0]->gl_Position[3] - newVertices[0]->gl_Position[index]) / (newVertices[2]->gl_Position[index] - sign * newVertices[2]->gl_Position[3] - newVertices[0]->gl_Position[index] + sign* newVertices[0]->gl_Position[3]);
			for (int i = 0; i < 4; ++i) {
				newTwo->gl_Position[i] = delta * newVertices[2]->gl_Position[i] + (1 - delta) * newTwo->gl_Position[i];
				// newTwo->data[i] = delta * newVertices[2]->data[i] + (1 - delta) * newTwo->data[i];
			}
			goesIn[0] = newVertices[1];
			goesIn[1] = newVertices[2];
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
			goesIn[0] = newVertices[1];
			goesIn[1] = newOne;
			goesIn[2] = newTwo;
			clip_triangle(state, goesIn, face + 1);
		}
	}
	else {
		// std::cout << "		All are inside." << std::endl;
		for (int i = 0; i < 3; ++i) {
			goesIn[i] = newVertices[i];
		}
		clip_triangle(state, goesIn, face + 1);
	}




    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
	
}





// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
	
	int image_index = 0;
	data_geometry* nonhomogeneous = new data_geometry[3];
	for (int k = 0; k < 3; ++k) {
		// std::cout << std::endl << "in[" << k << "]: " << in[k]->gl_Position[0] << ", " << in[k]->gl_Position[1] << ", " << in[k]->gl_Position[2] << ", " << in[k]->gl_Position[3] << std::endl;
		nonhomogeneous[k].gl_Position = in[k]->gl_Position / in[k]->gl_Position[3];
		// std::cout << "nonhomogeneous[" << k << "]: " << nonhomogeneous[k].gl_Position[0] << ", " << nonhomogeneous[k].gl_Position[1] << ", " << nonhomogeneous[k].gl_Position[2] << ", " << nonhomogeneous[k].gl_Position[3] << std::endl;
		int i = state.image_width / 2.0  * nonhomogeneous[k].gl_Position[0] + (state.image_width - 1) / 2.0;
		int j = state.image_height / 2.0 * nonhomogeneous[k].gl_Position[1] + (state.image_height - 1) / 2.0;
		// std::cout << "i: " << i << "; j: " << j << "; ";
		image_index = i + j * state.image_width;
		// std::cout << "image_index: " << image_index << std::endl;
		state.image_color[image_index] = make_pixel(255, 255, 255);
	}
	// std::cout << std::endl;
	
	// Transformation of triangle vertices.
	double ax = state.image_width / 2.0 * nonhomogeneous[0].gl_Position[0] + (state.image_width - 1) / 2.0;
	double ay = state.image_height / 2.0 * nonhomogeneous[0].gl_Position[1] + (state.image_height - 1) / 2.0;
	double bx = state.image_width / 2.0 * nonhomogeneous[1].gl_Position[0] + (state.image_width - 1) / 2.0;
	double by = state.image_height / 2.0 * nonhomogeneous[1].gl_Position[1] + (state.image_height - 1) / 2.0;
	double cx = state.image_width / 2.0 * nonhomogeneous[2].gl_Position[0] + (state.image_width - 1) / 2.0;
	double cy = state.image_height / 2.0 * nonhomogeneous[2].gl_Position[1] + (state.image_height - 1) / 2.0;

	double triangleArea = 0.5 * ((bx * cy  - cx * by) - (ax * cy - cx * ay) + (ax * by - bx * ay));

	double minx = std::min(ax, bx);
	minx = std::min(minx, cx);
	double maxx = std::max(ax, bx);
	maxx = std::max(maxx, cx);
	double miny = std::min(ay, by);
	miny = std::min(miny, cy);
	double maxy = std::max(ay, by);
	maxy = std::max(maxy, cy);

	
	for (int k = std::floor(minx); k < std::ceil(maxx); ++k) {
		for (int l = std::floor(miny); l < std::ceil(maxy); ++l) {
			double pbcArea = 0.5 * ((bx * cy - cx * by) - (k * cy - cx * l) + (k * by - bx * l));
			double apcArea = 0.5 * ((k * cy - cx * l) - (ax * cy - cx * ay) + (ax * l - k * ay));
			double abpArea = 0.5 * ((bx * l - k * by) - (ax * l - k * ay) + (ax * by - bx * ay));
			double alpha = pbcArea / triangleArea;
			double beta = apcArea / triangleArea;
			double gamma = abpArea / triangleArea;
			// double gamma = 1 - alpha - beta;
			


int flatCounter = 0;
int smoothCounter = 0;
int noperspectiveCounter = 0;
int invalidCounter = 0;


			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				// std::cout << "pbcArea: " << pbcArea << "; apcArea: " << apcArea << "; abpArea: " << abpArea << "; alpha: " << alpha << "; beta: " << beta << "; gamma: " << gamma << "; ";
				data_fragment frag;
				float frag_data[MAX_FLOATS_PER_VERTEX];
				frag.data = frag_data;
				// std::cout << "k: " << k << "; l: " << l << "; ";
				image_index = k + l * state.image_width;
				// std::cout << "image_index: " << image_index << std::endl;
				for (int j = 0; j < state.floats_per_vertex; ++j) {
					switch (state.interp_rules[j]) {
						case interp_type::flat: {
							// based on first data_geometry, first vertex's data.
							// std::cout << "flat ";
							// std::cout << in[0]->data[j] << "  ";
							frag.data[j] = in[0]->data[j];
							++flatCounter;
							break;
						}
						case interp_type::smooth: {
							++smoothCounter;
							// std::cout << "smooth ";
							double alphaPrime = (alpha / in[0]->gl_Position[3]) / (alpha / in[0]->gl_Position[3] + beta / in[1]->gl_Position[3] + gamma / in[2]->gl_Position[3]);
							double betaPrime = (beta / in[1]->gl_Position[3]) / (alpha / in[0]->gl_Position[3] + beta / in[1]->gl_Position[3] + gamma / in[2]->gl_Position[3]);
							double gammaPrime = (gamma / in[2]->gl_Position[3]) / (alpha / in[0]->gl_Position[3] + beta / in[1]->gl_Position[3] + gamma / in[2]->gl_Position[3]);
							frag.data[j] = in[0]->data[j] * alphaPrime + in[1]->data[j] * betaPrime + in[2]->data[j] * gammaPrime;
							break;
						}
						case interp_type::noperspective: {
							// frag.data[0] = v[0].data[0] * alpha + v[1].data[0] * beta + v[2].data[0]
							// in is v!
							frag.data[j] = in[0]->data[j] * alpha + in[1]->data[j] * beta + in[2]->data[j] * gamma;
							++noperspectiveCounter;
							// std::cout << frag.data[j] << " ";
							break;
							// std::cout << "noperspective ";
						}
						case interp_type::invalid: {
							++invalidCounter;
							// std::cout << "invalid ";
							break;
						}
					}
				}
				data_output d_o;
				state.fragment_shader(frag, d_o, state.uniform_data);
				int r = 255 * d_o.output_color[0];
				int g = 255 * d_o.output_color[1];
				int b = 255 * d_o.output_color[2];
				// std::cout << d_o.output_color[0] << " " << d_o.output_color[1] << " " << d_o.output_color[2] << " ";
				// std::cout << r << " " << g << " " << b << "  ";
				
				vec4 v0 = in[0]->gl_Position / in[0]->gl_Position[3];
				vec4 v1 = in[1]->gl_Position / in[1]->gl_Position[3];
				vec4 v2 = in[2]->gl_Position / in[2]->gl_Position[3];

				double z = alpha * v0[2] + beta * v1[2] + gamma * v2[2];
				// std::cout << "z: " << z << " ";
				// std::cout << "state.image_depth[" << image_index << "]: " << state.image_depth[image_index] << " ";
				if (z < state.image_depth[image_index]) {
					state.image_depth[image_index] = z;
					state.image_color[image_index] = make_pixel(r, g, b);
				}
			}
			/* else {
				// std::cout << "      alpha: " << alpha << "; beta: " << beta << "; gamma: " << gamma << std::endl;
			} */
			// std::cout << flatCounter << " " << smoothCounter << " " << noperspectiveCounter << " " << invalidCounter << "   ";
		}
	}
	// std::cout << "-------------------------------------------------------------------------------" << std::endl;
	/* std::cout << "ax: " << ax << "; ay: " << ay << "; bx: " << bx << "; by: " << by << "; cx: " << cx << "; cy: " << cy << std::endl;
	 std::cout << triangleArea << std::endl;
	 std::cout << "minx: " << minx << "; maxx: " << maxx << std::endl;
	 std::cout << "miny: " << miny << "; maxy: " << maxy << std::endl; */
}