#version 410 core

layout (location = 0) in vec4 position;
layout (location = 1) in vec3 normal;

layout (std140) uniform constants
{
    mat4 mv_matrix;
    mat4 proj_matrix;
	vec3 color;
};

out VS_OUT
{
	vec3 C;
} 
vs_out;

void main(void)
{
    vec4 P = mv_matrix * position;

    gl_Position = proj_matrix * P;
		
	vs_out.C = color;
}
