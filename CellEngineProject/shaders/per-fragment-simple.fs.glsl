#version 410 core

layout (location = 0) out vec4 color;

in VS_OUT
{
	vec3 C;
} 
fs_in;

void main(void)
{
	color = vec4(fs_in.C, 1.0);
}
