#version 450 core

layout (location = 0) in vec4 position;
layout (location = 1) in vec3 normal;

struct Particle
{
    mat4 mv_matrix;
    mat4 proj_matrix;
    vec3 color;
    float padding;
};

layout(std430, binding = 0) readonly buffer ParticleBuffer
{
    Particle particles[];
};

layout(std430, binding = 1) readonly buffer VisibleIndices
{
    uint VisibleList[];
};

out VS_OUT
{
    vec3 N;
    vec3 L;
    vec3 V;
	vec3 C;
} 
vs_out;

uniform vec3 light_pos = vec3(100.0, 100.0, 100.0);

void main(void)
{
//    uint ParticleIndex = VisibleList[gl_InstanceID];
//    Particle p = particles[ParticleIndex];

    Particle p = particles[gl_InstanceID];

    vec4 P = p.mv_matrix * position;

    vs_out.N = mat3(p.mv_matrix) * normal;

    vs_out.L = light_pos - P.xyz;

    vs_out.V = -P.xyz;

    gl_Position = p.proj_matrix * P;
	
	vs_out.C = p.color;
}
