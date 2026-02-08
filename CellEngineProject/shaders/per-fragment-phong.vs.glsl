#version 450 core

layout (location = 0) in vec4 position;
layout (location = 1) in vec3 normal;

out flat uint instanceID;

struct ParticleIn
{
    mat4 mv_matrix;
    vec3 color;
    uint ParticleIdL;
};

layout(std430, binding = 2) readonly buffer ParticleBufferIn
{
    ParticleIn ParticlesIn[];
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
uniform mat4 ProjectionMatrix;

void main(void)
{
    ParticleIn p = ParticlesIn[gl_InstanceID];

    vec4 P = p.mv_matrix * position;

    vs_out.N = mat3(p.mv_matrix) * normal;

    vs_out.L = light_pos - P.xyz;

    vs_out.V = -P.xyz;

    gl_Position = ProjectionMatrix * P;

	vs_out.C = p.color;

	instanceID = p.ParticleIdL;
}
