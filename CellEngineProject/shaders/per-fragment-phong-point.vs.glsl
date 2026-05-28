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

uniform float billboardDistance;
uniform vec2 screenSize;

void main(void)
{
    ParticleIn p = ParticlesIn[gl_InstanceID];
    instanceID = p.ParticleIdL;

    vec3 particleCenter = vec3(p.mv_matrix[3][0], p.mv_matrix[3][1], p.mv_matrix[3][2]);

    vec4 clipPos = ProjectionMatrix * vec4(particleCenter, 1.0);
    gl_Position = clipPos;

    gl_PointSize = 2;

    vs_out.C = p.color;
    vs_out.N = vec3(0, 0, 1);
    vs_out.L = vec3(0, 0, 1);
    vs_out.V = vec3(0, 0, -1);
}