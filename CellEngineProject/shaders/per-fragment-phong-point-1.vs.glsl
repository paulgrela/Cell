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

    float viewDistance = -particleCenter.z;

    //if (viewDistance > billboardDistance + 1000 - 400)
    if (viewDistance > billboardDistance)
    {
        vec4 clipPos = ProjectionMatrix * vec4(particleCenter, 1.0);
        gl_Position = clipPos;

        //float radius = length(vec3(p.mv_matrix[0][0], p.mv_matrix[1][0], p.mv_matrix[2][0]));

        //vec4 offsetClip = ProjectionMatrix * vec4(particleCenter.x + radius, particleCenter.y, particleCenter.z, 1.0);

        //vec2 centerNDC = clipPos.xy / clipPos.w;
        //vec2 offsetNDC = offsetClip.xy / offsetClip.w;

        //float radiusNDC = length(offsetNDC - centerNDC);
        //float radiusPixels = radiusNDC * screenSize.x * 0.5;

        //gl_PointSize = clamp(radiusPixels, 1.0, 64.0);
        gl_PointSize = 2;

        vs_out.C = p.color;
        vs_out.N = vec3(0, 0, 1);
        vs_out.L = vec3(0, 0, 1);
        vs_out.V = vec3(0, 0, -1);
    }
    else
    {
        vec4 P = p.mv_matrix * position;
        vs_out.N = mat3(p.mv_matrix) * normal;
        vs_out.L = light_pos - P.xyz;
        vs_out.V = -P.xyz;
        vs_out.C = p.color;

        gl_Position = ProjectionMatrix * P;
    }
}