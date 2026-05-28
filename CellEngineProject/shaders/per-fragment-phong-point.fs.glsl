#version 450 core

in flat uint instanceID;

layout (location = 0) out vec4 color;
layout (location = 1) out uint instanceIDOut;

in VS_OUT
{
    vec3 N;
    vec3 L;
    vec3 V;
    vec3 C;
}
fs_in;

uniform vec3 specular_albedo = vec3(0.7);
uniform float specular_power = 300.0;
uniform vec3 ambient = vec3(0.1, 0.1, 0.1);

void main(void)
{
    vec2 coord = gl_PointCoord * 2.0 - 1.0;
    float dist = length(coord);

    if (dist > 1.0)
        discard;

    vec3 normal;
    normal.xy = coord;
    normal.z = sqrt(max(0.0, 1.0 - dist * dist));

    vec3 L = normalize(vec3(0.3, 0.3, 1.0));
    float diffuse = max(dot(normal, L), 0.0);

    color = vec4(fs_in.C * (ambient + diffuse * vec3(0.8)), 1.0);

    instanceIDOut = instanceID;
}