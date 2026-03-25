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
    bool isBillboard = (fs_in.V.z < -0.99);
                                                                                                                        //color = vec4(34.0,234.0,123.0,1.0);
    if (isBillboard && gl_PointCoord != vec2(0,0))
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

        //color = vec4(ambient + diffuse, 1.0);
        //color = vec4(34.0,234.0,123.0,1.0);

        instanceIDOut = instanceID;
    }
    else
    {
        //if (fs_in.C.x < 0.0)
        //{
        //    color = vec4(vec3(0.0, 1.0, 1.0), 1.0);
        //    instanceIDOut = 0u;
        //}
        //else
        {
            vec3 diffuse_albedo = fs_in.C;
            vec3 N = normalize(fs_in.N);
            vec3 L = normalize(fs_in.L);
            vec3 V = normalize(fs_in.V);
            vec3 R = reflect(-L, N);

            vec3 diffuse = max(dot(N, L), 0.0) * diffuse_albedo;
            vec3 specular = pow(max(dot(R, V), 0.0), specular_power) * specular_albedo;

            color = vec4(ambient + diffuse + specular, 1.0);
                                                                                                                        //color = vec4(34.0,234.0,123.0,1.0);
            instanceIDOut = instanceID;
        }
    }
}