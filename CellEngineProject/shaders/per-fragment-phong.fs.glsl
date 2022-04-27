#version 410 core

layout (location = 0) out vec4 color;

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
	vec3 diffuse_albedo = fs_in.C;
	
    vec3 N = normalize(fs_in.N);
    vec3 L = normalize(fs_in.L);
    vec3 V = normalize(fs_in.V);

    vec3 R = reflect(-L, N);

    vec3 diffuse = max(dot(N, L), 0.0) * diffuse_albedo;
    
	vec3 specular = pow(max(dot(R, V), 0.0), specular_power) * specular_albedo;

    color = vec4(ambient + diffuse + specular, 1.0);
}
