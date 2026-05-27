#version 430 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

uniform int UseMainLineColor;
uniform vec3 MainLineColor;

uniform mat4 MVP;

out vec3 Color;

void main()
{
    gl_Position = MVP * vec4(aPos, 1.0);
    if (UseMainLineColor == 0)
        Color = aColor;
    else
        Color = MainLineColor;
}