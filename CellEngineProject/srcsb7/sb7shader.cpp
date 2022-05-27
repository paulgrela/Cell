
#define _CRT_SECURE_NO_WARNINGS 1

#include "GL/gl3w.h"

#include <cstdio>

namespace sb7
{
    namespace shader
    {
        extern GLuint Load(const char * FileName, GLenum ShaderType, bool CheckErrors)
        {
            GLuint Result = 0;
            FILE * ShaderFile;
            size_t FileSize;
            char* Data;

            ShaderFile = fopen(FileName, "rb");

            if (!ShaderFile)
                return 0;

            fseek(ShaderFile, 0, SEEK_END);
            FileSize = ftell(ShaderFile);
            fseek(ShaderFile, 0, SEEK_SET);

            Data = new char [FileSize + 1];

            if (!Data)
                goto FailDataAllocLabel;

            fread(Data, 1, FileSize, ShaderFile);
            Data[FileSize] = 0;
            fclose(ShaderFile);

            Result = glCreateShader(ShaderType);

            if (!Result)
                goto FailShaderAlloc;

            glShaderSource(Result, 1, &Data, nullptr);

            delete [] Data;

            glCompileShader(Result);

            if (CheckErrors)
            {
                GLint status = 0;
                glGetShaderiv(Result, GL_COMPILE_STATUS, &status);

                if (!status)
                {
                    char buffer[4096];
                    glGetShaderInfoLog(Result, 4096, nullptr, buffer);
                    #ifdef _WIN32
                    OutputDebugStringA(FileName);
                    OutputDebugStringA(":");
                    OutputDebugStringA(buffer);
                    OutputDebugStringA("\n");
                    #else
                    fprintf(stderr, "%s: %s\n", FileName, buffer);
                    #endif
                    goto FailCompileShader;
                }
            }

            return Result;

        FailCompileShader:
            glDeleteShader(Result);

        FailShaderAlloc:;
        FailDataAllocLabel:

            return Result;
        }

        GLuint FromString(const char* Source, GLenum shader_type, bool CheckErrors)
        {
            GLuint Shader;

            Shader = glCreateShader(shader_type);

            const char * Strings[] = { Source };
            glShaderSource(Shader, 1, Strings, nullptr);

            glCompileShader(Shader);

            if (CheckErrors)
            {
                GLint Status = 0;
                glGetShaderiv(Shader, GL_COMPILE_STATUS, &Status);

                if (!Status)
                {
                    char Buffer[4096];
                    glGetShaderInfoLog(Shader, 4096, nullptr, Buffer);
                    #ifdef _WIN32
                    OutputDebugStringA(Buffer);
                    OutputDebugStringA("\n");
                    #else
                    fprintf(stderr, "%s\n", Buffer);
                    #endif
                    goto FailCompileShaderLabel;
                }
            }

            return Shader;

        FailCompileShaderLabel:
            glDeleteShader(Shader);

            return 0;
        }

        }

        namespace program
        {
            GLuint LinkFromShaders(const GLuint* Shaders, int ShaderCount, bool DeleteShaders, bool CheckErrors)
            {
                int ShaderIndex;

                GLuint Program;

                Program = glCreateProgram();

                for (ShaderIndex = 0; ShaderIndex < ShaderCount; ShaderIndex++)
                    glAttachShader(Program, Shaders[ShaderIndex]);

                glLinkProgram(Program);

                if (CheckErrors)
                {
                    GLint Status;
                    glGetProgramiv(Program, GL_LINK_STATUS, &Status);

                    if (!Status)
                    {
                        char Buffer[4096];
                        glGetProgramInfoLog(Program, 4096, nullptr, Buffer);
                        #ifdef _WIN32
                        OutputDebugStringA(Buffer);
                        OutputDebugStringA("\n");
                        #endif
                        glDeleteProgram(Program);
                        return 0;
                    }
                }

                if (DeleteShaders)
                    for (ShaderIndex = 0; ShaderIndex < ShaderCount; ShaderIndex++)
                        glDeleteShader(Shaders[ShaderIndex]);

                return Program;
            }
    }
}
