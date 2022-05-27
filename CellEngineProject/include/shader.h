
#ifndef __SHADER_H__
#define __SHADER_H__

namespace sb7
{
    namespace shader
    {
        GLuint Load(const char* FileName, GLenum ShaderType = GL_FRAGMENT_SHADER,
        #ifdef _DEBUG
                                bool CheckErrors = true);
        #else
                                bool CheckErrors = false);
        #endif

        GLuint FromString(const char* Source, GLenum ShaderType,
        #ifdef _DEBUG
                                bool CheckErrors = true);
        #else
                                bool CheckErrors = false);
        #endif
    }

    namespace program
    {

        GLuint LinkFromShaders(const GLuint* Shaders, int ShaderCount, bool DeleteShaders,
        #ifdef _DEBUG
                                bool CheckErrors = true);
        #else
                                bool CheckErrors = false);
        #endif
    }
}

#endif
