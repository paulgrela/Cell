
#include <sb7.h>

sb7::OpenGLApplication* sb7::OpenGLApplication::OpenGLApplicationObject = 0;

#include <GL/glext.h>

#include <cstring>

GL3WglProc sb6GetProcAddress(const char* FunctionName)
{
    return gl3wGetProcAddress(FunctionName);
}

int sb6IsExtensionSupported(const char* ExtensionName)
{
    GLint numExtensions;

    glGetIntegerv(GL_NUM_EXTENSIONS, &numExtensions);

    for (GLint Extension = 0; Extension < numExtensions; Extension++)
    {
        const GLubyte* e = glGetStringi(GL_EXTENSIONS, Extension);

        if (!strcmp((const char*)e, ExtensionName))
            return 1;
    }

    return 0;
}

void APIENTRY sb7::OpenGLApplication::debug_callback(GLenum Source, GLenum Type, GLuint Id, GLenum Severity, GLsizei Length, const GLchar* Message, GLvoid* UserParam)
{
    reinterpret_cast<OpenGLApplication *>(UserParam)->OnDebugMessage(Source, Type, Id, Severity, Length, Message);
}