
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

void APIENTRY sb7::OpenGLApplication::debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, GLvoid* userParam)
{
    reinterpret_cast<OpenGLApplication*>(userParam)->onDebugMessage(source, type, id, severity, length, message);
}