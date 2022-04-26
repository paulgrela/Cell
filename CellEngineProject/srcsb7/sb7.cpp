
#include <sb7.h>

sb7::application * sb7::application::app = 0;

#include <GL/glext.h>

#include <string.h>

GL3WglProc sb6GetProcAddress(const char * funcname)
{
    return gl3wGetProcAddress(funcname);
}

int sb6IsExtensionSupported(const char * extname)
{
    GLint numExtensions;
    GLint i;

    glGetIntegerv(GL_NUM_EXTENSIONS, &numExtensions);

    for (i = 0; i < numExtensions; i++)
    {
        const GLubyte * e = glGetStringi(GL_EXTENSIONS, i);

        if (!strcmp((const char *)e, extname))
            return 1;
    }

    return 0;
}

void APIENTRY sb7::application::debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, GLvoid* userParam)
{
    reinterpret_cast<application *>(userParam)->onDebugMessage(source, type, id, severity, length, message);
}