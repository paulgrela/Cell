
#ifndef SB6EXT_H
#define SB6EXT_H

#include "GL/glext.h"

GL3WglProc sb6GetProcAddress(const char * FunctionName);
int sb6IsExtensionSupported(const char * ExtensionName);

#endif
