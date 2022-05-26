
#ifndef __SB6EXT_H__
#define __SB6EXT_H__

#include "GL/glext.h"

GL3WglProc sb6GetProcAddress(const char * FunctionName);
int sb6IsExtensionSupported(const char * ExtensionName);

#endif
