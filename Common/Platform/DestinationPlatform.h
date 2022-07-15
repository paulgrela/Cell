
#ifndef DESTINATION_PLATFORM_H
#define DESTINATION_PLATFORM_H

#define UNIX_PLATFORM

#ifdef WINDOWS_PLATFORM
#define OS_DIR_SEP "\\"
#else
#define OS_DIR_SEP "/"
#endif

#endif
