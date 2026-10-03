#ifndef FILE_UTILS_H_
#define FILE_UTILS_H_

#include <string>

namespace FileUtils
{
	void RewriteTextToFile(const char* FileName, std::string_view TextToAppend);
	void AppendTextToFile(const char* FileName, std::string_view TextToAppend);
	std::string ReadFirstLineFromFile(const char* FileName);
};

#endif
