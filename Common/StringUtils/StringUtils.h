#pragma once

#ifndef _STRING_UTILS_H_
#define _STRING_UTILS_H_

#include <string>

namespace string_utils
{
	std::string upper_case_str_transform(std::string s);
	std::string upper_case_str_for(std::string s);
	std::string align_str(std::string s, char c, uint64_t count);
	std::string trim_whitespace_surrounding(const std::string& s);
};

#endif
