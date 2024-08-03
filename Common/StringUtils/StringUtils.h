#pragma once

#ifndef _STRING_UTILS_H_
#define _STRING_UTILS_H_

#include <string>
#include <vector>
#include <cstdint>

namespace string_utils
{
    bool check_end_str(std::string_view FileName, std::string_view FileExtension);

	std::string upper_case_str_transform(std::string& s);
	std::string upper_case_str_for(std::string& s);
	std::string align_str(const std::string& s, char c, std::uint64_t count);
	std::string trim_whitespace_surrounding(const std::string& s);

    std::vector<std::string> split(std::string& s, std::string_view delimiter);
};

#endif
