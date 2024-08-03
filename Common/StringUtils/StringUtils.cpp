
#include <cstdint>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "StringUtils.h"

using namespace std;

namespace string_utils
{
    bool check_end_str(const string_view FileName, const string_view FileExtension)
    {
        return std::equal(FileExtension.rbegin(), FileExtension.rend(), string(FileName).rbegin());
    }

	string upper_case_str_transform(string& s)
	{
		transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}

	string upper_case_str_for(string& s)
	{
		for (auto &c : s) 
			c = toupper(c);
		return s;
	}

	string align_str(const string& s, char c, uint64_t count)
	{
		stringstream ss;
		ss << setfill(c) << setw(count) << s;
		return ss.str();
	}

	string trim_whitespace_surrounding(const string& s)
	{
		const char whitespace[]{" \t\n"};
		const size_t first(s.find_first_not_of(whitespace));
		if (string::npos == first) 
			return {}; 
		const size_t last(s.find_last_not_of(whitespace));
		return s.substr(first, (last - first + 1));
	}

    vector<string> split(string& s, string_view delimiter)
    {
        size_t pos_start = 0;
        size_t pos_end;
        size_t delim_len = delimiter.length();
        string token;
        vector<string> result_list_of_words;

        while ((pos_end = s.find(delimiter, pos_start)) != string::npos)
        {
            token = s.substr(pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            result_list_of_words.push_back(token);
        }

        result_list_of_words.push_back(s.substr(pos_start));
        return result_list_of_words;
    }
}