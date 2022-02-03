

#include <sstream>
#include <iomanip>
#include <algorithm>

#include "StringUtils.h"

using namespace std;

namespace string_utils
{
	string upper_case_str_transform(string s)
	{
		transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}

	string upper_case_str_for(string s)
	{
		for (auto &c : s) 
			c = toupper(c);
		return s;
	}

	string align_str(string s, char c, uint64_t count)
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
}