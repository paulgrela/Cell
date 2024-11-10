#pragma once

#ifndef _DATE_TIME_UTILS_H_
#define _DATE_TIME_UTILS_H_

#include <chrono>
#include <string>

std::string GetActualDateTimeStandardCPP(const char* sep1, const char* sep2, const char* sep3, const char* sep4, const char* sep5);
std::string GetActualDateTimeWindows(const char* sep1, const char* sep2, const char* sep3, const char* sep4, const char* sep5, const char* sep6);

std::string GetDurationTimeInOneLine(std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type duration_time, const char* TextToPrint, const char* ExceptionTextToPrint);
std::string GetDurationTimeInOneLineStr(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point stop_time, const char* TextToPrint, const char* ExceptionTextToPrint);

#endif