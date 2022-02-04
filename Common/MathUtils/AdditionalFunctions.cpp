
#include "AdditionalFunctions.h"

double DegToRad(double deg)
{
	return M_PI * deg / 180.0f;
}

double sinDeg(double deg)
{
	return sin(DegToRad(deg));
}

double cosDeg(double deg)
{
	return cos(DegToRad(deg));
}

float sinDegf(float deg)
{
	return (float)sinDeg(deg);
}

float cosDegf(float deg)
{
	return (float)cosDeg(deg);
}

