
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Logger.h"
#include "ExceptionsMacro.h"
#include "DestinationPlatform.h"

using namespace std;

void ReadChemicalReactionsFromFile()
{
    try
    {
        namespace pt = boost::property_tree;

        boost::property_tree::ptree root;

        boost::property_tree::read_json(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("iMB155.json"), root);

        LoggersManagerObject.Log(STREAM("REACTIONS READ FROM FILE"));
    }
    CATCH("reading chemical reactions from file")
}