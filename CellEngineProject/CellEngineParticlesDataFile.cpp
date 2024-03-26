
#include "DateTimeUtils.h"
#include "CellEngineParticlesDataFile.h"

using namespace std;

void CellEngineParticlesDataFile::ReadDataFromFile(const bool StartValuesBool, const bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        switch (Type)
        {
            case CellEngineConfigData::TypesOfFileToRead::BinaryFile : ReadParticlesFromBinaryFileAndPrepareData(StartValuesBool, UpdateParticleKindListOfVoxelsBool, Type); break;
            case CellEngineConfigData::TypesOfFileToRead::CIFFile : ReadDataFromCIFFile(); break;
            default : break;
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from file has taken time: ", "executing printing duration_time")));
    }
    CATCH("reading data from file")
}
