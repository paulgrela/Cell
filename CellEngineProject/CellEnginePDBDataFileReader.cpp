
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineConfigData.h"
#include "CellEnginePDBDataFileReader.h"
#include "CellEngineParticlesKindsManager.h"

using namespace std;
using namespace string_utils;

CellEngineAtom CellEnginePDBDataFileReader::ParseRecord(const char* LocalPDBRecord)
{
    CellEngineAtom CellEngineAtomObject{};

    try
    {
        string RecordStr = LocalPDBRecord;

        CellEngineAtomObject.EntityId = 1;

        string NameStr = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        string ResNameStr = trim_whitespace_surrounding(RecordStr.substr(17, 3));
        strncpy(CellEngineAtomObject.Name, NameStr.c_str(), NameStr.length() + 1);
        strncpy(CellEngineAtomObject.ResName, ResNameStr.c_str(), ResNameStr.length() + 1);

        CellEngineAtomObject.X = stof(RecordStr.substr(30, 8));
        CellEngineAtomObject.Y = stof(RecordStr.substr(38, 8));
        CellEngineAtomObject.Z = stof(RecordStr.substr(46, 8));
        
        auto AtomKindObjectIterator = ParticlesKindsManagerObject.GetGraphicAtomKindDataFromAtomName(CellEngineAtomObject.Name[0]);
        CellEngineAtomObject.AtomColor = AtomKindObjectIterator->Color;

        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXAtom = AtomKindObjectIterator->SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObjectIterator->SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObjectIterator->SizeZ;

        CellEngineAtomObject.SizeXParticle = 1;
        CellEngineAtomObject.SizeYParticle = 1;
        CellEngineAtomObject.SizeZParticle = 1;
        #endif
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

void CellEnginePDBDataFileReader::ReadDataFromFile(const bool StartValuesBool, const bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        GetMemoryForParticlesInSectors();

        string Line;
        vector<CellEngineAtom> ListOfAtoms;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(CellEngineConfigDataObject.CellStateFileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING OF PDB FILE"));

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == "ATOM" || Line.substr(0, 6) == "HETATM")
                ListOfAtoms.emplace_back(ParseRecord(Line.c_str()));
            else
            if (Line.substr(0, 3) == "END" )
            {
                vmath::vec3 CenterOfParticle(0.0, 0.0, 0.0);
                for (const CellEngineAtom& AtomObject : ListOfAtoms)
                    CenterOfParticle += AtomObject.Position();
                vector3_Real32 Center = { CenterOfParticle.X() / static_cast<RealType>(ListOfAtoms.size()), CenterOfParticle.Y() / static_cast<RealType>(ListOfAtoms.size()), CenterOfParticle.Z() / static_cast<RealType>(ListOfAtoms.size()) };

                Particle ParticleToInsert(1, 1, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                ParticleToInsert.ListOfAtoms = ListOfAtoms;
                ParticleToInsert.Center = Center;
                ParticleToInsert.ParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor());
                ParticleToInsert.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor());
                ParticleToInsert.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor());

                GetParticleFromIndex(1) = ParticleToInsert;

                EntityIdInt P = 1;
                ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::Lipid, P, "LipidOuterMembrane", "LipidOuterMembrane", "LipidOuterMembrane", -1, 0, "m", 1);
                ParticlesKindsManagerObject.GetParticleKind(1).GraphicData.Visible = true;
            }
        }

        LoggersManagerObject.Log(STREAM("FINISHED READING OF PDB FILE"));

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Reading data from pdb file has taken time: ","executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from PDB file")
}