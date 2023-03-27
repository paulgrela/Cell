
#include <regex>
#include <random>
#include <fstream>

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFile.h"

using namespace std;
using namespace string_utils;

CellEngineAtom CellEngineCIFDataFile::ParseRecord(const char* LocalCIFRecord)
{
    CellEngineAtom CellEngineAtomObject{};

    try
    {
        string RecordStr = LocalCIFRecord;

        vector<string> AtomFields = split(RecordStr, " ");

        CellEngineAtomObject.Serial = stoi(AtomFields[1]);

        strncpy(CellEngineAtomObject.Name, AtomFields[3].c_str(), AtomFields[3].length() + 1);
        strncpy(CellEngineAtomObject.ResName, AtomFields[5].c_str(), AtomFields[5].length() + 1);
        strncpy(CellEngineAtomObject.Chain, AtomFields[6].c_str(), AtomFields[6].length() + 1);
        CellEngineAtomObject.EntityId = stoi(AtomFields[7]);

        CellEngineAtomObject.X = stof(AtomFields[10]);
        CellEngineAtomObject.Y = stof(AtomFields[11]);
        CellEngineAtomObject.Z = stof(AtomFields[12]);

        auto AtomKindObjectIterator = CellEngineConfigDataObject.GetAtomKindDataForAtom(CellEngineAtomObject.Name[0]);
        CellEngineAtomObject.AtomColor = GetVector3FormVMathVec3(AtomKindObjectIterator->Color);
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXAtom = AtomKindObjectIterator->SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObjectIterator->SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObjectIterator->SizeZ;
        #endif

        auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(CellEngineAtomObject.EntityId)->second];
        CellEngineAtomObject.Visible = ParticleKindObject.Visible;
        CellEngineAtomObject.ParticleColor = GetVector3FormVMathVec3(ParticleKindObject.Color);
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXParticle = ParticleKindObject.SizeX;
        CellEngineAtomObject.SizeYParticle = ParticleKindObject.SizeY;
        CellEngineAtomObject.SizeZParticle = ParticleKindObject.SizeZ;
        #endif
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

int64_t XMin = 10000, XMax = -10000, YMin = 10000, YMax = -10000, ZMin = 10000, ZMax = -10000, FirstAccessToVoxel = 0;

void CellEngineCIFDataFile::ReadDataFromFile()
{
    try
    {
        string Line;
        std::vector<CellEngineAtom> LocalCellEngineParticlesCentersObject;
        std::vector<CellEngineAtom> LocalCellEngineAllAtomsObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(CellEngineConfigDataObject.CellStateFileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING FROM CIF FILE"));

        UnsignedIntType NumberOfAtoms = 0;
        UnsignedIntType NumberOfAtomsDNA = 0;

        CellEngineConfigDataObject.SelectRandomEngineForColors();

        smatch SMatchObject;

        vector<UnsignedIntType> AppliedMatrixesIds;
        regex RegexObject1(R"((\d+)-(\d+)|\((\d+)\))");

        vector<string> AppliedChainsNames;
        regex RegexObject2("([A-z]+[0-9]*)");

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == ". . ")
            {
                vector<string> AtomFields = split(Line, " ");

                auto ParticleKindObjectIterator = CellEngineConfigDataObject.ParticlesKindsXML.find(stoi(AtomFields[2]));
                if (ParticleKindObjectIterator == CellEngineConfigDataObject.ParticlesKindsXML.end())
                {
                    auto OthersParticleKindObjectIterator = CellEngineConfigDataObject.ParticlesKindsXML.find(10000);
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(ParticleKind{ stoi(AtomFields[2]), OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.Color, OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });
                }
                else
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(ParticleKind{ stoi(AtomFields[2]), ParticleKindObjectIterator->second.Visible, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.Color, ParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });

                CellEngineConfigDataObject.ParticlesKindsPos[stoi(AtomFields[2])] = CellEngineConfigDataObject.ParticlesKinds.size() - 1;
            }
            else
            if (Line.substr(0, 4) == "ATOM")
            {
                CellEngineAtom CellEngineAtomObject = ParseRecord(Line.c_str());
                ChainsNames[CellEngineAtomObject.Chain].emplace_back(CellEngineAtomObject);
            }
            else
            if (Line.find("point symmetry operation") != std::string::npos)
            {
                vector<string> MatrixFields = split(Line, " ");

                TransformationMatrix3x4 TransformationMatrix3x4Object{};

                UnsignedIntType TransformationMatrix3x4ObjectId = stoi(MatrixFields[0]);

                UnsignedIntType Shift = 6;
                TransformationMatrix3x4Object.Matrix[0][0] = stof(MatrixFields[Shift + 0]);
                TransformationMatrix3x4Object.Matrix[0][1] = stof(MatrixFields[Shift + 1]);
                TransformationMatrix3x4Object.Matrix[0][2] = stof(MatrixFields[Shift + 2]);
                TransformationMatrix3x4Object.Matrix[0][3] = stof(MatrixFields[Shift + 3]);
                TransformationMatrix3x4Object.Matrix[1][0] = stof(MatrixFields[Shift + 4]);
                TransformationMatrix3x4Object.Matrix[1][1] = stof(MatrixFields[Shift + 5]);
                TransformationMatrix3x4Object.Matrix[1][2] = stof(MatrixFields[Shift + 6]);
                TransformationMatrix3x4Object.Matrix[1][3] = stof(MatrixFields[Shift + 7]);
                TransformationMatrix3x4Object.Matrix[2][0] = stof(MatrixFields[Shift + 8]);
                TransformationMatrix3x4Object.Matrix[2][1] = stof(MatrixFields[Shift + 9]);
                TransformationMatrix3x4Object.Matrix[2][2] = stof(MatrixFields[Shift + 10]);
                TransformationMatrix3x4Object.Matrix[2][3] = stof(MatrixFields[Shift + 11]);

                TransformationsMatrixes[TransformationMatrix3x4ObjectId] = TransformationMatrix3x4Object;
            }
            else
            if (Line.substr(0, 4) == "1 '(")
            {
                AppliedMatrixesIds.clear();
                auto pos = Line.cbegin();
                auto end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject1); pos = SMatchObject.suffix().first)
                    if (SMatchObject.str(1).empty() == false && SMatchObject.str(2).empty() == false)
                        for (UnsignedIntType MatrixId = stoi(SMatchObject.str(1)); MatrixId <= stoi(SMatchObject.str(2)); MatrixId++)
                            AppliedMatrixesIds.push_back(MatrixId);
                    else
                        AppliedMatrixesIds.push_back(stoi(SMatchObject.str(3)));

                AppliedChainsNames.clear();
                pos = Line.cbegin();
                end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject2); pos = SMatchObject.suffix().first)
                    AppliedChainsNames.push_back(SMatchObject.str(1));

                vmath::vec3 ChainColor = CellEngineConfigDataObject.GetRandomColor();

                for (const auto& AppliedMatrixId : AppliedMatrixesIds)
                {
                    LocalCellEngineAllAtomsObject.clear();

                    if (AppliedChainsNames.empty() == true)
                        LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE AppliedChainsNames EMPTY = " << to_string(AppliedMatrixId)));

                    for (const auto& AppliedChainName : AppliedChainsNames)
                    {
                        auto AtomsForChainNameIterator = ChainsNames.find(AppliedChainName);
                        if (AtomsForChainNameIterator == ChainsNames.end())
                            LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE LACKS AppliedChainName: " << AppliedChainName));
                        else
                        {
                            if (AtomsForChainNameIterator->second.empty() == true)
                                LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE - ChainsNames.find(AppliedChainName)->second.size() == 0 " << AppliedChainName));

                            for (auto AppliedAtom : AtomsForChainNameIterator->second)
                            {
                                NumberOfAtoms++;

                                if (AppliedChainName == "NU01" || AppliedChainName == "NU11" || AppliedChainName == "NU02" || AppliedChainName == "NU12" || AppliedChainName == "NU03" || AppliedChainName == "NU13" || AppliedChainName == "NU04" || AppliedChainName == "NU14")
                                    NumberOfAtomsDNA++;

                                auto TransformationMatrixIterator = TransformationsMatrixes.find(AppliedMatrixId);

                                auto TransformationMatrix = glm::mat3();

                                TransformationMatrix[0][0] = TransformationMatrixIterator->second.Matrix[0][0];
                                TransformationMatrix[0][1] = TransformationMatrixIterator->second.Matrix[1][0];
                                TransformationMatrix[0][2] = TransformationMatrixIterator->second.Matrix[2][0];

                                TransformationMatrix[1][0] = TransformationMatrixIterator->second.Matrix[0][1];
                                TransformationMatrix[1][1] = TransformationMatrixIterator->second.Matrix[1][1];
                                TransformationMatrix[1][2] = TransformationMatrixIterator->second.Matrix[2][1];

                                TransformationMatrix[2][0] = TransformationMatrixIterator->second.Matrix[0][2];
                                TransformationMatrix[2][1] = TransformationMatrixIterator->second.Matrix[1][2];
                                TransformationMatrix[2][2] = TransformationMatrixIterator->second.Matrix[2][2];

                                auto Result = TransformationMatrix * glm::vec3(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z) + glm::vec3(TransformationMatrixIterator->second.Matrix[0][3], TransformationMatrixIterator->second.Matrix[1][3], TransformationMatrixIterator->second.Matrix[2][3]);

                                AppliedAtom.X = Result[0];
                                AppliedAtom.Y = Result[1];
                                AppliedAtom.Z = Result[2];

                                AppliedAtom.RandomParticleColor = GetVector3FormVMathVec3(ChainColor);

                                if (CellEngineConfigDataObject.VoxelWorld == true)
                                {
                                    int64_t SpaceX = static_cast<int64_t>(round(AppliedAtom.X / 4.0)) + 512;
                                    int64_t SpaceY = static_cast<int64_t>(round(AppliedAtom.Y / 4.0)) + 512;
                                    int64_t SpaceZ = static_cast<int64_t>(round(AppliedAtom.Z / 4.0)) + 512;

                                    XMin = min(SpaceX, XMin);
                                    XMax = max(SpaceX, XMax);
                                    YMin = min(SpaceY, YMin);
                                    YMax = max(SpaceY, YMax);
                                    ZMin = min(SpaceZ, YMin);
                                    ZMax = max(SpaceZ, YMax);

                                    if (Space[SpaceX][SpaceY][SpaceZ] == 0)
                                    {
                                        FirstAccessToVoxel++;

                                        AppliedAtom.X = (static_cast<float>(SpaceX - 512)) * 4;
                                        AppliedAtom.Y = (static_cast<float>(SpaceY - 512)) * 4;
                                        AppliedAtom.Z = (static_cast<float>(SpaceZ - 512)) * 4;

                                        LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                                    }

                                    Space[SpaceX][SpaceY][SpaceZ] = AppliedAtom.EntityId;
                                }
                                else
                                    LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                            }
                        }
                    }

                    LocalCellEngineParticlesCentersObject.emplace_back(CellEngineAtom(LocalCellEngineAllAtomsObject.front().X, LocalCellEngineAllAtomsObject.front().Y, LocalCellEngineAllAtomsObject.front().Z, AllAtoms.size(), LocalCellEngineParticlesCentersObject.size(), LocalCellEngineAllAtomsObject.front().Name, LocalCellEngineAllAtomsObject.front().ResName, LocalCellEngineAllAtomsObject.front().Chain, LocalCellEngineAllAtomsObject.front().ParticleColor));

                    AllAtoms.emplace_back(LocalCellEngineAllAtomsObject);
                }
            }
        }

        cout << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMin) << " ][ Xmax = " << to_string(XMax) << " ][ Ymin = " << to_string(YMin) << " ][ Ymax = " << to_string(YMax) << " ][ Zmin = " << to_string(ZMin) << " ][ Zmax = " << to_string(XMax) << " ] " << endl;
        uint64_t MaxRepVoxels = 0;
        uint64_t SumOfFullVoxels = 0;
        for (uint64_t X = 0; X < 1024; X++)
            for (uint64_t Y = 0; Y < 1024; Y++)
                for (uint64_t Z = 0; Z < 1024; Z++)
                    if (Space[X][Y][Z] != 0)
                    {
                        SumOfFullVoxels++;
                        MaxRepVoxels = max(Space[X][Y][Z], MaxRepVoxels);
                    }
        cout << "SumOfFullVoxels = " << to_string(SumOfFullVoxels) << " MaxRepVoxels = " << to_string(MaxRepVoxels) << " FirstAccessToVoxel = " << to_string(FirstAccessToVoxel) << endl;


        ParticlesCenters.emplace_back(LocalCellEngineParticlesCentersObject);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " | LocalCellEngineParticlesCentersObject.size() = " << LocalCellEngineParticlesCentersObject.size() << " | AllAtoms.size() = " << AllAtoms.size() << " | AllAtoms.back().size() = " << AllAtoms.back().size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | AtomsPositionsMatrixes.size() = " << TransformationsMatrixes.size() << " | " << endl));

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ", "executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}