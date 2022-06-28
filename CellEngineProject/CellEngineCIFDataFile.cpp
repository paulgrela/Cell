
#include <regex>
#include <fstream>

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
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

        auto AtomKindObjectIterator = CellEngineConfigDataObject.AtomsKinds.find(string(1, CellEngineAtomObject.Name[0]));
        if (AtomKindObjectIterator == CellEngineConfigDataObject.AtomsKinds.end())
            AtomKindObjectIterator = CellEngineConfigDataObject.AtomsKinds.find(string(1, 'E'));
        auto AtomKindObject = AtomKindObjectIterator->second;
        CellEngineAtomObject.AtomColor = AtomKindObject.Color;
        CellEngineAtomObject.SizeXAtom = AtomKindObject.SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObject.SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObject.SizeZ;

        //auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds.find(CellEngineAtomObject.EntityId)->second;
        auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(CellEngineAtomObject.EntityId)->second];
        CellEngineAtomObject.Visible = ParticleKindObject.Visible;
        CellEngineAtomObject.ParticleColor = ParticleKindObject.Color;
        CellEngineAtomObject.SizeXParticle = ParticleKindObject.SizeX;
        CellEngineAtomObject.SizeYParticle = ParticleKindObject.SizeY;
        CellEngineAtomObject.SizeZParticle = ParticleKindObject.SizeZ;
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

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

        smatch SMatchObject;

        vector<UnsignedIntType> AppliedMatrixesIds;
        regex RegexObject1("(\\d+)-(\\d+)");

        vector<string> AppliedChainsNames;
        regex RegexObject2("([A-z]+[0-9]*)");

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == ". . ")
            {
                vector<string> AtomFields = split(Line, " ");

                //TRZEBA WCZYTYWAC DO VECTORA i INDEKS W WEKTORZE ParticklesKindsPos - DZIEKI TEMU ULOZONE PO KOLEJI i TAK WSTAWIANE DO CHECKBOX_LIST
                //W CHECKBOX_LIST ZAZNACZENIE OZNACZA

                auto ParticleKindObjectIterator = CellEngineConfigDataObject.ParticlesKindsXML.find(stoi(AtomFields[2]));
                if (ParticleKindObjectIterator == CellEngineConfigDataObject.ParticlesKindsXML.end())
                {
                    auto OthersParticleKindObjectIterator = CellEngineConfigDataObject.ParticlesKindsXML.find(10000);
                    //CellEngineConfigDataObject.ParticlesKinds[stoi(AtomFields[2])] = ParticleKind{ OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.Color, OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) };
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(ParticleKind{ stoi(AtomFields[2]), OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.Color, OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });

//                    auto ParticleKindObjectIterator = CellEngineConfigDataObject.ParticlesKindsXML.find(stoi(AtomFields[2]));
//                    if (ParticleKindObjectIterator == CellEngineConfigDataObject.ParticlesKindsXML.end())
//                        ParticleKindObjectIterator->second.NameFromDataFile = AtomFields[5].substr(1, AtomFields[5].length() - 2);
                    CellEngineConfigDataObject.ParticlesKindsPos[stoi(AtomFields[2])] = CellEngineConfigDataObject.ParticlesKinds.size() - 1;
                }
                else
                {
                    //CellEngineConfigDataObject.ParticlesKinds[stoi(AtomFields[2])] = ParticleKind{ OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.Color, OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) };
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(ParticleKind{ stoi(AtomFields[2]), ParticleKindObjectIterator->second.Visible, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.Color, ParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });
                    CellEngineConfigDataObject.ParticlesKindsPos[stoi(AtomFields[2])] = CellEngineConfigDataObject.ParticlesKinds.size() - 1;
                }
                //else
//                if (ParticleKindObjectIterator == CellEngineConfigDataObject.ParticlesKindsXML.end())
//                    ParticleKindObjectIterator->second.NameFromDataFile = AtomFields[5].substr(1, AtomFields[5].length() - 2);
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

                TransformationMatrix3x4 TransofrmationMatrix3x4Object;

                uint64_t TransofrmationMatrix3x4ObjectId = stoi(MatrixFields[0]);

                UnsignedIntType Shift = 6;
                TransofrmationMatrix3x4Object.Matrix[0][0] = stof(MatrixFields[Shift + 0]);
                TransofrmationMatrix3x4Object.Matrix[0][1] = stof(MatrixFields[Shift + 1]);
                TransofrmationMatrix3x4Object.Matrix[0][2] = stof(MatrixFields[Shift + 2]);
                TransofrmationMatrix3x4Object.Matrix[0][3] = stof(MatrixFields[Shift + 3]);
                TransofrmationMatrix3x4Object.Matrix[1][0] = stof(MatrixFields[Shift + 4]);
                TransofrmationMatrix3x4Object.Matrix[1][1] = stof(MatrixFields[Shift + 5]);
                TransofrmationMatrix3x4Object.Matrix[1][2] = stof(MatrixFields[Shift + 6]);
                TransofrmationMatrix3x4Object.Matrix[1][3] = stof(MatrixFields[Shift + 7]);
                TransofrmationMatrix3x4Object.Matrix[2][0] = stof(MatrixFields[Shift + 8]);
                TransofrmationMatrix3x4Object.Matrix[2][1] = stof(MatrixFields[Shift + 9]);
                TransofrmationMatrix3x4Object.Matrix[2][2] = stof(MatrixFields[Shift + 10]);
                TransofrmationMatrix3x4Object.Matrix[2][3] = stof(MatrixFields[Shift + 11]);

                TransformationsMatrixes[TransofrmationMatrix3x4ObjectId] = TransofrmationMatrix3x4Object;
            }
            else
            if (Line.substr(0, 4) == "1 '(")
            {
                AppliedMatrixesIds.clear();
                auto pos = Line.cbegin();
                auto end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject1); pos = SMatchObject.suffix().first)
                    for (UnsignedIntType MatrixId = stoi(SMatchObject.str(1)); MatrixId <= stoi(SMatchObject.str(2)); MatrixId++)
                        AppliedMatrixesIds.push_back(MatrixId);

                AppliedChainsNames.clear();
                pos = Line.cbegin();
                end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject2); pos = SMatchObject.suffix().first)
                    AppliedChainsNames.push_back(SMatchObject.str(1));

                vmath::vec3 ChainColor = vmath::vec3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);

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

                                if (AppliedChainName == "BAR0" || AppliedChainName == "BAR1")
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

                                AppliedAtom.RandomParticleColor = ChainColor;

                                LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                            }
                        }
                    }

                    vmath::vec3 Center = GetCenter(LocalCellEngineAllAtomsObject);
                    LocalCellEngineParticlesCentersObject.emplace_back(CellEngineAtom(Center.X(), Center.Y(), Center.Z(), AllAtoms.size(), LocalCellEngineParticlesCentersObject.size(), LocalCellEngineAllAtomsObject.front().Name, LocalCellEngineAllAtomsObject.front().ResName, LocalCellEngineAllAtomsObject.front().Chain, LocalCellEngineAllAtomsObject.front().ParticleColor));

                    AllAtoms.emplace_back(LocalCellEngineAllAtomsObject);
                }
            }
        }

        ParticlesCenters.emplace_back(LocalCellEngineParticlesCentersObject);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " | LocalCellEngineParticlesCentersObject.size() = " << LocalCellEngineParticlesCentersObject.size() << " | AllAtoms.size() = " << AllAtoms.size() << " | AllAtoms.back().size() = " << AllAtoms.back().size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | AtomsPositionsMatrixes.size() = " << TransformationsMatrixes.size() << " | ChainsNames[\"BAF0\"].size() = " << ChainsNames["BAF0"].size() << " | " << endl));

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ", "executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}