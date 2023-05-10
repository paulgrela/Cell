
#include <regex>
#include <fstream>

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineUseful.h"
#include "CellEngineColors.h"
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

        CellEngineAtomObject.SetAtomPositionsData(stof(AtomFields[10]), stof(AtomFields[11]), stof(AtomFields[12]));

        auto AtomKindObjectIterator = CellEngineConfigDataObject.GetAtomKindDataForAtom(AtomFields[3][0]);
        CellEngineAtomObject.AtomColor = AtomKindObjectIterator->Color;
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXAtom = AtomKindObjectIterator->SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObjectIterator->SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObjectIterator->SizeZ;
        #endif

        auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(CellEngineAtomObject.EntityId)->second];
        CellEngineAtomObject.Visible = ParticleKindObject.Visible;
        CellEngineAtomObject.ParticleColor = ParticleKindObject.AtomColor;
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXParticle = ParticleKindObject.SizeX;
        CellEngineAtomObject.SizeYParticle = ParticleKindObject.SizeY;
        CellEngineAtomObject.SizeZParticle = ParticleKindObject.SizeZ;
        #endif
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

void CellEngineCIFDataFile::ReadDataFromFile()
{
    try
    {
        SetStartValues();

        ChainsNames.clear();
        TransformationsMatrixes.clear();

        string Line;
        std::vector<CellEngineAtom> LocalCellEngineAllAtomsObject;
        std::vector<CellEngineAtom> LocalCellEngineParticlesCentersObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(CellEngineConfigDataObject.CellStateFileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING FROM CIF FILE"));

        UnsignedInt NumberOfAtoms = 0;
        UnsignedInt NumberOfAtomsDNA = 0;
        UnsignedInt NumberOfNucleotidesInDNA = 0;
        UnsignedInt NumberOfParticles = 0;

        CellEngineColorsObject.SelectRandomEngineForColors();

        smatch SMatchObject;

        vector<UnsignedInt> AppliedMatrixesIds;
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
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(GraphicParticleKind{static_cast<UnsignedInt>(stoi(AtomFields[2])), OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.ParticleColor, OthersParticleKindObjectIterator->second.ParticleColor,CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });
                }
                else
                    CellEngineConfigDataObject.ParticlesKinds.emplace_back(GraphicParticleKind{static_cast<UnsignedInt>(stoi(AtomFields[2])), ParticleKindObjectIterator->second.Visible, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.ParticleColor, ParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });
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

                UnsignedInt TransformationMatrix3x4ObjectId = stoi(MatrixFields[0]);

                UnsignedInt Shift = 6;
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
                AppliedMatrixesIds.clear();
                auto pos = Line.cbegin();
                auto end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject1); pos = SMatchObject.suffix().first)
                    if (SMatchObject.str(1).empty() == false && SMatchObject.str(2).empty() == false)
                        for (UnsignedInt MatrixId = stoi(SMatchObject.str(1)); MatrixId <= stoi(SMatchObject.str(2)); MatrixId++)
                            AppliedMatrixesIds.emplace_back(MatrixId);
                    else
                        AppliedMatrixesIds.emplace_back(stoi(SMatchObject.str(3)));

                AppliedChainsNames.clear();
                pos = Line.cbegin();
                end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject2); pos = SMatchObject.suffix().first)
                    AppliedChainsNames.emplace_back(SMatchObject.str(1));

                vmath::vec3 ChainColor = CellEngineColorsObject.GetRandomColor();

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

                            vmath::vec3 UniqueParticleColor = CellEngineColorsObject.GetRandomColor();

                            NumberOfParticles++;
                            UniqueIdInt ParticleIndex = AddNewParticle(NumberOfParticles, Particle(NumberOfParticles, AtomsForChainNameIterator->second[0].EntityId, CellEngineUseful::GetChainIdFromChainName(AppliedChainName), 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(UniqueParticleColor)));

                            if (CellEngineUseful::IsNucleotide(AppliedChainName))
                                NumberOfNucleotidesInDNA++;

                            for (auto AppliedAtom : AtomsForChainNameIterator->second)
                            {
                                NumberOfAtoms++;

                                if (CellEngineUseful::IsNucleotide(AppliedChainName))
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

                                AppliedAtom.SetAtomPositionsData(Result[0], Result[1], Result[2]);

                                AppliedAtom.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(UniqueParticleColor);
                                AppliedAtom.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(ChainColor);

                                InsertAtom(LocalCellEngineAllAtomsObject, AppliedAtom, ParticleIndex);
                            }
                        }
                    }

                    InsertGroupOfAtoms(LocalCellEngineParticlesCentersObject, LocalCellEngineAllAtomsObject);
                }
            }
        }

        InsertParticlesCenters(LocalCellEngineParticlesCentersObject);

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintStatistics();

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " NumberOfParticles = " << NumberOfParticles << " | LocalCellEngineParticlesCentersObject.size() = " << LocalCellEngineParticlesCentersObject.size() << " | AllAtoms.size() = " << AllAtoms.size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | NumberOfNucleotidesInDNA = " << NumberOfNucleotidesInDNA <<" | AtomsPositionsMatrixes.size() = " << TransformationsMatrixes.size() << " | " << endl));

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ", "executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}