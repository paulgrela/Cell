
#include <regex>
#include <fstream>

#include "vmath.h"
#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineCIFDataFile.h"

using namespace std;
using namespace string_utils;

vmath::vec3 ChooseColor(const string& ChainName)
{
    vmath::vec3 ChosenColor;

    try
    {
        if (ChainName.substr(0, 3) == "BAF")
            ChosenColor = vmath::vec3(0.25f, 0.75f, 0.75f);
        else
        if(ChainName.substr(0, 3) == "BAE")
            ChosenColor = vmath::vec3(1.00f, 0.00f, 0.00f);
        else
        if(ChainName.substr(0, 3) == "ATP")
            ChosenColor = vmath::vec3(1.00f, 1.00f, 1.00f);
        else
        if(ChainName.substr(0, 3) == "BAR")
            ChosenColor = vmath::vec3(0.00f, 0.00f, 1.00f);
        else
        if(ChainName.substr(0, 1) == "A")
            ChosenColor = vmath::vec3(0.50f, 0.50f, 0.20f);
        else
            ChosenColor = vmath::vec3(0, 0, 0);
    }
    CATCH("chosing color for atom for cell visualization")

    return ChosenColor;
}

CellEngineAtom CellEngineCIFDataFile::ParseRecord(const char* LocalCIFRecord)
{
    CellEngineAtom CellEngineAtomObject{};

    try
    {
        string RecordStr = LocalCIFRecord;

        vector<string> AtomFields = split(RecordStr, " ");

        strncpy(CellEngineAtomObject.ResName, AtomFields[5].c_str(), AtomFields[5].length() + 1);
        strncpy(CellEngineAtomObject.Chain, AtomFields[6].c_str(), AtomFields[6].length() + 1);
        CellEngineAtomObject.EntityId = stoi(AtomFields[7]);

        CellEngineAtomObject.X = stof(AtomFields[10]);
        CellEngineAtomObject.Y = stof(AtomFields[11]);
        CellEngineAtomObject.Z = stof(AtomFields[12]);
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

CellEngineCIFDataFile::CellEngineCIFDataFile(const string_view FileName)
{
    DrawBondsBetweenParticlesCenters = false;
    DrawBondsBetweenAtoms = false;

    CheckAtomVisibility =  true;
    CutZ = 200;
    Distance = 1500;
    //Distance = 200;

    SizeX = 1;
    SizeY = 1;
    SizeZ = 1;

    SizeStep = 0.1;

    LoadOfAtomsStep = 100;

    ViewStep = 50;

    CameraXMoveStep = 100;
    CameraYMoveStep = 100;
    CameraZMoveStep = 100;

    ViewXMoveStep = 100;
    ViewYMoveStep = 100;
    ViewZMoveStep = 100;

    ShowDetailsInAtomScale = false;

    ChosenStructureIndex = 0;

    ReadDataFromFile(FileName);
}

void CellEngineCIFDataFile::ReadDataFromFile(const std::string_view FileName)
{
    try
    {
        string Line;
        std::vector<CellEngineAtom> LocalCellEngineParticlesCentersObject;
        std::vector<CellEngineAtom> LocalCellEngineAllAtomsObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(FileName).c_str(), std::ios_base::in);

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

                vmath::vec3 ChainColor;
                if (AppliedChainsNames.empty() == false)
                    ChainColor = ChooseColor(*AppliedChainsNames.begin());
                if (ChainColor.data[0] == 0.0 && ChainColor.data[1] == 0.0 && ChainColor.data[2] == 0.0)
                    ChainColor = vmath::vec3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);

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

                                if (AppliedChainName == "BAR0" || AppliedChainName == "NR0" || AppliedChainName == "NR1")
                                    NumberOfAtomsDNA++;

                                auto TransformationMatrixIterator = TransformationsMatrixes.find(AppliedMatrixId);

                                auto TransformationMatrix = vmath::mat4();

                                TransformationMatrix[0][0] = TransformationMatrixIterator->second.Matrix[0][0];
                                TransformationMatrix[0][1] = TransformationMatrixIterator->second.Matrix[0][1];
                                TransformationMatrix[0][2] = TransformationMatrixIterator->second.Matrix[0][2];

                                TransformationMatrix[1][0] = TransformationMatrixIterator->second.Matrix[1][0];
                                TransformationMatrix[1][1] = TransformationMatrixIterator->second.Matrix[1][1];
                                TransformationMatrix[1][2] = TransformationMatrixIterator->second.Matrix[1][2];

                                TransformationMatrix[2][0] = TransformationMatrixIterator->second.Matrix[2][0];
                                TransformationMatrix[2][1] = TransformationMatrixIterator->second.Matrix[2][1];
                                TransformationMatrix[2][2] = TransformationMatrixIterator->second.Matrix[2][2];

                                TransformationMatrix[3][0] = 0.0;
                                TransformationMatrix[3][1] = 0.0;
                                TransformationMatrix[3][2] = 0.0;
                                TransformationMatrix[3][3] = 1.0;

                                TransformationMatrix[0][3] = 0.0;
                                TransformationMatrix[1][3] = 0.0;
                                TransformationMatrix[2][3] = 0.0;

                                auto Result = vmath::vec4(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, 1) + vmath::vec4(TransformationMatrixIterator->second.Matrix[0][3], TransformationMatrixIterator->second.Matrix[1][3], TransformationMatrixIterator->second.Matrix[2][3], 1.0) * TransformationMatrix;

                                AppliedAtom.X = Result[0];
                                AppliedAtom.Y = Result[1];
                                AppliedAtom.Z = Result[2];

                                AppliedAtom.Color = ChainColor;

                                LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                            }
                        }
                    }

                    vmath::vec3 Center = GetCenter(LocalCellEngineAllAtomsObject);
                    LocalCellEngineParticlesCentersObject.emplace_back(CellEngineAtom(Center.X(), Center.Y(), Center.Z(), AllAtoms.size(), LocalCellEngineParticlesCentersObject.size(), LocalCellEngineAllAtomsObject.front().Name, LocalCellEngineAllAtomsObject.front().ResName, LocalCellEngineAllAtomsObject.front().Chain, LocalCellEngineAllAtomsObject.front().Color));

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