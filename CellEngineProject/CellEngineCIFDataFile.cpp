
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
    ViewStep = 50;

    ChosenStructureIndex = 0;

    ReadDataFromFile(FileName);
}

void CellEngineCIFDataFile::ReadDataFromFile(const std::string_view FileName)
{
    try
    {
        string Line;
        std::vector<CellEngineAtom> LocalCellEngineAtomsObject;
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
        //regex RegexObject2("([A-Z]+\\d+)");
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
                vector<string> AtomFields = split(Line, " ");

                AtomsPositionMatrix3x4 Matrix3x4Object;

                UnsignedIntType Shift = 6;
                Matrix3x4Object.Matrix[0][0] = stof(AtomFields[Shift + 0]);
                Matrix3x4Object.Matrix[0][1] = stof(AtomFields[Shift + 1]);
                Matrix3x4Object.Matrix[0][2] = stof(AtomFields[Shift + 2]);
                Matrix3x4Object.Matrix[0][3] = stof(AtomFields[Shift + 3]);
                Matrix3x4Object.Matrix[1][0] = stof(AtomFields[Shift + 4]);
                Matrix3x4Object.Matrix[1][1] = stof(AtomFields[Shift + 5]);
                Matrix3x4Object.Matrix[1][2] = stof(AtomFields[Shift + 6]);
                Matrix3x4Object.Matrix[1][3] = stof(AtomFields[Shift + 7]);
                Matrix3x4Object.Matrix[2][0] = stof(AtomFields[Shift + 8]);
                Matrix3x4Object.Matrix[2][1] = stof(AtomFields[Shift + 9]);
                Matrix3x4Object.Matrix[2][2] = stof(AtomFields[Shift + 10]);
                Matrix3x4Object.Matrix[2][3] = stof(AtomFields[Shift + 11]);

                AtomsPositionsMatrixes.push_back(Matrix3x4Object);
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

                for (const auto& AppliedMatrixId : AppliedMatrixesIds)
                {
                    LocalCellEngineAllAtomsObject.clear();
                    bool First = true;

                    //AtomsPositionsMatrixes[AppliedMatrixId - 2].Used = true;

                    vmath::vec3 ChainColor;
                    if (AppliedChainsNames.empty() == false)
                        ChainColor = ChooseColor(*AppliedChainsNames.begin());
                    if (ChainColor.data[0] == 0.0 && ChainColor.data[1] == 0.0 && ChainColor.data[2] == 0.0)
                        ChainColor = vmath::vec3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);

                    if (AppliedChainsNames.empty() == true)
                        LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE AppliedChainsNames EMPTY = " << to_string(AppliedMatrixId)));

                    for (const auto& AppliedChainName : AppliedChainsNames)
                    {
                        auto CNI = ChainsNames.find(AppliedChainName);
                        if (CNI == ChainsNames.end())
                            LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE LACKS AppliedChainName: " << AppliedChainName));
                        else
                        {
                            //AtomsPositionsMatrixes[AppliedMatrixId - 2].Used = true;

                            if (ChainsNames.find(AppliedChainName)->second.empty() == true)
                                LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE - ChainsNames.find(AppliedChainName)->second.size() == 0 " << AppliedChainName));

                            for (CellEngineAtom& AppliedAtom : ChainsNames.find(AppliedChainName)->second)
                            {
                                NumberOfAtoms++;

                                if (AppliedChainName == "BAR0" || AppliedChainName == "NR0" || AppliedChainName == "NR1")
                                    NumberOfAtomsDNA++;





    //                            auto RotationMatrix = vmath::mat3();
    //
    //                            RotationMatrix.data[0][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0];
    //                            RotationMatrix.data[1][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1];
    //                            RotationMatrix.data[2][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2];
    //
    //                            RotationMatrix.data[0][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0];
    //                            RotationMatrix.data[1][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1];
    //                            RotationMatrix.data[2][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2];
    //
    //                            RotationMatrix.data[0][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0];
    //                            RotationMatrix.data[1][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1];
    //                            RotationMatrix.data[2][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2];
    //
    //                            auto Result = vmath::vec3(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z) * RotationMatrix + vmath::vec3(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3]);
    //
    //                            AppliedAtom.X = Result.data[0];
    //                            AppliedAtom.Y = Result.data[1];
    //                            AppliedAtom.Z = Result.data[2];


                                auto TransformationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);

                                TransformationMatrix.data[0][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0];
                                TransformationMatrix.data[1][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1];
                                TransformationMatrix.data[2][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2];
                                TransformationMatrix.data[3][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];

                                TransformationMatrix.data[0][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0];
                                TransformationMatrix.data[1][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1];
                                TransformationMatrix.data[2][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2];
                                TransformationMatrix.data[3][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];

                                TransformationMatrix.data[0][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0];
                                TransformationMatrix.data[1][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1];
                                TransformationMatrix.data[2][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2];
                                TransformationMatrix.data[3][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];

                                TransformationMatrix.data[0][3] = 1.0;
                                TransformationMatrix.data[1][3] = 1.0;
                                TransformationMatrix.data[2][3] = 1.0;
                                TransformationMatrix.data[3][3] = 1.0;

//                                TransformationMatrix.data[0][3] = 0.0;
//                                TransformationMatrix.data[1][3] = 0.0;
//                                TransformationMatrix.data[2][3] = 0.0;
//                                TransformationMatrix.data[3][3] = 1.0;

                                auto Result1 = vmath::vec4(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, 1) * TransformationMatrix;

                                AppliedAtom.X = Result1.data[0];
                                AppliedAtom.Y = Result1.data[1];
                                AppliedAtom.Z = Result1.data[2];

                                // * vmath::scale(vmath::vec3(1.5, 1.5, 1.5));

    //                            auto RotationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);
    //
    //                            RotationMatrix.data[0][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0];
    //                            RotationMatrix.data[1][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1];
    //                            RotationMatrix.data[2][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2];
    //                            RotationMatrix.data[0][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0];
    //                            RotationMatrix.data[1][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1];
    //                            RotationMatrix.data[2][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2];
    //                            RotationMatrix.data[0][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0];
    //                            RotationMatrix.data[1][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1];
    //                            RotationMatrix.data[2][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2];

    //                            RotationMatrix.data[0][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0];
    //                            RotationMatrix.data[1][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0];
    //                            RotationMatrix.data[2][0] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0];
    //                            RotationMatrix.data[0][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1];
    //                            RotationMatrix.data[1][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1];
    //                            RotationMatrix.data[2][1] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1];
    //                            RotationMatrix.data[0][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2];
    //                            RotationMatrix.data[1][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2];
    //                            RotationMatrix.data[2][2] = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2];

                                //auto TranslationMatrix = vmath::translate(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3]);
                                //auto Result = vmath::vec4(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, 1) * RotationMatrix * TranslationMatrix;
                                //auto Result = vmath::vec4(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, 1) * TranslationMatrix * RotationMatrix.transpose();



                                    bool aaa = false;

                                    if(aaa)
                                        LoggersManagerObject.Log(STREAM("MATRIX = [ " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2]) << " " << to_string(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3]) << " ]"));
                                    if(aaa)
                                        LoggersManagerObject.Log(STREAM("B1: " << to_string(AppliedAtom.X) << " " << to_string(AppliedAtom.Y) << " " << to_string(AppliedAtom.Z) << " " << AppliedAtom.Chain));

    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];

    //                                AppliedAtom.X = AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];
    //                                AppliedAtom.Y = AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];
    //                                AppliedAtom.Z = AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];

    //                            AppliedAtom.X = AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3] * AppliedAtom.X;
    //                            AppliedAtom.Y = AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3] * AppliedAtom.Y;
    //                            AppliedAtom.Z = AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3] * AppliedAtom.Z;

    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.Z;
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Z;
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z;
    //                            AppliedAtom.X = AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];
    //                            AppliedAtom.Y = AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];
    //                            AppliedAtom.Z = AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];
    //
    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.Z;
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Z;
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z;

                                if(aaa)
                                        LoggersManagerObject.Log(STREAM("B2: " << to_string(AppliedAtom.X) << " " << to_string(AppliedAtom.Y) << " " << to_string(AppliedAtom.Z) << " " << AppliedAtom.Chain));

    //                                AppliedAtom.X = AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];
    //                                AppliedAtom.Y = AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];
    //                                AppliedAtom.Z = AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];

    //                            AppliedAtom.X = AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3] * AppliedAtom.X;
    //                            AppliedAtom.Y = AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3] * AppliedAtom.Y;
    //                            AppliedAtom.Z = AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3] * AppliedAtom.Z;

    //                                AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.Z;
    //                                AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Z;
    //                                AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z;

                                if(aaa)
                                        LoggersManagerObject.Log(STREAM("B3: " << to_string(AppliedAtom.X) << " " << to_string(AppliedAtom.Y) << " " << to_string(AppliedAtom.Z) << " " << AppliedAtom.Chain << endl));
                                if(aaa)
                                        getchar();


    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3];
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3];
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3];

    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3] + AppliedAtom.X;
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3] + AppliedAtom.Y;
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3] + AppliedAtom.Z;


    //                            AppliedAtom.X = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[3][0];
    //                            AppliedAtom.Y = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[3][1];
    //                            AppliedAtom.Z = AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.X + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Y + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z + AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[3][2];

                                if(aaa)
                                        LoggersManagerObject.Log(STREAM("B4: " << to_string(AppliedAtom.X) << " " << to_string(AppliedAtom.Y) << " " << to_string(AppliedAtom.Z) << " " << AppliedAtom.Chain << endl));
                                if(aaa)
                                        getchar();
                                //AtomsPositionsMatrixes[AppliedMatrixId - 2].Used = true;


                                AppliedAtom.Color = ChainColor;

                                if (First == true)
                                {
                                    First = false;
                                    //AtomsPositionsMatrixes[AppliedMatrixId - 2].Used = true;
                                    //LocalCellEngineAtomsObject.emplace_back(CellEngineAtom(AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[0][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[1][3], AtomsPositionsMatrixes[AppliedMatrixId - 2].Matrix[2][3], AllAtoms.size() , LocalCellEngineAtomsObject.size(), (char*)("H"), (char*)("MTR"), (char*)AppliedChainName.c_str()));
                                    //LocalCellEngineAtomsObject.emplace_back(CellEngineAtom(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, LocalCellEngineAtomsObject.size(), LocalCellEngineAtomsObject.size(), (char*)("H"), (char*)("MTR"), (char*)AppliedChainName.c_str()));
                                    LocalCellEngineAtomsObject.emplace_back(CellEngineAtom(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z, AllAtoms.size(), LocalCellEngineAtomsObject.size(), (char*)("H"), (char*)("MTR"), (char*)AppliedChainName.c_str(), AppliedAtom.Color));
                                }

                                LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                            }
                        }
                    }
                    AllAtoms.emplace_back(LocalCellEngineAllAtomsObject);
                }
            }
        }

        LoggersManagerObject.Log(STREAM("Number of matrixes used = " << to_string(count_if(AtomsPositionsMatrixes.begin(), AtomsPositionsMatrixes.end(), [](AtomsPositionMatrix3x4& AtomsPositionMatrix3x4Object){ return AtomsPositionMatrix3x4Object.Used == 1; }))));

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " | LocalCellEngineAtomsObject.size() = " << LocalCellEngineAtomsObject.size() << " | AllAtoms.size() = " << AllAtoms.size() << " | AllAtoms.back().size() = " << AllAtoms.back().size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | AtomsPositionsMatrixes.size() = " << AtomsPositionsMatrixes.size() << " | ChainsNames[\"BAF0\"].size() = " << ChainsNames["BAF0"].size() << " | " << endl));

        Atoms.emplace_back(LocalCellEngineAtomsObject);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ", "executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}

IntType CellEngineCIFDataFile::GetNumberOfAtoms() const
{
    return Atoms[ChosenStructureIndex].size();
}

const CellEngineAtom& CellEngineCIFDataFile::GetAtom(IntType DataRawIndex) const
{
    return Atoms[ChosenStructureIndex][DataRawIndex];
}

FloatVectorType CellEngineCIFDataFile::MassCenter()
{
    FloatVectorType MassCenter(0.0, 0.0, 0.0);

    try
    {
        for (const CellEngineAtom& AtomObject : Atoms[ChosenStructureIndex])
            MassCenter += AtomObject.Position();

        MassCenter /= Atoms[ChosenStructureIndex].size();
    }
    CATCH_AND_THROW("counting mass center")

    return MassCenter;
}
