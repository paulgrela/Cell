
#include <regex>
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
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
        regex RegexObject2("([A-Z]+\\d*)");

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

                Matrix3x4Object.MatrixId = stoi(AtomFields[0]);

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

                AtomsPositiosnMatrixes.push_back(Matrix3x4Object);
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
                    for (const auto& AppliedChainName : AppliedChainsNames)
                    {
                        LocalCellEngineAllAtomsObject.clear();

                        bool First = true;
                        for (CellEngineAtom& AppliedAtom : ChainsNames.find(AppliedChainName)->second)
                        {
                            if (First == true)
                            {
                                First = false;
                                LocalCellEngineAtomsObject.emplace_back(CellEngineAtom(AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[0][3], AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[1][3], AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[2][3], LocalCellEngineAtomsObject.size(), LocalCellEngineAtomsObject.size(), (char*)("H"), (char*)("MTR"), (char*)AppliedChainName.c_str()));
                            }

                            NumberOfAtoms++;
                            if (AppliedChainName == "BAR0")
                                NumberOfAtomsDNA++;

                            AppliedAtom.X = AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[0][0] * AppliedAtom.X + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[0][1] * AppliedAtom.Y + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[0][2] * AppliedAtom.Z + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[0][3];
                            AppliedAtom.Y = AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[1][0] * AppliedAtom.X + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[1][1] * AppliedAtom.Y + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[1][2] * AppliedAtom.Z + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[1][3];
                            AppliedAtom.Z = AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[2][0] * AppliedAtom.X + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[2][1] * AppliedAtom.Y + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[2][2] * AppliedAtom.Z + AtomsPositiosnMatrixes[AppliedMatrixId - 2].Matrix[2][3];

                            LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
                        }

                        AllAtoms.emplace_back(LocalCellEngineAllAtomsObject);
                    }
            }
        }

        LoggersManagerObject.Log(STREAM(NumberOfAtoms << " " << LocalCellEngineAtomsObject.size() << " " << NumberOfAtomsDNA << " " << AtomsPositiosnMatrixes.size() << " " << ChainsNames["BAF0"].size() << endl));

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
