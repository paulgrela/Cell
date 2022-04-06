
#include <regex>
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineCIFDataFile.h"

using namespace std;
using namespace string_utils;

Atom::Atom(const char* CIFRecord, UnsignedIntType AtomIndex)
{
    try
    {
        ParseRecord(CIFRecord);
        this->AtomIndex = AtomIndex;
    }
    CATCH("initation of Atom")
}

void Atom::ParseRecord(const char* LocalCIFRecord)
{
    try
    {
        string RecordStr = LocalCIFRecord;

        vector<string> AtomFields = split(RecordStr, " ");

        ResName = AtomFields[5];
        Chain = AtomFields[6];
        EntityId = stoi(AtomFields[7]);

        X = stod(AtomFields[10]);
        Y = stod(AtomFields[11]);
        Z = stod(AtomFields[12]);
    }
    CATCH("parsing atom record")
}

DoubleVectorType Atom::Position() const
{
    return DoubleVectorType(X, Y, Z);
}

CIFDataFile::CIFDataFile(const string_view FileName)
{
    ChosenStructureIndex = 0;

    ReadDataFromFile(FileName);
}

void CIFDataFile::ReadDataFromFile(const std::string_view FileName)
{
    try
    {
        string Line;
        std::vector<Atom> LocalAtomsObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(FileName).c_str(), std::ios_base::in);

        std::vector<string> ApplyEntities;
        std::vector<UnsignedIntType> MatrixesId;

        LoggersManagerObject.Log(STREAM("STARTED READING FROM CIF FILE"));

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == "ATOM")
            {
                Atom AtomObject(Line.c_str(), LocalAtomsObject.size());
                ChainsNames[AtomObject.Chain].push_back(AtomObject);
            }
            else
            if (Line.find("point symmetry operation") != std::string::npos)
            {
                vector<string> AtomFields = split(Line, " ");

                Matrix MatrixObject;

                MatrixObject.MatrixId = stoi(AtomFields[0]);

                UnsignedIntType Shift = 6;
                MatrixObject.Matrix[0][0] = stod(AtomFields[Shift + 0]);
                MatrixObject.Matrix[0][1] = stod(AtomFields[Shift + 1]);
                MatrixObject.Matrix[0][2] = stod(AtomFields[Shift + 2]);
                MatrixObject.Matrix[0][3] = stod(AtomFields[Shift + 3]);
                MatrixObject.Matrix[1][0] = stod(AtomFields[Shift + 4]);
                MatrixObject.Matrix[1][1] = stod(AtomFields[Shift + 5]);
                MatrixObject.Matrix[1][2] = stod(AtomFields[Shift + 6]);
                MatrixObject.Matrix[1][3] = stod(AtomFields[Shift + 7]);
                MatrixObject.Matrix[2][0] = stod(AtomFields[Shift + 8]);
                MatrixObject.Matrix[2][1] = stod(AtomFields[Shift + 9]);
                MatrixObject.Matrix[2][2] = stod(AtomFields[Shift + 10]);
                MatrixObject.Matrix[2][3] = stod(AtomFields[Shift + 11]);

                Matrixes.push_back(MatrixObject);

                //NAJPIERW NIECH PRZEPISZE do atomow tylko same matrixes
                //i dopiero jak sie zblizam to przepisuje dla najblizszych punktow dodatkowe atomy
            }
            else
            if (Line.substr(0, 4) == "1 '(")
            {
                vector<UnsignedIntType> AppliedMatrixesIds;
                smatch SMatchObject;
                regex RegexObject("(\\d+)-(\\d+)");
                auto pos = Line.cbegin();
                auto end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject); pos = SMatchObject.suffix().first)
                    for (UnsignedIntType MatrixId = stoi(SMatchObject.str(1)); MatrixId <= stoi(SMatchObject.str(2)); MatrixId++)
                        AppliedMatrixesIds.push_back(MatrixId);

                vector<string> AppliedChainsNames;
                regex RegexObject1("([A-Z]+\\d*)");
                pos = Line.cbegin();
                end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject1); pos = SMatchObject.suffix().first)
                    AppliedChainsNames.push_back(SMatchObject.str(1));

                for (const auto& AppliedMatrixId : AppliedMatrixesIds)
                    for (const auto& AppliedChainName : AppliedChainsNames)
                        for (const Atom& Atom : ChainsNames.find(AppliedChainName)->second)
                        {
                            Matrixes[AppliedMatrixId].Matrix[1][1] * Atom.X;
                            LocalAtomsObject.push_back(Atom);
                        }
            }
        }

        Atoms.push_back(LocalAtomsObject);

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ","executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}

IntType CIFDataFile::GetNumberOfAtoms() const
{
    return Atoms[ChosenStructureIndex].size();
}

const Atom& CIFDataFile::GetAtom(IntType DataRawIndex) const
{
    return Atoms[ChosenStructureIndex][DataRawIndex];
}

DoubleVectorType CIFDataFile::MassCenter() const
{
    DoubleVectorType MassCenter(0.0, 0.0, 0.0);

    try
    {
        for (const Atom& AtomObject : Atoms[ChosenStructureIndex])
            MassCenter += AtomObject.Position();

        MassCenter /= Atoms[ChosenStructureIndex].size();
    }
    CATCH_AND_THROW("counting mass center")

    return MassCenter;
}
