
#include <regex>
#include <fstream>

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "DestinationPlatform.h"

#include "CellEngineUseful.h"
#include "CellEngineColors.h"
#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFileReader.h"

using namespace std;
using namespace string_utils;

CellEngineAtom CellEngineCIFDataFileReader::ParseRecord(const char* LocalCIFRecord)
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

        auto AtomKindObjectIterator = ParticlesKindsManagerObject.GetGraphicAtomKindDataFromAtomName(AtomFields[3][0]);
        CellEngineAtomObject.AtomColor = AtomKindObjectIterator->Color;
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXAtom = AtomKindObjectIterator->SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObjectIterator->SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObjectIterator->SizeZ;
        #endif

        auto ParticleKindObject = ParticlesKindsManagerObject.GetGraphicParticleKind(CellEngineAtomObject.EntityId);
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

void CopyFieldsToMatrix(TransformationMatrix3x4& TransformationMatrix3x4Object, vector<string>& MatrixFields)
{
    try
    {
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
    }
    CATCH("copying fields to matrix")
}

glm::vec3 CountResultPositionsFromTransformationMatrix(std::unordered_map<UnsignedInt, TransformationMatrix3x4>::iterator& TransformationMatrixIterator, CellEngineAtom& AppliedAtom)
{
    glm::vec3 Result;

    try
    {
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

        Result = TransformationMatrix * glm::vec3(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z) + glm::vec3(TransformationMatrixIterator->second.Matrix[0][3], TransformationMatrixIterator->second.Matrix[1][3], TransformationMatrixIterator->second.Matrix[2][3]);
    }
    CATCH("copying from matrix to matrix")

    return Result;
}

void CellEngineCIFDataFileReader::ReadDataFromFile()
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

                auto ParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(stoi(AtomFields[2]));
                if (ParticleKindObjectIterator == ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.end())
                {
                    auto OthersParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(10000);
                    ParticlesKindsManagerObject.ParticlesKinds.emplace_back(ParticleKindGraphicData{static_cast<EntityIdInt>(stoi(AtomFields[2])), OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.ParticleColor, OthersParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), OthersParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });
                }
                else
                    ParticlesKindsManagerObject.ParticlesKinds.emplace_back(ParticleKindGraphicData{static_cast<EntityIdInt>(stoi(AtomFields[2])), ParticleKindObjectIterator->second.Visible, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.ParticleColor, ParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectIterator->second.NameFromXML, AtomFields[5].substr(1, AtomFields[5].length() - 2) });

                ParticlesKindsManagerObject.ParticlesKindsPos[stoi(AtomFields[2])] = ParticlesKindsManagerObject.ParticlesKinds.size() - 1;
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

                CopyFieldsToMatrix(TransformationMatrix3x4Object, MatrixFields);

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
                            UniqueIdInt ParticleIndex = AddNewParticle(Particle(NumberOfParticles, AtomsForChainNameIterator->second[0].EntityId, CellEngineUseful::GetChainIdFromChainName(AppliedChainName), -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(UniqueParticleColor)));

                            if (CellEngineUseful::IsNucleotide(AppliedChainName))
                                NumberOfNucleotidesInDNA++;

                            for (auto AppliedAtom : AtomsForChainNameIterator->second)
                            {
                                NumberOfAtoms++;

                                if (CellEngineUseful::IsNucleotide(AppliedChainName))
                                    NumberOfAtomsDNA++;

                                auto TransformationMatrixIterator = TransformationsMatrixes.find(AppliedMatrixId);

                                auto Result = CountResultPositionsFromTransformationMatrix(TransformationMatrixIterator, AppliedAtom);

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

        PreprocessData();

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintStatistics();

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " NumberOfParticles = " << NumberOfParticles << " | LocalCellEngineParticlesCentersObject.size() = " << LocalCellEngineParticlesCentersObject.size() << " | AllAtoms.size() = " << AllAtoms.size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | NumberOfNucleotidesInDNA = " << NumberOfNucleotidesInDNA <<" | AtomsPositionsMatrixes.size() = " << TransformationsMatrixes.size() << " | " << endl));

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ", "executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void AddParticlesKinds()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 12, "Polymerase", "POL", 0, 0 });

        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", 0, 0 });

        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", "?", 0, 0 });

        LoggersManagerObject.Log(STREAM("ADDED PARTICLES KINDS"));
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineParticlesDataFileReader::ReadParticlesFromFile()
{
    try
    {
        try
        {
            SetStartValues();

            LoggersManagerObject.Log(STREAM("START OF READING PARTICLES FROM BINARY FILE"));

            string ParticlesDataFileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("particles") + OS_DIR_SEP + string("ParticlesDataFile.dat");
            ifstream ParticlesDataFile(ParticlesDataFileName, ios_base::in | ios_base::binary);
            //ifstream ParticlesDataFile(CellEngineConfigDataObject.CellStateFileName, ios_base::in | ios_base::binary);

            UnsignedInt ParticlesSize;
            ParticlesDataFile.read((char*)&ParticlesSize, sizeof(ParticlesSize));
            LoggersManagerObject.Log(STREAM("Number of Particles to be read = " << ParticlesSize));

            for (UnsignedInt ParticleObjectIndex = 1; ParticleObjectIndex <= ParticlesSize; ParticleObjectIndex++)
            {
                Particle ParticleObject;

                ParticlesDataFile.read((char*)&ParticleObject.EntityId, sizeof(ParticleObject.EntityId));
                ParticlesDataFile.read((char*)&ParticleObject.ChainId, sizeof(ParticleObject.ChainId));


                auto ParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(ParticleObject.EntityId);
                if (ParticleKindObjectIterator == ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.end())
                {
                    auto OthersParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(10000);
                    ParticlesKindsManagerObject.ParticlesKinds.emplace_back(ParticleKindGraphicData{static_cast<EntityIdInt>(ParticleObject.EntityId), OthersParticleKindObjectIterator->second.Visible, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.ParticleColor, OthersParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), OthersParticleKindObjectIterator->second.NameFromXML, "NAME" });
                }
                else
                    ParticlesKindsManagerObject.ParticlesKinds.emplace_back(ParticleKindGraphicData{static_cast<EntityIdInt>(ParticleObject.EntityId), ParticleKindObjectIterator->second.Visible, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.ParticleColor, ParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectIterator->second.NameFromXML, "NAME" });

                ParticlesKindsManagerObject.ParticlesKindsPos[ParticleObject.EntityId] = ParticlesKindsManagerObject.ParticlesKinds.size() - 1;



                ParticlesDataFile.read((char*)&ParticleObject.Index, sizeof(ParticleObject.Index));
                ParticlesDataFile.read((char*)&ParticleObject.GenomeIndex, sizeof(ParticleObject.GenomeIndex));
                ParticlesDataFile.read((char*)&ParticleObject.ElectricCharge, sizeof(ParticleObject.ElectricCharge));

                ParticlesDataFile.read((char*)&ParticleObject.Center, sizeof(ParticleObject.Center));

                ParticlesDataFile.read((char*)&ParticleObject.UniqueColor, sizeof(ParticleObject.UniqueColor));

                ParticlesDataFile.read((char*)&ParticleObject.SelectedForReaction, sizeof(ParticleObject.SelectedForReaction));

                UnsignedInt ParticleListOfVoxelsSize;
                ParticlesDataFile.read((char*)&ParticleListOfVoxelsSize, sizeof(ParticleListOfVoxelsSize));
                ParticleObject.ListOfVoxels.clear();
                for (UnsignedInt VoxelObjectIndex = 1; VoxelObjectIndex <= ParticleListOfVoxelsSize; VoxelObjectIndex++)
                {
                    vector3_16 VoxelObject;
                    ParticlesDataFile.read((char*)&VoxelObject, sizeof(VoxelObject));
                    ParticleObject.ListOfVoxels.push_back(VoxelObject);
                }

                ParticlesDataFile.read((char*)&ParticleObject.PrevTemporary, sizeof(ParticleObject.PrevTemporary));
                ParticlesDataFile.read((char*)&ParticleObject.NextTemporary, sizeof(ParticleObject.NextTemporary));
                ParticlesDataFile.read((char*)&ParticleObject.PairedNucleotideTemporary, sizeof(ParticleObject.PairedNucleotideTemporary));

                UnsignedInt LinkedParticlesPointersSize;
                ParticlesDataFile.read((char*)&LinkedParticlesPointersSize, sizeof(LinkedParticlesPointersSize));
                for (UnsignedInt LinkedParticlesPointerIndex = 1; LinkedParticlesPointerIndex <= LinkedParticlesPointersSize; LinkedParticlesPointerIndex++)
                {
                    UniqueIdInt LinkedParticlesPointerObject;
                    ParticlesDataFile.read((char*)&LinkedParticlesPointerObject, sizeof(LinkedParticlesPointerObject));
                    ParticleObject.LinkedParticlesPointersListTemporary.push_back(LinkedParticlesPointerObject);
                }

                AddNewParticle(ParticleObject);
            }

            ParticlesDataFile.close();

            LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM BINARY FILE"));

            PreprocessData();

            LoggersManagerObject.Log(STREAM("END OF PREPROCESSING DATA"));
        }
        CATCH("reading particles from file")
    }
    CATCH("reading data from particles file")
}

//void CellEngineParticlesDataFileReader::PrepareParticlesAfterReadingFromFile()
//{
//    try
//    {
//        LoggersManagerObject.Log(STREAM("START OF PREPARING PARTICLES"));
//
//        for (auto& ParticleObject : CellEngineVoxelSimulationSpaceObjectPointer->Particles)
//        {
//            for (const auto& VoxelObject : ParticleObject.second.ListOfVoxels)
//                GetSpaceVoxel(VoxelObject.X, VoxelObject.Y, VoxelObject.Z) = ParticleObject.second.Index;
//
//            ParticleObject.second.Prev = ParticleObject.second.PrevTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PrevTemporary) : nullptr;
//            ParticleObject.second.Next = ParticleObject.second.NextTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.NextTemporary) : nullptr;
//            ParticleObject.second.PairedNucleotide = ParticleObject.second.PairedNucleotideTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PairedNucleotideTemporary) : nullptr;
//
//            ParticleObject.second.LinkedParticlesPointersList.clear();
//            for (auto& LinkedParticlesPointerObjectTemporary : ParticleObject.second.LinkedParticlesPointersListTemporary)
//                ParticleObject.second.LinkedParticlesPointersList.emplace_back(&GetParticleFromIndex(LinkedParticlesPointerObjectTemporary));
//        }
//
//        LoggersManagerObject.Log(STREAM("END OF PREPARING PARTICLES"));
//    }
//    CATCH("preparing particles after reading from file")
//};

void CellEngineParticlesDataFileReader::ReadParticlesFromFileAndPrepareData()
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING PARTICLES FROM FILE AND PREPARING DATA"));

        AddParticlesKinds();

        CellEngineVoxelSimulationSpaceObjectPointer->Particles.clear();

        ReadParticlesFromFile();

        //CellEngineVoxelSimulationSpaceObjectPointer->SetValueToVoxelsForSelectedSpace(nullptr, 0, 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension);

        CellEngineVoxelSimulationSpaceObjectPointer->PrepareParticlesAfterReadingFromFile();

        CellEngineVoxelSimulationSpaceObjectPointer->InitiateFreeParticleIndexes();

        CellEngineConfigDataObject.GenomeReadFromFile = true;

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM FILE AND PREPARING DATA"));
    }
    CATCH("reading particles from file")
};

void CellEngineParticlesDataFileReader::ReadDataFromFile()
{
    SetStartValues();
    AddParticlesKinds();
    CellEngineVoxelSimulationSpaceObjectPointer->AddChemicalReactions();

    ReadParticlesFromFileAndPrepareData();
}