
#include <fstream>

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "DestinationPlatform.h"
#include "CellEngineParticlesBinaryDataFileReaderWriter.h"

using namespace std;

void SavePointerToBinaryFile(ofstream& ParticlesDataFile, const Particle* PointerToParticle)
{
    try
    {
        if (PointerToParticle != nullptr)
            ParticlesDataFile.write((char*)&PointerToParticle->Index, sizeof(UniqueIdInt));
        else
        {
            UniqueIdInt ValueToWrite = 0;
            ParticlesDataFile.write((char*)&ValueToWrite, sizeof(UniqueIdInt));
        }
    }
    CATCH("saving pointer to binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::SaveParticlesKindsToBinaryFile(ofstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF SAVING PARTICLES KINDS TO BINARY FILE"));

        UnsignedInt ParticlesKindsSize = ParticlesKindsManagerObject.ParticlesKinds.size();
        ParticlesDataFile.write((char*)&ParticlesKindsSize, sizeof(ParticlesKindsSize));
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            ParticlesDataFile.write((char*)&ParticleKindObject.EntityId, sizeof(ParticleKindObject.EntityId));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.SizeX, sizeof(ParticleKindObject.GraphicData.SizeX));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.SizeY, sizeof(ParticleKindObject.GraphicData.SizeY));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.SizeZ, sizeof(ParticleKindObject.GraphicData.SizeZ));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.AtomColor, sizeof(ParticleKindObject.GraphicData.AtomColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.ParticleColor, sizeof(ParticleKindObject.GraphicData.ParticleColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.GraphicData.RandomParticleColor, sizeof(ParticleKindObject.GraphicData.RandomParticleColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.XSizeDiv2, sizeof(ParticleKindObject.XSizeDiv2));
            ParticlesDataFile.write((char*)&ParticleKindObject.YSizeDiv2, sizeof(ParticleKindObject.YSizeDiv2));
            ParticlesDataFile.write((char*)&ParticleKindObject.ZSizeDiv2, sizeof(ParticleKindObject.ZSizeDiv2));

            UniqueIdInt ParticleKindObjectNameLength = ParticleKindObject.Name.length();
            ParticlesDataFile.write((char*)&ParticleKindObjectNameLength, sizeof(ParticleKindObjectNameLength));
            ParticlesDataFile.write((char*)ParticleKindObject.Name.c_str(), ParticleKindObjectNameLength);

            UniqueIdInt ParticleKindObjectSymbolLength = ParticleKindObject.Symbol.length();
            ParticlesDataFile.write((char*)&ParticleKindObjectSymbolLength, sizeof(ParticleKindObjectSymbolLength));
            ParticlesDataFile.write((char*)ParticleKindObject.Symbol.c_str(), ParticleKindObjectSymbolLength);

            UniqueIdInt ParticleKindObjectNameFromXMLLength = ParticleKindObject.GraphicData.NameFromXML.length();
            ParticlesDataFile.write((char*)&ParticleKindObjectNameFromXMLLength, sizeof(ParticleKindObjectNameFromXMLLength));
            ParticlesDataFile.write((char*)ParticleKindObject.GraphicData.NameFromXML.c_str(), ParticleKindObjectNameFromXMLLength);

            UniqueIdInt ParticleKindObjectNameFromDataFileLength = ParticleKindObject.GraphicData.NameFromDataFile.length();
            ParticlesDataFile.write((char*)&ParticleKindObjectNameFromDataFileLength, sizeof(ParticleKindObjectNameFromDataFileLength));
            ParticlesDataFile.write((char*)ParticleKindObject.GraphicData.NameFromDataFile.c_str(), ParticleKindObjectNameFromDataFileLength);

            UnsignedInt ParticleKindListOfVoxelsSize = ParticleKindObject.ListOfVoxels.size();
            ParticlesDataFile.write((char*)&ParticleKindListOfVoxelsSize, sizeof(ParticleKindListOfVoxelsSize));
            for (const auto& VoxelObject : ParticleKindObject.ListOfVoxels)
                ParticlesDataFile.write((char*)&VoxelObject, sizeof(VoxelObject));
        }

        LoggersManagerObject.Log(STREAM("END OF SAVING PARTICLES KINDS TO BINARY FILE"));
    }
    CATCH("saving particles to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveParticlesToBinaryFile(ofstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF SAVING PARTICLES TO BINARY FILE"));

        UnsignedInt ParticlesSize = Particles->size();
        LoggersManagerObject.Log(STREAM("Number of Particles to be saved = " << ParticlesSize));
        ParticlesDataFile.write((char*)&ParticlesSize, sizeof(ParticlesSize));

        for (const auto& ParticleObject : *Particles)
        {
            ParticlesDataFile.write((char*)&ParticleObject.second.EntityId, sizeof(ParticleObject.second.EntityId));
            ParticlesDataFile.write((char*)&ParticleObject.second.ChainId, sizeof(ParticleObject.second.ChainId));
            ParticlesDataFile.write((char*)&ParticleObject.second.Index, sizeof(ParticleObject.second.Index));
            ParticlesDataFile.write((char*)&ParticleObject.second.GenomeIndex, sizeof(ParticleObject.second.GenomeIndex));
            ParticlesDataFile.write((char*)&ParticleObject.second.ElectricCharge, sizeof(ParticleObject.second.ElectricCharge));

            ParticlesDataFile.write((char*)&ParticleObject.second.Center, sizeof(ParticleObject.second.Center));

            ParticlesDataFile.write((char*)&ParticleObject.second.UniqueColor, sizeof(ParticleObject.second.UniqueColor));

            ParticlesDataFile.write((char*)&ParticleObject.second.SelectedForReaction, sizeof(ParticleObject.second.SelectedForReaction));

            UnsignedInt ParticleListOfVoxelsSize = ParticleObject.second.ListOfVoxels.size();
            ParticlesDataFile.write((char*)&ParticleListOfVoxelsSize, sizeof(ParticleListOfVoxelsSize));
            for (const auto& VoxelObject : ParticleObject.second.ListOfVoxels)
                ParticlesDataFile.write((char*)&VoxelObject, sizeof(VoxelObject));

            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Prev);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Next);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.PairedNucleotide);

            UnsignedInt LinkedParticlesPointersSize = ParticleObject.second.LinkedParticlesPointersList.size();
            ParticlesDataFile.write((char*)&LinkedParticlesPointersSize, sizeof(LinkedParticlesPointersSize));
            for (const auto& LinkedParticlePointerObject : ParticleObject.second.LinkedParticlesPointersList)
                SavePointerToBinaryFile(ParticlesDataFile, LinkedParticlePointerObject);
        }

        LoggersManagerObject.Log(STREAM("END OF SAVING PARTICLES TO BINARY FILE"));
    }
    CATCH("saving particles to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveParticlesKindsAndParticlesToBinaryFile()
{
    try
    {
        LoggersManagerObject.Log(STREAM("SAVING OF SAVING DATA TO BINARY FILE"));

        string ParticlesDataFileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("particles") + OS_DIR_SEP + string("ParticlesDataFile.dat");
        ofstream ParticlesDataFile(ParticlesDataFileName, ios_base::out | ios_base::trunc | ios_base::binary);

        SaveParticlesKindsToBinaryFile(ParticlesDataFile);
        SaveParticlesToBinaryFile(ParticlesDataFile);

        ParticlesDataFile.close();

        LoggersManagerObject.Log(STREAM("END OF SAVING DATA TO BINARY FILE"));
    }
    CATCH("saving particles to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveDataToFile()
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        SaveParticlesKindsAndParticlesToBinaryFile();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of saving data to file has taken time: ", "executing printing duration_time")));
    }
    CATCH("saving data to file")
};








void CellEngineParticlesBinaryDataFileReaderWriter::PrepareParticlesAfterReadingFromBinaryFile()
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF PREPARING PARTICLES"));

        for (auto& ParticleObject : *Particles)
        {
            CellEngineVoxelSimulationSpaceObjectPointer->SetAllVoxelsInListOfVoxelsToValueOut(ParticleObject.second.ListOfVoxels, ParticleObject.second.Index);

            ParticleObject.second.Prev = ParticleObject.second.PrevTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PrevTemporary) : nullptr;
            ParticleObject.second.Next = ParticleObject.second.NextTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.NextTemporary) : nullptr;
            ParticleObject.second.PairedNucleotide = ParticleObject.second.PairedNucleotideTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PairedNucleotideTemporary) : nullptr;

            ParticleObject.second.LinkedParticlesPointersList.clear();
            for (auto& LinkedParticlesPointerObjectTemporary : ParticleObject.second.LinkedParticlesPointersListTemporary)
                ParticleObject.second.LinkedParticlesPointersList.emplace_back(&GetParticleFromIndex(LinkedParticlesPointerObjectTemporary));
        }

        LoggersManagerObject.Log(STREAM("END OF PREPARING PARTICLES"));
    }
    CATCH("preparing particles after reading from file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesKindsFromBinaryFile(ifstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING PARTICLES KINDS FROM BINARY FILE"));

        UnsignedInt ParticlesKindsSize;
        ParticlesDataFile.read((char*)&ParticlesKindsSize, sizeof(ParticlesKindsSize));
        LoggersManagerObject.Log(STREAM("Number of Particles Kinds to be read = " << ParticlesKindsSize));

        for (UnsignedInt ParticleKindObjectIndex = 1; ParticleKindObjectIndex <= ParticlesKindsSize; ParticleKindObjectIndex++)
        {
            ParticleKind ParticleKindObject;

            ParticlesDataFile.read((char*)&ParticleKindObject.EntityId, sizeof(ParticleKindObject.EntityId));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeX, sizeof(ParticleKindObject.GraphicData.SizeX));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeY, sizeof(ParticleKindObject.GraphicData.SizeY));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeZ, sizeof(ParticleKindObject.GraphicData.SizeZ));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.AtomColor, sizeof(ParticleKindObject.GraphicData.AtomColor));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.ParticleColor, sizeof(ParticleKindObject.GraphicData.ParticleColor));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.RandomParticleColor, sizeof(ParticleKindObject.GraphicData.RandomParticleColor));
            ParticlesDataFile.read((char*)&ParticleKindObject.XSizeDiv2, sizeof(ParticleKindObject.XSizeDiv2));
            ParticlesDataFile.read((char*)&ParticleKindObject.YSizeDiv2, sizeof(ParticleKindObject.YSizeDiv2));
            ParticlesDataFile.read((char*)&ParticleKindObject.ZSizeDiv2, sizeof(ParticleKindObject.ZSizeDiv2));

            UniqueIdInt ParticleKindNameLength = 0;
            ParticlesDataFile.read((char*)&ParticleKindNameLength, sizeof(ParticleKindNameLength));
            ParticleKindObject.Name.resize(ParticleKindNameLength);
            ParticlesDataFile.read((char*)ParticleKindObject.Name.c_str(), ParticleKindNameLength);

            UniqueIdInt ParticleKindSymbolLength = 0;
            ParticlesDataFile.read((char*)&ParticleKindSymbolLength, sizeof(ParticleKindSymbolLength));
            ParticleKindObject.Symbol.resize(ParticleKindSymbolLength);
            ParticlesDataFile.read((char*)ParticleKindObject.Symbol.c_str(), ParticleKindSymbolLength);

            UniqueIdInt ParticleKindNameFromXMLLength = 0;
            ParticlesDataFile.read((char*)&ParticleKindNameFromXMLLength, sizeof(ParticleKindNameFromXMLLength));
            ParticleKindObject.GraphicData.NameFromXML.resize(ParticleKindNameFromXMLLength);
            ParticlesDataFile.read((char*)ParticleKindObject.GraphicData.NameFromXML.c_str(), ParticleKindNameFromXMLLength);

            UniqueIdInt ParticleKindNameFromDataFileLength = 0;
            ParticlesDataFile.read((char*)&ParticleKindNameFromDataFileLength, sizeof(ParticleKindNameFromDataFileLength));
            ParticleKindObject.GraphicData.NameFromDataFile.resize(ParticleKindNameFromDataFileLength);
            ParticlesDataFile.read((char*)ParticleKindObject.GraphicData.NameFromDataFile.c_str(), ParticleKindNameFromDataFileLength);

            UnsignedInt ParticleKindListOfVoxelsSize;
            ParticlesDataFile.read((char*)&ParticleKindListOfVoxelsSize, sizeof(ParticleKindListOfVoxelsSize));
            ParticleKindObject.ListOfVoxels.clear();
            for (UnsignedInt VoxelObjectIndex = 1; VoxelObjectIndex <= ParticleKindListOfVoxelsSize; VoxelObjectIndex++)
            {
                vector3_16 VoxelObject;
                ParticlesDataFile.read((char*)&VoxelObject, sizeof(VoxelObject));
                ParticleKindObject.ListOfVoxels.push_back(VoxelObject);
            }
            ParticleKindObject.GraphicData.Visible = true;

            ParticlesKindsManagerObject.ParticlesKinds.emplace_back(ParticleKindObject);
            ParticlesKindsManagerObject.ParticlesKindsPos[ParticleKindObject.EntityId] = ParticlesKindsManagerObject.ParticlesKinds.size() - 1;
        }

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES KINDS FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}


void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesFromBinaryFile(ifstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING PARTICLES FROM BINARY FILE"));

        UnsignedInt ParticlesSize;
        ParticlesDataFile.read((char*)&ParticlesSize, sizeof(ParticlesSize));
        LoggersManagerObject.Log(STREAM("Number of Particles to be read = " << ParticlesSize));

        for (UnsignedInt ParticleObjectIndex = 1; ParticleObjectIndex <= ParticlesSize; ParticleObjectIndex++)
        {
            Particle ParticleObject;

            ParticlesDataFile.read((char*)&ParticleObject.EntityId, sizeof(ParticleObject.EntityId));
            ParticlesDataFile.read((char*)&ParticleObject.ChainId, sizeof(ParticleObject.ChainId));
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
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesKindsAndParticlesFromBinaryFile()
{
    try
    {
        Particles->clear();
        ParticlesKindsManagerObject.ParticlesKinds.clear();
        ParticlesKindsManagerObject.ParticlesKindsPos.clear();

        LoggersManagerObject.Log(STREAM("START OF READING DATA FROM BINARY FILE"));

        ifstream ParticlesDataFile(CellEngineConfigDataObject.CellStateFileName, ios_base::in | ios_base::binary);

        ReadParticlesKindsFromBinaryFile(ParticlesDataFile);
        ReadParticlesFromBinaryFile(ParticlesDataFile);

        ParticlesDataFile.close();

        LoggersManagerObject.Log(STREAM("END OF READING DATA FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesFromBinaryFileAndPrepareData(const bool StartValuesBool)
{
    try
    {
        if (StartValuesBool == true)
            SetStartValues();

        LoggersManagerObject.Log(STREAM("START OF READING PARTICLES FROM FILE AND PREPARING DATA"));

        ReadParticlesKindsAndParticlesFromBinaryFile();

        PrepareParticlesAfterReadingFromBinaryFile();

        //CellEngineVoxelSimulationSpaceObjectPointer->PreprocessData(false);
        PreprocessData();

        CellEngineConfigDataObject.GenomeReadFromFile = true;

        LoggersManagerObject.Log(STREAM("END OF PREPROCESSING DATA"));

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM FILE AND PREPARING DATA"));
    }
    CATCH("reading particles from file")
};
