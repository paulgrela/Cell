
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

void SaveStringToBinaryFile(ofstream& ParticlesDataFile, const string& StringToBeSaved)
{
    try
    {
        UniqueIdInt Length = StringToBeSaved.length();
        ParticlesDataFile.write((char*)&Length, sizeof(Length));
        ParticlesDataFile.write((char*)StringToBeSaved.c_str(), Length);
    }
    CATCH("saving string to binary file")
}

template <class TElement>
void SaveVectorToBinaryFile(ofstream& ParticlesDataFile, const vector<TElement>& VectorToBeSaved)
{
    try
    {
        UnsignedInt Size = VectorToBeSaved.size();
        ParticlesDataFile.write((char*)&Size, sizeof(Size));
        for (const auto& Object : VectorToBeSaved)
            if constexpr(std::is_pointer_v<TElement>)
                SavePointerToBinaryFile(ParticlesDataFile, Object);
            else
                ParticlesDataFile.write((char*)&Object, sizeof(Object));
    }
    CATCH("saving vector to binary file")
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

            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.Name);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.Symbol);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromXML);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromDataFile);

            SaveVectorToBinaryFile<vector3_16>(ParticlesDataFile, ParticleKindObject.ListOfVoxels);
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

            SaveVectorToBinaryFile<vector3_16>(ParticlesDataFile, ParticleObject.second.ListOfVoxels);

            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Prev);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Next);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.PairedNucleotide);

            SaveVectorToBinaryFile<Particle*>(ParticlesDataFile, ParticleObject.second.LinkedParticlesPointersList);
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

void ReadStringFromBinaryFile(ifstream& ParticlesDataFile, string& StringToRead)
{
    try
    {
        UniqueIdInt Length = 0;
        ParticlesDataFile.read((char*)&Length, sizeof(Length));
        StringToRead.resize(Length);
        ParticlesDataFile.read((char*)StringToRead.c_str(), Length);
    }
    CATCH("reading string from binary file")
}

template <class TElement>
void ReadVectorFromBinaryFile(ifstream& ParticlesDataFile, vector<TElement>& VectorToBeRead)
{
    try
    {
        ParticlesDataFile.clear();

        UnsignedInt Size;
        ParticlesDataFile.read((char*)&Size, sizeof(Size));
        for (UnsignedInt Index = 1; Index <= Size; Index++)
        {
            TElement Object{};
            ParticlesDataFile.read((char*)&Object, sizeof(Object));
            VectorToBeRead.emplace_back(Object);
        }
    }
    CATCH("reading vector from binary file")
}

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

            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.Name);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.Symbol);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromXML);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromDataFile);

            ReadVectorFromBinaryFile<vector3_16>(ParticlesDataFile, ParticleKindObject.ListOfVoxels);

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

            ReadVectorFromBinaryFile<vector3_16>(ParticlesDataFile, ParticleObject.ListOfVoxels);

            ParticlesDataFile.read((char*)&ParticleObject.PrevTemporary, sizeof(ParticleObject.PrevTemporary));
            ParticlesDataFile.read((char*)&ParticleObject.NextTemporary, sizeof(ParticleObject.NextTemporary));
            ParticlesDataFile.read((char*)&ParticleObject.PairedNucleotideTemporary, sizeof(ParticleObject.PairedNucleotideTemporary));

            ReadVectorFromBinaryFile<UniqueIdInt>(ParticlesDataFile, ParticleObject.LinkedParticlesPointersListTemporary);

            AddNewParticle(ParticleObject);
        }

        ParticlesDataFile.close();

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesKindsAndParticlesFromBinaryFile(CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        Particles->clear();
        ParticlesKindsManagerObject.ParticlesKinds.clear();
        ParticlesKindsManagerObject.ParticlesKindsPos.clear();

        LoggersManagerObject.Log(STREAM("START OF READING DATA FROM BINARY FILE"));

        string OpenFileName = (Type == CellEngineConfigData::TypesOfFileToRead::CIFFile ? CellEngineConfigDataObject.CellStateFileName : CellEngineConfigDataObject.CellStateFileNameBackup);

        ifstream ParticlesDataFile(OpenFileName, ios_base::in | ios_base::binary);

        ReadParticlesKindsFromBinaryFile(ParticlesDataFile);
        ReadParticlesFromBinaryFile(ParticlesDataFile);

        ParticlesDataFile.close();

        LoggersManagerObject.Log(STREAM("END OF READING DATA FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesFromBinaryFileAndPrepareData(const bool StartValuesBool, CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        if (StartValuesBool == true)
            SetStartValues();

        LoggersManagerObject.Log(STREAM("START OF READING PARTICLES FROM FILE AND PREPARING DATA"));

        ReadParticlesKindsAndParticlesFromBinaryFile(Type);

        PrepareParticlesAfterReadingFromBinaryFile();

        //CellEngineVoxelSimulationSpaceObjectPointer->PreprocessData(false);
        PreprocessData();

        CellEngineConfigDataObject.GenomeReadFromFile = true;

        LoggersManagerObject.Log(STREAM("END OF PREPROCESSING DATA"));

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM FILE AND PREPARING DATA"));
    }
    CATCH("reading particles from file")
};
