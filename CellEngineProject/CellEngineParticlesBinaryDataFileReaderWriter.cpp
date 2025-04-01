

#include <fstream>

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "DestinationPlatform.h"
#include "CellEngineAminoAcids.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineParticlesBinaryDataFileReaderWriter.h"

#include "CellEngineImGuiMenu.h"

using namespace std;

constexpr bool AdditionalInfoPrinting = false;
constexpr bool AddReactionsOfSuddenAppearanceAndDisappearance = false;

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
        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            ParticlesDataFile.write((char*)&ParticleKindObject.second.EntityId, sizeof(ParticleKindObject.second.EntityId));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GeneId, sizeof(ParticleKindObject.second.GeneId));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.ElectricCharge, sizeof(ParticleKindObject.second.ElectricCharge));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.Counter, sizeof(ParticleKindObject.second.Counter));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.SizeX, sizeof(ParticleKindObject.second.GraphicData.SizeX));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.SizeY, sizeof(ParticleKindObject.second.GraphicData.SizeY));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.SizeZ, sizeof(ParticleKindObject.second.GraphicData.SizeZ));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.AtomColor, sizeof(ParticleKindObject.second.GraphicData.AtomColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.ParticleColor, sizeof(ParticleKindObject.second.GraphicData.ParticleColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.GraphicData.RandomParticleColor, sizeof(ParticleKindObject.second.GraphicData.RandomParticleColor));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.XSizeDiv2, sizeof(ParticleKindObject.second.XSizeDiv2));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.YSizeDiv2, sizeof(ParticleKindObject.second.YSizeDiv2));
            ParticlesDataFile.write((char*)&ParticleKindObject.second.ZSizeDiv2, sizeof(ParticleKindObject.second.ZSizeDiv2));

            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.IdStr);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.Name);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.Formula);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.Compartment);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.GraphicData.NameFromXML);
            SaveStringToBinaryFile(ParticlesDataFile, ParticleKindObject.second.GraphicData.NameFromDataFile);

            UnsignedInt ParticleKindSpecialDataSectorSize = ParticleKindObject.second.ParticleKindSpecialDataSector.size();
            ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorSize, sizeof(ParticleKindSpecialDataSectorSize));
            for (const auto& ParticleKindSpecialDataSectorObject : ParticleKindObject.second.ParticleKindSpecialDataSector)
            {
                SaveStringToBinaryFile(ParticlesDataFile, ParticleKindSpecialDataSectorObject.Description);
                SaveStringToBinaryFile(ParticlesDataFile, ParticleKindSpecialDataSectorObject.AddedParticle);
                ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorObject.GeneId, sizeof(ParticleKindSpecialDataSectorObject.GeneId));
                ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorObject.IsProtein, sizeof(ParticleKindSpecialDataSectorObject.IsProtein));
                ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorObject.ParticleType, sizeof(ParticleKindSpecialDataSectorObject.ParticleType));
                ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorObject.CleanProductOfTranscription, sizeof(ParticleKindSpecialDataSectorObject.CleanProductOfTranscription));
                ParticlesDataFile.write((char*)&ParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation, sizeof(ParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation));
            }

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                SaveVectorToBinaryFile<vector3_16>(ParticlesDataFile, ParticleKindObject.second.ListOfVoxels);
            else
                SaveVectorToBinaryFile<CellEngineAtom>(ParticlesDataFile, ParticleKindObject.second.ListOfAtoms);
        }

        LoggersManagerObject.Log(STREAM("END OF SAVING PARTICLES KINDS TO BINARY FILE - SAVED = " << ParticlesKindsManagerObject.ParticlesKinds.size()));
    }
    CATCH("saving particles to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveParticlesToBinaryFile(ofstream& ParticlesDataFile) const
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF SAVING PARTICLES TO BINARY FILE"));

        UnsignedInt ParticlesSize = 0;
        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
            ParticlesSize++;
        LoggersManagerObject.Log(STREAM("Number of Particles to be saved = " << ParticlesSize));
        ParticlesDataFile.write((char*)&ParticlesSize, sizeof(ParticlesSize));

        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
        {
            ParticlesDataFile.write((char*)&ParticleObject.second.EntityId, sizeof(ParticleObject.second.EntityId));
            ParticlesDataFile.write((char*)&ParticleObject.second.ChainId, sizeof(ParticleObject.second.ChainId));
            ParticlesDataFile.write((char*)&ParticleObject.second.Index, sizeof(ParticleObject.second.Index));
            ParticlesDataFile.write((char*)&ParticleObject.second.GenomeIndex, sizeof(ParticleObject.second.GenomeIndex));
            ParticlesDataFile.write((char*)&ParticleObject.second.ElectricCharge, sizeof(ParticleObject.second.ElectricCharge));
            ParticlesDataFile.write((char*)&ParticleObject.second.Center, sizeof(ParticleObject.second.Center));
            ParticlesDataFile.write((char*)&ParticleObject.second.UniqueParticleColor, sizeof(ParticleObject.second.UniqueParticleColor));
            ParticlesDataFile.write((char*)&ParticleObject.second.SelectedForReaction, sizeof(ParticleObject.second.SelectedForReaction));

            SaveStringToBinaryFile(ParticlesDataFile, ParticleObject.second.SequenceStr);

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                SaveVectorToBinaryFile<vector3_16>(ParticlesDataFile, ParticleObject.second.ListOfVoxels);
            else
                SaveVectorToBinaryFile<CellEngineAtom>(ParticlesDataFile, ParticleObject.second.ListOfAtoms);

            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Prev);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.Next);
            SavePointerToBinaryFile(ParticlesDataFile, ParticleObject.second.PairedNucleotidePtr);

            SaveVectorToBinaryFile<Particle*>(ParticlesDataFile, ParticleObject.second.LinkedParticlesPointersList);
        }

        LoggersManagerObject.Log(STREAM("END OF SAVING PARTICLES TO BINARY FILE - SAVED = " << ParticlesSize));
    }
    CATCH("saving particles to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveChemicalReactionsToBinaryFile(ofstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF SAVING CHEMICAL REACTIONS TO BINARY FILE"));

        UnsignedInt ChemicalReactionsSize = ChemicalReactionsManagerObject.ChemicalReactions.size();
        LoggersManagerObject.Log(STREAM("Number of chemical reactions to be saved = " << ChemicalReactionsSize));
        ParticlesDataFile.write((char*)&ChemicalReactionsSize, sizeof(ChemicalReactionsSize));

        for (const auto& ChemicalReactionObject : ChemicalReactionsManagerObject.ChemicalReactions)
            if (ChemicalReactionObject.SpecialReactionFunction == nullptr)
            {
                ParticlesDataFile.write((char*)&ChemicalReactionObject.ReactionIdNum, sizeof(ChemicalReactionObject.ReactionIdNum));
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactionIdStr);
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactionName);
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactantsStr);
                ParticlesDataFile.write((char*)&ChemicalReactionObject.Duration, sizeof(ChemicalReactionObject.Duration));
                ParticlesDataFile.write((char*)&ChemicalReactionObject.Reversible, sizeof(ChemicalReactionObject.Reversible));
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.UpperFluxBound);
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.LowerFluxBound);
                ParticlesDataFile.write((char*)&ChemicalReactionObject.AdditionalParameter1, sizeof(ChemicalReactionObject.AdditionalParameter1));
                ParticlesDataFile.write((char*)&ChemicalReactionObject.AdditionalParameter2, sizeof(ChemicalReactionObject.AdditionalParameter2));
                SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionObject.Comment);

                UnsignedInt ChemicalReactionReactantsSize = ChemicalReactionObject.Reactants.size();
                ParticlesDataFile.write((char*)&ChemicalReactionReactantsSize, sizeof(ChemicalReactionReactantsSize));
                for (const auto& ChemicalReactionReactantObject : ChemicalReactionObject.Reactants)
                {
                    ParticlesDataFile.write((char*)&ChemicalReactionReactantObject.EntityId, sizeof(ChemicalReactionReactantObject.EntityId));
                    ParticlesDataFile.write((char*)&ChemicalReactionReactantObject.Counter, sizeof(ChemicalReactionReactantObject.Counter));
                    ParticlesDataFile.write((char*)&ChemicalReactionReactantObject.ToRemoveInReaction, sizeof(ChemicalReactionReactantObject.ToRemoveInReaction));
                    SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionReactantObject.SequenceStr);
                    SaveVectorToBinaryFile<ChainIdInt>(ParticlesDataFile, ChemicalReactionReactantObject.Sequence);
                    SaveVectorToBinaryFile<UniqueIdInt>(ParticlesDataFile, ChemicalReactionReactantObject.LinkedParticleTypes);
                }

                UnsignedInt ChemicalReactionProductsSize = ChemicalReactionObject.Products.size();
                ParticlesDataFile.write((char*)&ChemicalReactionProductsSize, sizeof(ChemicalReactionProductsSize));
                for (const auto& ChemicalReactionProductObject : ChemicalReactionObject.Products)
                {
                    ParticlesDataFile.write((char*)&ChemicalReactionProductObject.EntityId, sizeof(ChemicalReactionProductObject.EntityId));
                    ParticlesDataFile.write((char*)&ChemicalReactionProductObject.Counter, sizeof(ChemicalReactionProductObject.Counter));
                    ParticlesDataFile.write((char*)&ChemicalReactionProductObject.ToRemoveInReaction, sizeof(ChemicalReactionProductObject.ToRemoveInReaction));
                    SaveStringToBinaryFile(ParticlesDataFile, ChemicalReactionProductObject.SequenceStr);
                    SaveVectorToBinaryFile<ChainIdInt>(ParticlesDataFile, ChemicalReactionProductObject.Sequence);
                    SaveVectorToBinaryFile<UniqueIdInt>(ParticlesDataFile, ChemicalReactionProductObject.LinkedParticleTypes);
                }
            }

        LoggersManagerObject.Log(STREAM("END OF SAVING CHEMICAL REACTIONS TO BINARY FILE"));
    }
    CATCH("saving chemical reactions to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveGenesToBinaryFile(ofstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF SAVING GENES TO BINARY FILE"));

        UnsignedInt GenesSize = ParticlesKindsManagerObject.Genes.size();
        LoggersManagerObject.Log(STREAM("Number of genes to be saved = " << GenesSize));
        ParticlesDataFile.write((char*)&GenesSize, sizeof(GenesSize));

        for (const auto& GeneObject : ParticlesKindsManagerObject.Genes)
        {
            ParticlesDataFile.write((char*)&GeneObject.second.NumId, sizeof(GeneObject.second.NumId));
            SaveStringToBinaryFile(ParticlesDataFile, GeneObject.second.StrId);
            SaveStringToBinaryFile(ParticlesDataFile, GeneObject.second.Description);
            SaveStringToBinaryFile(ParticlesDataFile, GeneObject.second.ProteinId);
            SaveStringToBinaryFile(ParticlesDataFile, GeneObject.second.Sequence);
            ParticlesDataFile.write((char*)&GeneObject.second.StartPosInGenome, sizeof(GeneObject.second.StartPosInGenome));
            ParticlesDataFile.write((char*)&GeneObject.second.EndPosInGenome, sizeof(GeneObject.second.EndPosInGenome));
        }

        LoggersManagerObject.Log(STREAM("END OF SAVING GENES TO BINARY FILE"));
    }
    CATCH("saving chemical genes to file")
};

void CellEngineParticlesBinaryDataFileReaderWriter::SaveParticlesKindsAndParticlesAndChemicalReactionsAndGenesToBinaryFile() const
{
    try
    {
        LoggersManagerObject.Log(STREAM("SAVING OF SAVING DATA TO BINARY FILE"));

        string ParticlesDataFileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("binary") + OS_DIR_SEP + string("ParticlesDataFile.dat");
        ofstream ParticlesDataFile(ParticlesDataFileName, ios_base::out | ios_base::trunc | ios_base::binary);

        SaveChemicalReactionsToBinaryFile(ParticlesDataFile);
        SaveParticlesKindsToBinaryFile(ParticlesDataFile);
        SaveParticlesToBinaryFile(ParticlesDataFile);
        SaveGenesToBinaryFile(ParticlesDataFile);

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

        SaveParticlesKindsAndParticlesAndChemicalReactionsAndGenesToBinaryFile();

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

        FOR_EACH_PARTICLE_IN_SECTORS_XYZ
            if (CellEngineUseful::IsDNA(ParticleObject.second.EntityId) == false)
            {
                if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                CellEngineVoxelSimulationSpaceObjectPointer->SetAllVoxelsInListOfVoxelsToValueForOuterClass(ParticleObject.second.ListOfVoxels, ParticleObject.second.Index);

                ParticleObject.second.Prev = ParticleObject.second.PrevTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PrevTemporary) : nullptr;
                ParticleObject.second.Next = ParticleObject.second.NextTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.NextTemporary) : nullptr;
                ParticleObject.second.PairedNucleotidePtr = ParticleObject.second.PairedNucleotideTemporary != 0 ? &GetParticleFromIndex(ParticleObject.second.PairedNucleotideTemporary) : nullptr;

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
        VectorToBeRead.clear();

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
            ParticlesDataFile.read((char*)&ParticleKindObject.GeneId, sizeof(ParticleKindObject.GeneId));
            ParticlesDataFile.read((char*)&ParticleKindObject.ElectricCharge, sizeof(ParticleKindObject.ElectricCharge));
            ParticlesDataFile.read((char*)&ParticleKindObject.Counter, sizeof(ParticleKindObject.Counter));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeX, sizeof(ParticleKindObject.GraphicData.SizeX));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeY, sizeof(ParticleKindObject.GraphicData.SizeY));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.SizeZ, sizeof(ParticleKindObject.GraphicData.SizeZ));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.AtomColor, sizeof(ParticleKindObject.GraphicData.AtomColor));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.ParticleColor, sizeof(ParticleKindObject.GraphicData.ParticleColor));
            ParticlesDataFile.read((char*)&ParticleKindObject.GraphicData.RandomParticleColor, sizeof(ParticleKindObject.GraphicData.RandomParticleColor));
            ParticleKindObject.GraphicData.ParticleColor = ParticlesKindsManagerObject.GetParticleKindGraphicDataFromConfigXMLData(ParticleKindObject.EntityId);

            float XSizeDiv2, YSizeDiv2, ZSizeDiv2;
            ParticlesDataFile.read((char*)&XSizeDiv2, sizeof(XSizeDiv2));
            ParticlesDataFile.read((char*)&YSizeDiv2, sizeof(YSizeDiv2));
            ParticlesDataFile.read((char*)&ZSizeDiv2, sizeof(ZSizeDiv2));
            ParticleKindObject.XSizeDiv2 = XSizeDiv2;
            ParticleKindObject.YSizeDiv2 = YSizeDiv2;
            ParticleKindObject.ZSizeDiv2 = ZSizeDiv2;

            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.IdStr);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.Name);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.Formula);
            ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.Compartment);
	        ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromXML);
	        ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindObject.GraphicData.NameFromDataFile);

	        UnsignedInt ParticleKindSpecialDataSectorSize;
	        ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorSize, sizeof(ParticleKindSpecialDataSectorSize));
	        for (UnsignedInt ParticleKindSpecialDataSectorIndex = 0; ParticleKindSpecialDataSectorIndex < ParticleKindSpecialDataSectorSize; ParticleKindSpecialDataSectorIndex++)
            {
                ParticleKindSpecialData ParticleKindSpecialDataSectorObject;
                ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindSpecialDataSectorObject.Description);
                ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindSpecialDataSectorObject.AddedParticle);
                ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorObject.GeneId, sizeof(ParticleKindSpecialDataSectorObject.GeneId));
                ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorObject.IsProtein, sizeof(ParticleKindSpecialDataSectorObject.IsProtein));
                ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorObject.ParticleType, sizeof(ParticleKindSpecialDataSectorObject.ParticleType));
                ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorObject.CleanProductOfTranscription, sizeof(ParticleKindSpecialDataSectorObject.CleanProductOfTranscription));
                ParticlesDataFile.read((char*)&ParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation, sizeof(ParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation));
                ParticleKindObject.ParticleKindSpecialDataSector.emplace_back(ParticleKindSpecialDataSectorObject);
            }

            if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true || CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                ReadVectorFromBinaryFile<vector3_16>(ParticlesDataFile, ParticleKindObject.ListOfVoxels);
            else
            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                ReadVectorFromBinaryFile<CellEngineAtom>(ParticlesDataFile, ParticleKindObject.ListOfAtoms);

            ParticleKindObject.GraphicData.Visible = true;

            ParticlesKindsManagerObject.ParticlesKinds[ParticleKindObject.EntityId] = ParticleKindObject;
        }

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES KINDS FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}

void ReadVoxelsVectorDividedByStepsFromBinaryFile(ifstream& ParticlesDataFile, ListOfVoxelsType& VectorToBeRead, const RealType DivideFactor)
{
    try
    {
        const PositionInt Step = DivideFactor * DivideFactor * DivideFactor;

        VectorToBeRead.clear();

        UnsignedInt Size;
        ParticlesDataFile.read((char*)&Size, sizeof(Size));
        for (UnsignedInt Index = 0; Index < Size; Index++)
        {
            vector3_16 Object{};
            ParticlesDataFile.read((char*)&Object, sizeof(Object));
            if (Index % Step == 0)
                VectorToBeRead.emplace_back(vector3_16{ static_cast<PositionInt>(static_cast<RealType>(Object.X) / DivideFactor), static_cast<PositionInt>(static_cast<RealType>(Object.Y) / DivideFactor), static_cast<PositionInt>(static_cast<RealType>(Object.Z) / DivideFactor) });
        }
    }
    CATCH("reading vector from binary file")
}

void CheckCenterForSector(const Particle& ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector)
{
    const UnsignedInt X = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z).SectorPosX;
    const UnsignedInt Y = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z).SectorPosY;
    const UnsignedInt Z = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z).SectorPosZ;

    if (ParticlesInSector[X][Y][Z].Particles.find(ParticleObject.Index)->second.Center != ParticleObject.Center)
        cout << "Error for particle XYZ = " << ParticleObject.EntityId << " " << ParticleObject.Index << endl;
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesFromBinaryFile(ifstream& ParticlesDataFile)
{
    try
    {
        UnsignedInt BadParticlesCenters = 0;

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

            vector3_Real32 CenterReadObject{};
            ParticlesDataFile.read((char*)&CenterReadObject, sizeof(CenterReadObject));
            ParticleObject.Center = { static_cast<RealType>(CenterReadObject.X), static_cast<RealType>(CenterReadObject.Y), static_cast<RealType>(CenterReadObject.Z) };

            ParticleObject.Center.X /= CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles;
            ParticleObject.Center.Y /= CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles;
            ParticleObject.Center.Z /= CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles;

            if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false && CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            {
                if (ParticleObject.Center.X > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.Center.Y > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.Center.Z > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                {
                    LoggersManagerObject.Log(STREAM(ParticleObject.EntityId << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << endl));
                    BadParticlesCenters++;
                }
            }

            ParticlesDataFile.read((char*)&ParticleObject.UniqueParticleColor, sizeof(ParticleObject.UniqueParticleColor));
            ParticlesDataFile.read((char*)&ParticleObject.SelectedForReaction, sizeof(ParticleObject.SelectedForReaction));

            ReadStringFromBinaryFile(ParticlesDataFile, ParticleObject.SequenceStr);

            if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true || CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            {
                ReadVoxelsVectorDividedByStepsFromBinaryFile(ParticlesDataFile, ParticleObject.ListOfVoxels, CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles);
                if (ParticleObject.ListOfVoxels.empty() == true)
                    cout << "Error for particle type = " << ParticleObject.EntityId << " " << ParticleObject.Index << endl;
            }
            else
            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                ReadVectorFromBinaryFile<CellEngineAtom>(ParticlesDataFile, ParticleObject.ListOfAtoms);

            ParticlesDataFile.read((char*)&ParticleObject.PrevTemporary, sizeof(ParticleObject.PrevTemporary));
            ParticlesDataFile.read((char*)&ParticleObject.NextTemporary, sizeof(ParticleObject.NextTemporary));
            ParticlesDataFile.read((char*)&ParticleObject.PairedNucleotideTemporary, sizeof(ParticleObject.PairedNucleotideTemporary));

            ReadVectorFromBinaryFile<UniqueIdInt>(ParticlesDataFile, ParticleObject.LinkedParticlesPointersListTemporary);

            if (ParticleObject.Index != 0)
            {
                ParticleObject.ParticleColor = ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).ParticleColor;
                ParticleObject.RandomParticleKindColor = ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).RandomParticleColor;
            }

            if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false)
            {
                auto ParticleSectorPos = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);

                if (AdditionalInfoPrinting == true)
                    cout << ParticleSectorPos.SectorPosX << " " << ParticleSectorPos.SectorPosY << " " << ParticleSectorPos.SectorPosZ << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " " << endl;

                if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                    SetCurrentSectorPos(ParticleSectorPos);

                if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true)
                {
                    if (Particles[ParticleSectorPos.SectorPosX][ParticleSectorPos.SectorPosY][ParticleSectorPos.SectorPosZ].MPIProcessIndex == MPIProcessDataObject.CurrentMPIProcessIndex)
                        AddNewParticle(ParticleObject);
                }
                else
                    AddNewParticle(ParticleObject);

                if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false)
                    if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                        CheckCenterForSector(ParticleObject, GetParticles());
            }

            if (ParticleObject.Index == 0)
                LoggersManagerObject.Log(STREAM("Particle with Index 0 = " << ParticleObject.Index << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z));
        }

        LoggersManagerObject.Log(STREAM("END OF READING PARTICLES FROM BINARY FILE - bad particles = " << BadParticlesCenters));

        ParticlesSize = 0;
        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
            ParticlesSize++;
        LoggersManagerObject.Log(STREAM("Number of Particles read = " << ParticlesSize));
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadChemicalReactionsFromBinaryFile(ifstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING CHEMICAL REACTIONS FROM BINARY FILE"));

        UnsignedInt ChemicalReactionsSize;
        ParticlesDataFile.read((char*)&ChemicalReactionsSize, sizeof(ChemicalReactionsSize));
        LoggersManagerObject.Log(STREAM("Number of chemical reactions to be read = " << ChemicalReactionsSize));

        for (UnsignedInt ChemicalReactionObjectIndex = 1; ChemicalReactionObjectIndex <= ChemicalReactionsSize; ChemicalReactionObjectIndex++)
        {
            ChemicalReaction ChemicalReactionObject;

            ParticlesDataFile.read((char*)&ChemicalReactionObject.ReactionIdNum, sizeof(ChemicalReactionObject.ReactionIdNum));
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactionIdStr);
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactionName);
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.ReactantsStr);
            ParticlesDataFile.read((char*)&ChemicalReactionObject.Duration, sizeof(ChemicalReactionObject.Duration));
            ParticlesDataFile.read((char*)&ChemicalReactionObject.Reversible, sizeof(ChemicalReactionObject.Reversible));
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.UpperFluxBound);
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.LowerFluxBound);
            ParticlesDataFile.read((char*)&ChemicalReactionObject.AdditionalParameter1, sizeof(ChemicalReactionObject.AdditionalParameter1));
            ParticlesDataFile.read((char*)&ChemicalReactionObject.AdditionalParameter2, sizeof(ChemicalReactionObject.AdditionalParameter2));
            ReadStringFromBinaryFile(ParticlesDataFile, ChemicalReactionObject.Comment);

            UnsignedInt ChemicalReactionReactantsSize;
            ParticlesDataFile.read((char*)&ChemicalReactionReactantsSize, sizeof(ChemicalReactionReactantsSize));
            for (UnsignedInt ChemicalReactionReactantIndex = 0; ChemicalReactionReactantIndex < ChemicalReactionReactantsSize; ChemicalReactionReactantIndex++)
            {
                ParticleKindForChemicalReaction ParticleKindForChemicalReactionObject;

                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.EntityId, sizeof(ParticleKindForChemicalReactionObject.EntityId));
                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.Counter, sizeof(ParticleKindForChemicalReactionObject.Counter));
                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.ToRemoveInReaction, sizeof(ParticleKindForChemicalReactionObject.ToRemoveInReaction));
                ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindForChemicalReactionObject.SequenceStr);
                ReadVectorFromBinaryFile<ChainIdInt>(ParticlesDataFile, ParticleKindForChemicalReactionObject.Sequence);
                ReadVectorFromBinaryFile<UniqueIdInt>(ParticlesDataFile, ParticleKindForChemicalReactionObject.LinkedParticleTypes);

                ChemicalReactionObject.Reactants.emplace_back(ParticleKindForChemicalReactionObject);
            }

            UnsignedInt ChemicalReactionProductsSize;
            ParticlesDataFile.read((char*)&ChemicalReactionProductsSize, sizeof(ChemicalReactionProductsSize));
            for (UnsignedInt ChemicalReactionProductIndex = 0; ChemicalReactionProductIndex < ChemicalReactionProductsSize; ChemicalReactionProductIndex++)
            {
                ParticleKindForChemicalReaction ParticleKindForChemicalReactionObject;

                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.EntityId, sizeof(ParticleKindForChemicalReactionObject.EntityId));
                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.Counter, sizeof(ParticleKindForChemicalReactionObject.Counter));
                ParticlesDataFile.read((char*)&ParticleKindForChemicalReactionObject.ToRemoveInReaction, sizeof(ParticleKindForChemicalReactionObject.ToRemoveInReaction));
                ReadStringFromBinaryFile(ParticlesDataFile, ParticleKindForChemicalReactionObject.SequenceStr);
                ReadVectorFromBinaryFile<ChainIdInt>(ParticlesDataFile, ParticleKindForChemicalReactionObject.Sequence);
                ReadVectorFromBinaryFile<UniqueIdInt>(ParticlesDataFile, ParticleKindForChemicalReactionObject.LinkedParticleTypes);

                ChemicalReactionObject.Products.emplace_back(ParticleKindForChemicalReactionObject);
            }

            if (ChemicalReactionObject.ReactionName.substr(0,2) != "EX" || AddReactionsOfSuddenAppearanceAndDisappearance == true)
                ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReactionObject);
        }

        LoggersManagerObject.Log(STREAM("END OF READING CHEMICAL REACTIONS FROM BINARY FILE"));
    }
    CATCH("reading chemical reactions from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadGenesFromBinaryFile(ifstream& ParticlesDataFile)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING GENES FROM BINARY FILE"));

        UnsignedInt GenesSize;
        ParticlesDataFile.read((char*)&GenesSize, sizeof(GenesSize));
        LoggersManagerObject.Log(STREAM("Number of genes to be read = " << GenesSize));

        for (UnsignedInt GeneObjectIndex = 1; GeneObjectIndex <= GenesSize; GeneObjectIndex++)
        {
            Gene GeneObject;

            ParticlesDataFile.read((char*)&GeneObject.NumId, sizeof(GeneObject.NumId));
            ReadStringFromBinaryFile(ParticlesDataFile, GeneObject.StrId);
            ReadStringFromBinaryFile(ParticlesDataFile, GeneObject.Description);
            ReadStringFromBinaryFile(ParticlesDataFile, GeneObject.ProteinId);
            ReadStringFromBinaryFile(ParticlesDataFile, GeneObject.Sequence);
            ParticlesDataFile.read((char*)&GeneObject.StartPosInGenome, sizeof(GeneObject.StartPosInGenome));
            ParticlesDataFile.read((char*)&GeneObject.EndPosInGenome, sizeof(GeneObject.EndPosInGenome));

            ParticlesKindsManagerObject.Genes[GeneObject.NumId] = GeneObject;

            if (GeneObject.ProteinId.starts_with("30S"))
                ParticlesKindsManagerObject.Ribosomes30SProteinsList.emplace_back(GeneObject.NumId);
            else
            if (GeneObject.ProteinId.starts_with("50S"))
                ParticlesKindsManagerObject.Ribosomes50SProteinsList.emplace_back(GeneObject.NumId);
        }

        LoggersManagerObject.Log(STREAM("END OF READING GENES FROM BINARY FILE"));
    }
    CATCH("reading genes from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadParticlesKindsAndParticlesAndChemicalReactionsAndGenesFromBinaryFile(CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        LoggersManagerObject.Log(STREAM("START OF READING DATA FROM BINARY FILE"));

        Particles[0][0][0].Particles.clear();
        ParticlesKindsManagerObject.Genes.clear();
        ParticlesKindsManagerObject.ParticlesKinds.clear();
        ChemicalReactionsManagerObject.ClearAllChemicalReactions();

        string OpenFileName = (Type == CellEngineConfigData::TypesOfFileToRead::CIFFile ? CellEngineConfigDataObject.CellStateFileName : CellEngineConfigDataObject.CellStateFileNameBackup);
        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true)
            OpenFileName = CellEngineConfigDataObject.CellStateFileNameBackup;

        ifstream ParticlesDataFile(OpenFileName, ios_base::in | ios_base::binary);

        ReadChemicalReactionsFromBinaryFile(ParticlesDataFile);
        ReadParticlesKindsFromBinaryFile(ParticlesDataFile);
        ReadParticlesFromBinaryFile(ParticlesDataFile);
        ReadGenesFromBinaryFile(ParticlesDataFile);

        ParticlesDataFile.close();

        LoggersManagerObject.Log(STREAM("END OF READING DATA FROM BINARY FILE"));
    }
    CATCH("reading particles from binary file")
}

void CellEngineParticlesBinaryDataFileReaderWriter::PreprocessLinkAndAssociateEveryParticleKindWithProperChemicalReaction()
{
    try
    {
        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            ParticleKindObject.second.AssociatedChemicalReactions.clear();

            for (const auto &ChemicalReactionObject : ChemicalReactionsManagerObject.ChemicalReactions)
                for (const auto& ReactantObject : ChemicalReactionObject.Reactants)
                    if (ReactantObject.EntityId == ParticleKindObject.second.EntityId)
                    {
                        ParticleKindObject.second.AssociatedChemicalReactions.insert(ChemicalReactionObject.ReactionIdNum);
                        break;
                    }
        }

        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            string ListOfAssociatedReactions;
            for (const auto& AssociatedChemicalReactionObject : ParticleKindObject.second.AssociatedChemicalReactions)
                ListOfAssociatedReactions += to_string(AssociatedChemicalReactionObject) + " " + ChemicalReactionsManagerObject.GetReactionFromNumId(AssociatedChemicalReactionObject).ReactantsStr + " | ";

            LoggersManagerObject.Log(STREAM("PARTICLE KIND = " << to_string(ParticleKindObject.second.EntityId) << " " << ParticleKindObject.second.IdStr << " " << ParticleKindObject.second.Name));
            LoggersManagerObject.Log(STREAM("ListOfAssociatedReactions = " << ListOfAssociatedReactions));
        }
    }
    CATCH("preprocessing linking and associating every particle kind with proper chemical reaction")
}

void CellEngineParticlesBinaryDataFileReaderWriter::FindNucleotidesIdentifiers()
{
    try
    {
        CellEngineConfigDataObject.ATP_ID = ParticlesKindsManagerObject.GetParticleKindFromStrId("M_atp_c")->EntityId;
        CellEngineConfigDataObject.CTP_ID = ParticlesKindsManagerObject.GetParticleKindFromStrId("M_ctp_c")->EntityId;
        CellEngineConfigDataObject.GTP_ID = ParticlesKindsManagerObject.GetParticleKindFromStrId("M_gtp_c")->EntityId;
        CellEngineConfigDataObject.TTP_ID = ParticlesKindsManagerObject.GetParticleKindFromStrId("M_ttp_c")->EntityId;
        CellEngineConfigDataObject.UTP_ID = ParticlesKindsManagerObject.GetParticleKindFromStrId("M_utp_c")->EntityId;
    }
    CATCH("finding nucleotides identifiers")
}

void CellEngineParticlesBinaryDataFileReaderWriter::FindGenomeParameters() const
{
    try
    {
        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false && CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
        {
            CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeDataFromFile(CellEngineConfigDataObject.DNAPaired);
            CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeSequenceFromFile(CellEngineConfigDataObject.DNAPaired);
            CellEngineVoxelSimulationSpaceObjectPointer->TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms();
        }
        else
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
        {
            CellEngineFullAtomSimulationSpaceObjectPointer->ReadGenomeSequenceFromFile(CellEngineConfigDataObject.DNAPaired);
            CellEngineFullAtomSimulationSpaceObjectPointer->TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms();
        }
    }
    CATCH("finding genome parameters")
}

void CellEngineParticlesBinaryDataFileReaderWriter::ReadAllDataFromBinaryFileAndPrepareData(const bool StartValuesBool, const bool UpdateParticleKindListOfVoxelsBool, const CellEngineConfigData::TypesOfFileToRead Type)
{
    try
    {
        if (StartValuesBool == true)
            SetStartValues();

        LoggersManagerObject.Log(STREAM("START OF READING ALL DATA FROM BINARY FILE AND PREPARING DATA"));

        ReadParticlesKindsAndParticlesAndChemicalReactionsAndGenesFromBinaryFile(Type);

        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false)
        PrepareParticlesAfterReadingFromBinaryFile();

        PreprocessData(UpdateParticleKindListOfVoxelsBool);

        ChemicalReactionsManagerObject.PreprocessAllChemicalReactions();

        PreprocessLinkAndAssociateEveryParticleKindWithProperChemicalReaction();

        CellEngineVoxelSimulationSpaceObjectPointer->AddSpecialParticlesKinds();

        FindNucleotidesIdentifiers();

        CellEngineAminoAcidsManagerObject.MapAminoAcidsForProperCodons();

        CellEngineAminoAcidsManagerObject.MapAminoAcidsIdToAminoAcidsObject();

        FindGenomeParameters();

        CellEngineConfigDataObject.GenomeReadFromFile = true;

        LoggersManagerObject.Log(STREAM("END OF PREPROCESSING DATA"));

        LoggersManagerObject.Log(STREAM("END OF READING ALL DATA FROM BINARY FILE AND PREPARING DATA"));
    }
    CATCH("reading particles from file")
};
