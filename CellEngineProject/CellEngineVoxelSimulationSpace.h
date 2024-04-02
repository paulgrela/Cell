#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include <stack>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
#include "CellEngineConfigData.h"
#include "CellEngineChemicalReactions.h"

using SimulationSpaceVoxel = UniqueIdInt;

constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024 = 1024;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048 = 2048;

using Space_1024_1024_1024 = SimulationSpaceVoxel[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024];
using Space_2048_2048_2048 = SimulationSpaceVoxel[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048];

enum class ComparisonType { ByVectorLoop, ByString };

#define SIMULATION_DETAILED_LOG

static void SwitchOffLogs()
{
#ifdef SIMULATION_DETAILED_LOG
    LoggersManagerObject.InitializePrintingParameters(false, false, false, false, false, false, false, false, false, false, false, false, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);
#endif
}

static void SwitchOnLogs()
{
#ifdef SIMULATION_DETAILED_LOG
    LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);
#endif
}

class BasicVoxelsOperations
{
protected:
    void* SpacePointer;
protected:
    inline SimulationSpaceVoxel& GetSpaceVoxel(UnsignedInt x, UnsignedInt y, UnsignedInt z)
    {
        return CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? (*static_cast<Space_2048_2048_2048*>(SpacePointer))[x][y][z] : (*static_cast<Space_1024_1024_1024*>(SpacePointer))[x][y][z];
    }
};

class CellEngineParticlesVoxelsOperations : public BasicVoxelsOperations
{
public:
    inline void SetAllVoxelsInListOfVoxelsToValue(std::vector<vector3_16>& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue);
    void SetAllVoxelsInListOfVoxelsToValueOut(std::vector<vector3_16>& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue);
public:
    inline void MoveParticleByVector(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ);
    inline void MoveParticleNearOtherParticle(Particle& ParticleObject, const Particle& NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ);
    inline bool CheckFreeSpaceForParticleMovedByVector(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ);
    inline static bool CheckBoundsForParticleMovedByVector(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    inline bool CheckBoundsAndFreeSpaceForParticleMovedByVector(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    inline bool MoveParticleByVectorIfVoxelSpaceIsEmpty(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ);
    inline bool MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(Particle& ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    inline bool MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(Particle& ParticleObject, const Particle& NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ);
    inline void MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(Particle& ParticleObject, const Particle& NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ);
};

class CellEngineBasicParticlesOperations : public CellEngineParticlesVoxelsOperations
{
protected:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    UnsignedInt MaxParticleIndex{};
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
public:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
public:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return Particles[ParticleIndex];
    }
public:
    void InitiateFreeParticleIndexes();
public:
    inline UniqueIdInt GetNewFreeIndexOfParticle();
    [[nodiscard]] UniqueIdInt GetFreeIndexesOfParticleSize() const;
public:
    UniqueIdInt AddNewParticle(const Particle& ParticleParam);
    virtual void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) = 0;
public:
    void PreprocessData(bool UpdateParticleKindListOfVoxelsBool);
public:
    void SetStartValuesForSpaceMinMax();
    void GetMinMaxCoordinatesForAllParticles(bool UpdateParticleKindListOfVoxelsBool);
    static void GetMinMaxOfCoordinates(UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam);
    static void GetMinMaxCoordinatesForParticle(Particle& ParticleObject, bool UpdateParticleKindListOfVoxels);
    static void UpdateParticleKindListOfVoxels(Particle& ParticleObject, UnsignedInt ParticleXMin, UnsignedInt ParticleXMax, UnsignedInt ParticleYMin, UnsignedInt ParticleYMax, UnsignedInt ParticleZMin, UnsignedInt ParticleZMax);
public:
    std::vector<UniqueIdInt> GetAllParticlesWithChosenEntityId(UniqueIdInt EntityId);
    UnsignedInt GetNumberOfParticlesWithChosenEntityId(UniqueIdInt EntityId);
public:
   explicit CellEngineBasicParticlesOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
    }
};

class CellEngineParticlesVoxelsShapesGenerator : virtual public CellEngineBasicParticlesOperations
{
public:
    virtual void ClearVoxelSpaceAndParticles() = 0;
public:
    typedef bool (CellEngineParticlesVoxelsShapesGenerator::*CheckFreeSpaceForSelectedSpaceType)(UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UniqueIdInt );
    typedef void (CellEngineParticlesVoxelsShapesGenerator::*SetValueToVoxelsForSelectedSpaceType)(std::vector<vector3_16>*, UniqueIdInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt);
public:
    inline void GenerateParticleVoxelsWhenSelectedSpaceIsFree(UnsignedInt LocalNewParticleIndex, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToVoxelsForSelectedSpaceType SetValueToVoxelsForSelectedSpace);
public:
    inline void SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(std::vector <vector3_16> *FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ);
public:
    inline bool CheckFreeSpaceInCuboidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UniqueIdInt ValueToCheck);
    inline bool CheckFreeSpaceForSphereSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
    inline bool CheckFreeSpaceForEllipsoidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
public:
    inline void SetValueToVoxelsForCuboidSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    inline void SetValueToVoxelsForSphereSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
    inline void SetValueToVoxelsForEllipsoidSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
public:
    explicit CellEngineParticlesVoxelsShapesGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineRealRandomParticlesGenerator : public CellEngineParticlesVoxelsShapesGenerator
{
public:
    void GenerateAllRealRandomParticles();
    void GenerateRealRandomMembraneParticles();
    void GenerateRealRandomRibosomesParticles();
public:
    explicit CellEngineRealRandomParticlesGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineParticlesVoxelsShapesGenerator(ParticlesParam)
    {
    }
};

class CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator : public CellEngineRealRandomParticlesGenerator
{
public:
    std::mt19937_64 mt64R{ std::random_device{}() };
protected:
    std::vector<std::string> GenomesLines;
    std::vector<std::vector<UniqueIdInt>> Genomes;
public:
    void GenerateOneStrand(EntityIdInt EntityId, std::string_view Sequence, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleStepX, UnsignedInt ParticleStepY, UnsignedInt ParticleStepZ);
    Particle* GenerateNucleotideParticle(Particle* ParticlePrev, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeThread, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, bool AddToGenome, std::vector<UniqueIdInt>& Genome, vector3_16 UniqueColorParam, bool LinkWithPreviousNucleotide);
    std::tuple<Particle*, Particle*> GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt AddSizeX, UnsignedInt AddSizeY, UnsignedInt AddSizeZ, vector3_16 UniqueColorParam, bool Linked, bool LinkWithPreviousNucleotide);
public:
    void EraseAllDNAParticles();
public:
    void GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5);
protected:
    static void UpdateRandomPositions(UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, UnsignedInt Size);
    static bool TestFormerForbiddenPositions(std::unordered_set<std::string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt Size);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> EraseLastRandomDNAParticle(std::vector<UniqueIdInt>& Genome);
protected:
    void GetMinMaxCoordinatesForDNA();
public:
    void SaveGenomeDataToFile(UnsignedInt ParticleSize);
    void ReadGenomeDataFromFile(bool Paired);
    void ReadGenomeSequenceFromFile();
    void TestGeneratedGenomeCorrectness(UnsignedInt ParticleSize);
public:
    explicit CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam), CellEngineRealRandomParticlesGenerator(ParticlesParam)
    {
    }
};



class CellEngineChemicalReactionsInVoxelSpace : virtual public CellEngineBasicParticlesOperations
{
protected:
    std::map<EntityIdInt, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdInt> ParticlesSortedByCapacityFoundInProximity;
protected:
    std::vector<UniqueIdInt> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFoundInProximity;
public:
    static bool CompareFitnessOfParticle(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, std::vector <vector3_16> &Centers);
public:
    explicit CellEngineChemicalReactionsInVoxelSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineNucleicAcidsBasicOperations
{
public:
    inline std::string GetPairedSequenceStr(std::string SequenceStr);
    inline std::tuple<Particle*, Particle*, UnsignedInt, std::string, std::vector<ChainIdInt>> GetNucleotidesSequence(Particle* Particle::*Direction, UnsignedInt LengthOfSequence, Particle& ParticleObjectForReaction,bool ToString, bool ToVector, bool (*Predicate)(const Particle*));
    inline void CutDNAPrev(Particle* NucleotideObjectForReactionPtr);
    inline void CutDNANext(Particle* NucleotideObjectForReactionPtr);
    inline void LinkDNA(Particle* NucleotideObjectForReactionPtr1, Particle* NucleotideObjectForReactionPtr2);
    inline void SeparateTwoPairedDNANucleotides(Particle* ParticleNucleotide);
    inline void SeparateDNAStrands(Particle* Particle::*Direction, Particle* NucleotidePtr, UnsignedInt LengthOfStrand);
    inline void JoinDNAStrands(Particle* Particle::*Direction, Particle* Strand1, Particle* Strand2);
    inline void DeleteLinkedParticlesPointersList(Particle& ParticleObject);
};

class CellEngineNucleicAcidsChemicalReactionsInVoxelSpace : public CellEngineChemicalReactionsInVoxelSpace, public CellEngineNucleicAcidsBasicOperations
{
public:
    void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) override;
public:
    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
public:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject);
public:
    void MakingZeroSizeForContainersForFoundParticlesInProximity();
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex);
    void PrintInformationAboutFoundParticlesInProximity();
public:
    bool FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor);
    bool FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    bool CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
public:
    explicit CellEngineNucleicAcidsChemicalReactionsInVoxelSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInVoxelSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};




class CellEngineVoxelSimulationSpace : public CellEngineChemicalReactions, public CellEngineNucleicAcidsChemicalReactionsInVoxelSpace, public CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator
{
    friend class CellEngineParticlesDataFile;
public:
    std::mt19937_64 mt64R{ std::random_device{}() };
private:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    UnsignedInt SumOfNotEmptyVoxels{};
public:
    void CountStatisticsOfVoxelSimulationSpace();
    void SetAtomInVoxelSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
public:
    SimulationSpaceVoxel GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexInSimulationSpaceVoxel(UniqueIdInt ParticleIndex);
public:
    void ClearSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    static void SaveParticlesToFile();
    void ReadParticlesFromFileAndPrepareData();
public:
    static void AddParticlesKinds();
    void AddChemicalReactions();
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void GeneratePlanedEllipsoidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GeneratePlanedCuboidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor);
public:
    void GenerateRandomReactionsForAllParticles();
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
public:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReaction();
    void FindAndExecuteChosenReaction(UnsignedInt ReactionId);
public:
    void GenerateRandomReactionForParticle(Particle& ParticleObject);
    void GenerateRandomReactionForSelectedVoxelSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateChosenReactionForSelectedVoxelSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) override;
    bool IsChemicalReactionPossible(const Reaction& ReactionObject) override;
    bool MakeChemicalReaction(Reaction& ReactionObject) override;
public:
    explicit CellEngineVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam);
    ~CellEngineVoxelSimulationSpace();
};

#endif
