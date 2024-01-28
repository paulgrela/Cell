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

class CellEngineVoxelSimulationSpace : public CellEngineChemicalReactions
{
public:
    std::mt19937_64 mt64R{ std::random_device{}() };
private:
    std::vector<std::string> GenomesLines;
private:
    std::vector<std::vector<UniqueIdInt>> Genomes;
private:
    UnsignedInt MaxParticleIndex{};
    std::unordered_map<UniqueIdInt, Particle> Particles;
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
private:
    std::map<EntityIdInt, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdInt> ParticlesSortedByCapacityFoundInProximity;
private:
    std::vector<UniqueIdInt> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFoundInProximity;
private:
    void* SpacePointer;
private:
    inline SimulationSpaceVoxel& GetSpaceVoxel(UnsignedInt x, UnsignedInt y, UnsignedInt z);
    inline Particle& GetParticleFromIndex(UniqueIdInt ParticleIndex);
private:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    void SetStartValuesForSpaceMinMax();
public:
    void SetAllVoxelsInListOfVoxelsToValue(std::vector<vector3_16>& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue);
    void MakeZeroVoxelsForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    static void GetMinMaxOfCoordinates(UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    UnsignedInt SumOfNotEmptyVoxels{};
public:
    UniqueIdInt GetNewFreeIndexOfParticle();
    UniqueIdInt GetFreeIndexesOfParticleSize();
public:
    void CountStatisticsOfVoxelSimulationSpace();
    void SetAtomInVoxelSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
public:
    SimulationSpaceVoxel GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexInSimulationSpaceVoxel(UniqueIdInt ParticleIndex);
public:
    UniqueIdInt AddNewParticle(const Particle& ParticleParam);
public:
    void AddBasicParticlesKindsAndReactions();
public:
    void FillSquareParticle(UnsignedInt RPosX, UnsignedInt RPosY, UnsignedInt RPosZ, UniqueIdInt LocalNewParticleIndex, UniqueIdInt Size);
public:
    void GeneratePlanedParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
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
    bool FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor);
    bool FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    bool CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool ElongationStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool ElongationContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
public:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject);
    static bool CompareFitnessOfParticle(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticlesChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, std::vector<vector3_16>& Centers);
public:
    void GenerateOneStrand(EntityIdInt EntityId, std::string_view Sequence, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleStepX, UnsignedInt ParticleStepY, UnsignedInt ParticleStepZ);
    Particle* GenerateNucleotideParticle(Particle* ParticlePrev, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeThread, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, bool AddToGenome, std::vector<UniqueIdInt>& Genome, vector3_16 UniqueColorParam);
    std::tuple<Particle*, Particle*> GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt AddSizeX, UnsignedInt AddSizeY, UnsignedInt AddSizeZ, vector3_16 UniqueColorParam, bool Linked);
public:
    void GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5);
    void EraseAllDNAParticles();
    static void UpdateRandomPositions(UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, UnsignedInt Size);
    static bool TestFormerForbiddenPositions(std::unordered_set<std::string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt Size);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> EraseLastRandomDNAParticle(std::vector<UniqueIdInt>& Genome);
public:
    void GetMinMaxCoordinatesForDNA();
public:
    void SaveGenomeDataToFile(UnsignedInt ParticleSize);
    void ReadGenomeDataFromFile(bool Paired);
    void ReadGenomeSequenceFromFile();
    void TestGeneratedGenomeCorrectness(UnsignedInt ParticleSize);
public:
    void PreprocessData();
    void InitiateFreeParticleIndexes();
    void GetMinMaxCoordinatesForAllParticles();
    static void GetMinMaxCoordinatesForParticle(Particle& ParticleObject);
public:
    CellEngineVoxelSimulationSpace();
    ~CellEngineVoxelSimulationSpace();
};

#endif
