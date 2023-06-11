#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include <stack>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
#include "CellEngineConfigData.h"

using SimulationSpaceVoxel = UniqueIdInt;

constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024 = 1024;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048 = 2048;

using Space_1024_1024_1024 = SimulationSpaceVoxel[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst1024];
using Space_2048_2048_2048 = SimulationSpaceVoxel[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048];

class CellEngineVoxelSimulationSpace
{
public:
    std::mt19937_64 mt64R{ std::random_device{}() };
private:
    std::string GenomeLine;
private:
    std::vector<UniqueIdInt> Genome1;
    std::vector<UniqueIdInt> Genome2;
private:
    UnsignedInt MaxParticleIndex{};
    std::unordered_map<UniqueIdInt, Particle> Particles;
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
private:
    std::vector<Reaction> Reactions;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
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
    void AddReaction(const Reaction& ReactionParam);
public:
    void AddBasicParticlesKindsAndReactions();
public:
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, std::map<EntityIdInt, UnsignedInt>& ParticlesKinds);
    void GenerateRandomReactionForParticle(Particle& ParticleObject);
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateRandomReactionsForAllParticles();
    void GenerateRandomReactionsInWholeVoxelSimulationSpace();
public:
    void FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor, std::map<EntityIdInt, UnsignedInt>& ParticlesKindsFoundInParticlesProximity, std::map<UnsignedInt, UniqueIdInt>& ParticlesSortedByCapacityFoundInParticlesProximity);
    std::vector<UniqueIdInt> ChooseParticlesForReactionFromAllParticlesInProximity(Reaction& ReactionObject, const std::map<UnsignedInt, UniqueIdInt>& ParticlesSortedByCapacityFoundInParticlesProximity);
    bool CompareFitnessOfDNASequenceByString(EntityIdInt ReactantEntityId, const Particle& ParticleObjectForReaction);
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(EntityIdInt ReactantEntityId, const Particle& ParticleObjectForReaction);
    void EraseParticlesChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, std::vector<vector3_16>& Centers);
public:
    void GenerateParticle(EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, std::vector<UniqueIdInt>& Genome, vector3_16 UniqueColorParam);
    void GenerateTwoPairedNucleotides(EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt AddSizeX, UnsignedInt AddSizeY, UnsignedInt AddSizeZ, vector3_16 UniqueColorParam);
public:
    void GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5);
    void EraseAllDNAParticles();
    static void UpdateRandomPositions(UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, UnsignedInt Size);
    static bool TestFormerForbiddenPositions(std::unordered_set<std::string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt Size);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> EraseLastRandomParticle(std::vector<UniqueIdInt>& Genome);
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
