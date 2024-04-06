
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H
#define CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H

#include <set>

#include "CellEngineReaction.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"

//class CellEngineNucleicAcidsComplexOperations : public CellEngineChemicalReactionsInBasicSimulationSpace, public CellEngineNucleicAcidsBasicOperations
//{
//public:
//    bool CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//    bool PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
//protected:
//    explicit CellEngineNucleicAcidsComplexOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
//    {
//    }
//};
//
//class CellEngineChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsComplexOperations
//{
//protected:
//    enum class ComparisonType { ByVectorLoop, ByString };
//protected:
//    void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) override;
//protected:
//    void MakingZeroSizeForContainersForFoundParticlesInProximity();
//    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex);
//    void PrintInformationAboutFoundParticlesInProximity();
//protected:
//    bool FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor);
//
//    bool FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
//
//    virtual void aaa(std::set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;
//
//    void ccc(UniqueIdInt ParticleIndex, std::set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides);
//protected:
//    explicit CellEngineChemicalReactionsInSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsComplexOperations(ParticlesParam)
//    {
//    }
//};
//
//class CellEngineNucleicAcidsChemicalReactionsInSimulationSpace : public CellEngineChemicalReactionsInSimulationSpace
//{
//protected:
//    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
//    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
//protected:
//    explicit CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
//    {
//    }
//};
//
//class CellEngineNucleicAcidsChemicalReactionsInVoxelSimulationSpace : public CellEngineNucleicAcidsChemicalReactionsInSimulationSpace
//{
//protected:
//    void aaa(set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
//protected:
//    explicit CellEngineNucleicAcidsChemicalReactionsInVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
//    {
//    }
//};

class CellEngineNucleicAcidsComplexOperations : public CellEngineChemicalReactionsInBasicSimulationSpace, public CellEngineNucleicAcidsBasicOperations
{
public:
    bool CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
    bool PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject);
protected:
    explicit CellEngineNucleicAcidsComplexOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineNucleicAcidsChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsComplexOperations
{
protected:
    enum class ComparisonType { ByVectorLoop, ByString };
protected:
    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
protected:
    explicit CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsComplexOperations(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsChemicalReactionsInSimulationSpace
{
protected:
    void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) override;
protected:
    void MakingZeroSizeForContainersForFoundParticlesInProximity();
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex);
    void PrintInformationAboutFoundParticlesInProximity();
protected:
    bool FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor);

    bool FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);

    virtual void aaa(std::set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;

    void ccc(UniqueIdInt ParticleIndex, std::set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides);
protected:
    explicit CellEngineChemicalReactionsInSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesParam)
    {
    }
};

class CellEngineNucleicAcidsChemicalReactionsInVoxelSimulationSpace : public CellEngineChemicalReactionsInSimulationSpace
{
protected:
    void aaa(set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
protected:
    explicit CellEngineNucleicAcidsChemicalReactionsInVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
