
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H
#define CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H

#include "CellEngineReaction.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineChemicalReactionsInVoxelSpace.h"

class CellEngineNucleicAcidsChemicalReactionsInVoxelSpace : public CellEngineChemicalReactionsInVoxelSpace, public CellEngineNucleicAcidsBasicOperations
{
protected:
    enum class ComparisonType { ByVectorLoop, ByString };
protected:
    void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) override;
protected:
    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
protected:
    void MakingZeroSizeForContainersForFoundParticlesInProximity();
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex);
    void PrintInformationAboutFoundParticlesInProximity();
protected:
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
protected:
    explicit CellEngineNucleicAcidsChemicalReactionsInVoxelSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInVoxelSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
