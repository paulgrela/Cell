
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_COMPLEX_OPERATIONS_H
#define CELL_ENGINE_NUCLEIC_ACIDS_COMPLEX_OPERATIONS_H

#include <set>

#include "CellEngineChemicalReaction.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"

bool CheckIfThisIsPromoter(UnsignedInt Box10Position);

class CellEngineNucleicAcidsComplexOperations : public CellEngineChemicalReactionsInBasicSimulationSpace, public CellEngineNucleicAcidsBasicOperations
{
protected:
    virtual void MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ) = 0;
public:
    bool CutDNAInChosenPlaceSpecialReactionFunction(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNASingleStrandInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNABothStrandsInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNAInAnyPlaceSpecialReactionFunction(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNASingleStrandInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNABothStrandsInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool CutDNACrisperInChosenPlaceSpecialReactionFunction(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNASingleStrandCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNABothStrandsCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool PolymeraseRNATranscriptionStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool PolymeraseRNATranscriptionContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool RibosomeTranslationStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool RibosomeTranslationContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
protected:
    explicit CellEngineNucleicAcidsComplexOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    ~CellEngineNucleicAcidsComplexOperations() override = default;
};


#endif
