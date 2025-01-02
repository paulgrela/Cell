
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
    bool CutDNAInChosenPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNASingleStrandInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNABothStrandsInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNAInAnyPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNASingleStrandInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNABothStrandsInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool CutDNACrisperInChosenPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNASingleStrandCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool CutDNABothStrandsCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool PolymeraseRNATranscriptionStart(bool FullPromoter, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool PolymeraseRNATranscriptionStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool PolymeraseRNATranscriptionFullStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    void CheckEndingByHairpin(Particle& ParticleObject);
    static void CheckEndingByCodonStop(Particle& ParticleObject, const string& SequenceOfLettersToCheckFinishSequence);
    bool PolymeraseRNATranscriptionContinue(bool EndingByHairpin, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool PolymeraseRNATranscriptionContinueEndedByHairpinSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool PolymeraseRNATranscriptionContinueEndedByCodonStopSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
public:
    bool RibosomeTranslationStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
    bool RibosomeTranslationContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject);
protected:
    explicit CellEngineNucleicAcidsComplexOperations(ParticlesContainer<Particle>& ParticlesParam) : CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    ~CellEngineNucleicAcidsComplexOperations() override = default;
};


#endif
