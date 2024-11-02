
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
    bool CutDNAInChosenPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool CutDNASingleStrandInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool CutDNABothStrandsInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool LinkDNAInAnyPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const  CurrentThreadPosType& CurrentThreadPos);
    bool LinkDNASingleStrandInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool LinkDNABothStrandsInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool CutDNACrisperInChosenPlace(bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool CutDNASingleStrandCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool CutDNABothStrandsCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool PolymeraseRNATranscriptionStart(bool FullPromoter, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool PolymeraseRNATranscriptionStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool PolymeraseRNATranscriptionFullStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    void CheckEndingByHairpin(Particle& ParticleObject);
    static void CheckEndingByCodonStop(Particle& ParticleObject, const string& SequenceOfLettersToCheckFinishSequence);
    bool PolymeraseRNATranscriptionContinue(bool EndingByHairpin, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool PolymeraseRNATranscriptionContinueEndedByHairpinSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool PolymeraseRNATranscriptionContinueEndedByCodonStopSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
public:
    bool RibosomeTranslationStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
    bool RibosomeTranslationContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
protected:
    explicit CellEngineNucleicAcidsComplexOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    ~CellEngineNucleicAcidsComplexOperations() override = default;
};


#endif
