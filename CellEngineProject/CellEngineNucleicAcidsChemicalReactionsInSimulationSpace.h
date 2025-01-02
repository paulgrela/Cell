
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_SIMULATION_SPACE_H
#define CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_SIMULATION_SPACE_H

#include <set>

#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineNucleicAcidsComplexOperations.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"

class CellEngineNucleicAcidsChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsComplexOperations
{
protected:
    enum class ComparisonType { ByVectorLoop, ByString };
protected:
    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
    std::tuple<std::vector<ChainIdInt>, std::string> GetNucleotidesSequenceFromRNAInOneParticle(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, UnsignedInt SizeOfLoop);
public:
    bool CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectTestedForReaction);
protected:
    explicit CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineNucleicAcidsComplexOperations(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
