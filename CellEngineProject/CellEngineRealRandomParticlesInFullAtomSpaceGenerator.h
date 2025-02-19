
#ifndef CELL_ENGINE_REAL_RANDOM_PARTICLES_IN_FULL_ATOM_SPACE_GENERATOR_H
#define CELL_ENGINE_REAL_RANDOM_PARTICLES_IN_FULL_ATOM_SPACE_GENERATOR_H

#include "CellEngineParticleKind.h"
#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesFullAtomShapesGenerator.h"

class CellEngineRealRandomParticlesInFullAtomSpaceGenerator : public CellEngineParticlesFullAtomShapesGenerator, virtual public CellEngineRandomDeviceEngine, virtual public CellEngineBasicParticlesOperations
{
private:
    const UnsignedInt MaxNumberOfTriesToInsertNewParticle = 1000;
    const UnsignedInt MaxNumberOfTriesToFindARandomPositionInCell = 1000;
private:
    UnsignedInt TotalNumberOfAllParticles = 0;
public:
    void GenerateAllRealRandomParticles();
public:
    void UpdateSequence(ParticlesTypes ParticleTypeParam) const;
    static void ShowParticlesKindsData(ParticlesTypes ParticleTypeParam);
protected:
    [[nodiscard]] UnsignedInt GetNumberOfRealParticlesOfKind(ParticlesTypes ParticleTypeParam) const;
    std::tuple<UnsignedInt, UnsignedInt> GetNumberOfParticlesKind(ParticlesTypes ParticleTypeParam, bool AddToTotalNumberOfAllParticles);
    std::tuple<float, float, float> GetRandomPositionInsideSphere(float Radius, float RadiusSize);
public:
    void InsertNewRandomParticlesForType(ParticlesTypes ParticleTypeParam, bool ModifyParticleKindIdForRNA, UnsignedInt Radius, UnsignedInt RadiusSize);
    bool TryToGenerateRandomParticlesForType(const std::pair<const EntityIdInt, ParticleKind>& ParticleKindObject, EntityIdInt EntityId, UnsignedInt Radius, UnsignedInt RadiusSize, UnsignedInt& NumberOfErrors, GeneIdInt GeneNumId, const std::string& GeneSequence, UnsignedInt GeneSequenceLength);
public:
    void RemoveParticlesWithChosenEntityId(EntityIdInt EntityId, UnsignedInt NumberOfParticlesToBeRemoved);
    void RemoveAllParticlesWithChosenEntityId(EntityIdInt EntityId);
    void RemoveParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam, UnsignedInt NumberOfParticlesToBeRemoved);
    void RemoveAllParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam);
public:
    void RemoveAllRNAParticles();
    void RemoveAllmRNAParticles();
    void RemoveAlltRNAParticles();
    void RemoveAllrRNAParticles();
public:
    void PrintNumberOfParticlesForAllMainTypesOfParticles();
    static void PrintInformationAboutRibosomesProteins();
protected:
    explicit CellEngineRealRandomParticlesInFullAtomSpaceGenerator(ParticlesContainer<Particle>& ParticlesParam) : CellEngineParticlesFullAtomShapesGenerator(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
