
#ifndef CELL_ENGINE_DATA_FILE_H
#define CELL_ENGINE_DATA_FILE_H

#include <memory>

#include "CellEngineMacros.h"

#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"
#include "CellEngineFilmOfStructures.h"
#include "CellEngineVoxelSimulationSpace.h"
#include "CellEngineFullAtomSimulationSpace.h"

class CellEngineDataFile : public CellEngineFilmOfStructures
{
public:
    CellEngineDataFile() = default;
public:
    ~CellEngineDataFile() override = default;
public:
    SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace> CellEngineSimulationSpaceForThreadsObjectsPointer;
    std::unique_ptr<CellEngineFullAtomSimulationSpace> CellEngineFullAtomSimulationSpaceObjectPointer;
    std::unique_ptr<CellEngineVoxelSimulationSpace> CellEngineVoxelSimulationSpaceObjectPointer;
public:
    bool FilmOfStructuresActive = false;
protected:
    ParticlesContainer<Particle> Particles;
protected:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return Particles[0][0][0][ParticleIndex];
    }
public:
    ParticlesContainerInternal<Particle>::iterator GetParticleIteratorFromIndex(const UniqueIdInt ParticleIndex)
    {
        return Particles[0][0][0].find(ParticleIndex);
    }
    [[nodiscard]] auto GetParticleEnd() const
    {
        return Particles[0][0][0].end();
    }
public:
    void GetMemoryForParticlesInSectors()
    {
        try
        {
            Particles.clear();
            Particles.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
            for (auto& ParticlesInSectorsXPos : Particles)
            {
                ParticlesInSectorsXPos.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInY);
                for (auto& ParticlesInSectorsYPos : ParticlesInSectorsXPos)
                {
                    ParticlesInSectorsYPos.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
                    for (auto& ParticlesInSectorsZPos : ParticlesInSectorsYPos)
                        ParticlesInSectorsZPos.clear();
                }
            }
        }
        CATCH("getting memory for particles in sectors")
    }
public:
    virtual void ReadDataFromFile(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type) = 0;
public:
    virtual void SaveDataToFile()
    {
    }
public:
    static vmath::vec3 GetCenter(const std::vector<CellEngineAtom>& AtomsParam)
    {
        vmath::vec3 Center(0.0, 0.0, 0.0);

        try
        {
            for (const CellEngineAtom& AtomObject : AtomsParam)
                Center += AtomObject.Position();
            Center /= static_cast<float>(AtomsParam.size());
        }
        CATCH_AND_THROW("counting mass center")

        return Center;
    }
public:
    vmath::vec3 GetCenterForAllParticles()
    {
        vmath::vec3 Center(0.0, 0.0, 0.0);

        try
        {
            float NumberOfParticles = 0;

            FOR_EACH_PARTICLE_IN_XYZ
            {
                vmath::vec3 CenterOfParticle(0.0, 0.0, 0.0);
                for (const CellEngineAtom& AtomObject : ParticleObject.second.ListOfAtoms)
                {
                    Center += AtomObject.Position();
                    CenterOfParticle += AtomObject.Position();
                    NumberOfParticles++;
                }
                ParticleObject.second.Center = { CenterOfParticle.X() / static_cast<float>(ParticleObject.second.ListOfAtoms.size()), CenterOfParticle.Y() / static_cast<float>(ParticleObject.second.ListOfAtoms.size()), CenterOfParticle.Z() / static_cast<float>(ParticleObject.second.ListOfAtoms.size()) };
            }
            Center /= NumberOfParticles;
        }
        CATCH_AND_THROW("counting mass center")

        //return Center;
        return { 0.0, 0.0, 0.0 };
    }
public:
    [[nodiscard]] ParticlesContainer<Particle>& GetParticles()
    {
        return Particles;
    }
    [[nodiscard]] UnsignedInt GetNumberOfStructures() override
    {
        return 0;
    }
};

inline std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;

#endif