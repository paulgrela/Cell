
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
    inline Particle& GetParticleFromIndex(const UniqueIdUnsignedInt ParticleIndex)
    {
        return Particles[0][0][0].Particles[ParticleIndex];
    }
public:
    ParticlesDetailedContainer<Particle>::iterator GetParticleIteratorFromIndex(const UniqueIdUnsignedInt ParticleIndex)
    {
        return Particles[0][0][0].Particles.find(ParticleIndex);
    }
    [[nodiscard]] auto GetParticleEnd() const
    {
        return Particles[0][0][0].Particles.end();
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
                        ParticlesInSectorsZPos.Particles.clear();
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
    [[nodiscard]] vmath::vec3 GetCenterForAllParticles()
    {
        vmath::vec3 Center(0.0, 0.0, 0.0);

        try
        {
            float NumberOfParticles = 0;

            FOR_EACH_SECTOR_IN_XYZ_ONLY
                for (auto& ParticleObject : Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles | views::values)
                {
                    ParticleObject.OrderInParticleIndex = static_cast<UnsignedInt>(NumberOfParticles);

                    vmath::vec3 CenterOfParticle(0.0, 0.0, 0.0);

                    for (const CellEngineAtom& AtomObject : ParticleObject.ListOfAtoms)
                        CenterOfParticle += AtomObject.Position();

                    Center += vmath::vec3{ CenterOfParticle.X() / static_cast<float>(ParticleObject.ListOfAtoms.size()), CenterOfParticle.Y() / static_cast<float>(ParticleObject.ListOfAtoms.size()), CenterOfParticle.Z() / static_cast<float>(ParticleObject.ListOfAtoms.size()) };

                    NumberOfParticles++;
                }

            Center /= NumberOfParticles;
        }
        CATCH_AND_THROW("counting mass center")

        return Center;
    }
public:
    [[nodiscard]] UnsignedInt GetNumberOfAllParticles() const
    {
        UnsignedInt NumberOfParticles = 0;

        try
        {
            FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
                NumberOfParticles++;
        }
        CATCH_AND_THROW("getting number of particles")

        return NumberOfParticles;
    }
public:
    [[nodiscard]] UnsignedInt GetNumberOfAllParticlesInSector(const SectorPosType& SectorPos) const
    {
        UnsignedInt NumberOfParticles = 0;

        try
        {
            for (const auto& ParticleObject : Particles[SectorPos.SectorPosX][SectorPos.SectorPosY][SectorPos.SectorPosZ].Particles)
                NumberOfParticles++;
        }
        CATCH_AND_THROW("getting number of particles")

        return NumberOfParticles;
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