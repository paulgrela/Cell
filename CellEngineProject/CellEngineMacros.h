#ifndef CELL_ENGINE_MACROS_H
#define CELL_ENGINE_MACROS_H

#define FOR_EACH_SECTOR_IN_XYZ_ONLY \
        for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++) \
            for (UnsignedInt ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++) \
                for (UnsignedInt ParticleSectorZIndex = 0; ParticleSectorZIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex++)

#define FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST \
        for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++) \
            for (UnsignedInt ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++) \
                for (UnsignedInt ParticleSectorZIndex = 0; ParticleSectorZIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex++) \
                    for (const auto& ParticleObject : Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)

#define FOR_EACH_PARTICLE_IN_SECTORS_XYZ \
        for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++) \
            for (UnsignedInt ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++) \
                for (UnsignedInt ParticleSectorZIndex = 0; ParticleSectorZIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex++) \
                    for (auto& ParticleObject : Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)

#endif
