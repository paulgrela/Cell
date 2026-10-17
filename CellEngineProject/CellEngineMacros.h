#ifndef CELL_ENGINE_MACROS_H
#define CELL_ENGINE_MACROS_H

#include "../Compilation/ConditionalCompilationConstants.h"

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

#ifdef SIMULATION_DETAILED_DEBUG_LOG
    #define DEBUG_LOG(...) __VA_ARGS__
#else
    #define DEBUG_LOG(...)
#endif

#ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
    #define DEBUG_LOG_DIFFUSION(...) __VA_ARGS__
#else
    #define DEBUG_LOG_DIFFUSION(...)
#endif

#ifdef SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION
    #define DEBUG_LOG_CHEMICAL_REACTION(...) __VA_ARGS__
#else
    #define DEBUG_LOG_CHEMICAL_REACTION(...)
#endif

#if defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION) || defined(SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION)
    #define DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(...) __VA_ARGS__
#else
    #define DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(...)
#endif

#endif
