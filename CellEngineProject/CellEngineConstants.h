
#ifndef CELL_ENGINE_CONSTANTS_H
#define CELL_ENGINE_CONSTANTS_H

#include <string>
#include "CellEngineTypes.h"

#define EXTENDED_RAM_MEMORY_
#define RNA_IN_ONE_PARTICLE

#define USING_MODULES_

constexpr UnsignedInt MaxLengthOfGene = 16384;

constexpr EntityIdInt StartParticleKindId = 100000;
constexpr EntityIdInt StartReactionId = 10000;

constexpr UniqueIdInt ParticleIndexesCreatorFactor = 10'000'000;

constexpr bool AdditionalSortParticlesInProximityByCapacity = false;


constexpr UnsignedInt GenomeLength = 543380;

constexpr ChainIdInt ACode = 1;
constexpr ChainIdInt CCode = 2;
constexpr ChainIdInt GCode = 3;
constexpr ChainIdInt TCode = 4;
constexpr ChainIdInt UCode = 5;

constexpr ChainIdInt TCodeAdd = 11;
constexpr ChainIdInt GCodeAdd = 12;
constexpr ChainIdInt CCodeAdd = 13;
constexpr ChainIdInt ACodeAdd = 14;
constexpr ChainIdInt UCodeAdd = 15;

constexpr ChainIdInt RCode = 21;
constexpr ChainIdInt YCode = 22;
constexpr ChainIdInt SCode = 23;
constexpr ChainIdInt WCode = 24;
constexpr ChainIdInt KCode = 25;
constexpr ChainIdInt MCode = 26;
constexpr ChainIdInt BCode = 27;
constexpr ChainIdInt DCode = 28;
constexpr ChainIdInt HCode = 29;
constexpr ChainIdInt VCode = 30;

constexpr ChainIdInt NCode = 31;

constexpr std::string RNAStartSequence = "ACGUACGU";

#endif
