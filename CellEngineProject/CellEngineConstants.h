
#ifndef CELL_ENGINE_CONSTANTS_H
#define CELL_ENGINE_CONSTANTS_H

#include <string>
#include "CellEngineTypes.h"

#include "../Common/Compilation/ConditionalCompilationConstants.h"

constexpr UnsignedInt MaxLengthOfGene = 16384;

constexpr UniqueIdUInt32Type UnknownParticleKindId = 99999;

constexpr UniqueIdUInt32Type StartParticleKindId = 100000;
constexpr UniqueIdUInt32Type StartChemicalReactionId = 10000;

constexpr RealType ParticleMoveSizeInDiffusionRealType = 10.0;

constexpr bool AdditionalSortParticlesInProximityByCapacity = false;

constexpr UnsignedInt GenomeLength = 543380;

constexpr ChainIdUInt16Type ACode = 1;
constexpr ChainIdUInt16Type CCode = 2;
constexpr ChainIdUInt16Type GCode = 3;
constexpr ChainIdUInt16Type TCode = 4;
constexpr ChainIdUInt16Type UCode = 5;

constexpr ChainIdUInt16Type TCodeAdd = 11;
constexpr ChainIdUInt16Type GCodeAdd = 12;
constexpr ChainIdUInt16Type CCodeAdd = 13;
constexpr ChainIdUInt16Type ACodeAdd = 14;
constexpr ChainIdUInt16Type UCodeAdd = 15;

constexpr ChainIdUInt16Type RCode = 21;
constexpr ChainIdUInt16Type YCode = 22;
constexpr ChainIdUInt16Type SCode = 23;
constexpr ChainIdUInt16Type WCode = 24;
constexpr ChainIdUInt16Type KCode = 25;
constexpr ChainIdUInt16Type MCode = 26;
constexpr ChainIdUInt16Type BCode = 27;
constexpr ChainIdUInt16Type DCode = 28;
constexpr ChainIdUInt16Type HCode = 29;
constexpr ChainIdUInt16Type VCode = 30;

constexpr ChainIdUInt16Type NCode = 31;

constexpr std::string RNAStartSequence = "ACGUACGU";

constexpr float AtomRadius = 1.0;

constexpr std::string JCVISYN3APredStr = "JCVISYN3A_";

constexpr UnsignedInt NumberOfAllNeighbors = 6;

constexpr UnsignedInt MaxMPIMessageSize = 1024 * 1024;

constexpr RealType UnusedValueRealType = 0.0;

#endif
