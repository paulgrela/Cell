
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

#include <cmath>
#include <memory>

#include "CellEngineMacros.h"
#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"

constexpr bool PrintInfoWarning = false;

struct hash_pair
{
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        const size_t hash1 = std::hash<T1>{}(p.first);
        const size_t hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};

class CellEngineParticlesFullAtomOperations
{
protected:
    static inline std::unordered_map<std::pair<UnsignedInt, UnsignedInt>, std::unique_ptr<std::mutex>, hash_pair> InsertExtractParticleMutexObject;
public:
    static void SetProperThreadIndexForEveryParticlesSector(ParticlesContainer<Particle>& Particles)
    {
        try
        {
            UnsignedInt ThreadIndexPos = 1;
            FOR_EACH_PARTICLE_IN_XYZ_ONLY
            {
                Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].CurrentThreadPos = { ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation + 1, ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation + 1, ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation + 1 };
                Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].CurrentThreadIndex = ThreadIndexPos;
                ThreadIndexPos++;
            }
        }
        CATCH("setting proper thread index for every particles sector")
    }
public:
    static void GetMemoryForMutexesForParallelExecutionInParticlesSectors(ParticlesContainer<Particle>& Particles)
    {
        try
        {
            UnsignedInt CreatedMutexesCounter = 0;

            for (SignedInt ParticleSectorXIndex1 = 0; ParticleSectorXIndex1 < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex1++)
                for (SignedInt ParticleSectorYIndex1 = 0; ParticleSectorYIndex1 < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex1++)
                    for (SignedInt ParticleSectorZIndex1 = 0; ParticleSectorZIndex1 < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex1++)
                        for (SignedInt ParticleSectorXIndex2 = 0; ParticleSectorXIndex2 < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex2++)
                            for (SignedInt ParticleSectorYIndex2 = 0; ParticleSectorYIndex2 < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex2++)
                                for (SignedInt ParticleSectorZIndex2 = 0; ParticleSectorZIndex2 < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex2++)
                                    if (Particles[ParticleSectorXIndex1][ParticleSectorYIndex1][ParticleSectorZIndex1].CurrentThreadPos != Particles[ParticleSectorXIndex2][ParticleSectorYIndex2][ParticleSectorZIndex2].CurrentThreadPos)
                                        if (std::abs(ParticleSectorXIndex2 - ParticleSectorXIndex1) <= 1 && std::abs(ParticleSectorYIndex2 - ParticleSectorYIndex1) <= 1 && std::abs(ParticleSectorZIndex2 - ParticleSectorZIndex1) <= 1)
                                        {
                                            InsertExtractParticleMutexObject[std::make_pair(Particles[ParticleSectorXIndex1][ParticleSectorYIndex1][ParticleSectorZIndex1].CurrentThreadIndex, Particles[ParticleSectorXIndex2][ParticleSectorYIndex2][ParticleSectorZIndex2].CurrentThreadIndex)] = std::make_unique<std::mutex>();

                                            if (PrintInfoWarning == true)
                                                LoggersManagerObject.Log(STREAM("S1 = " << Particles[ParticleSectorXIndex1][ParticleSectorYIndex1][ParticleSectorZIndex1].CurrentThreadIndex << "(" << ParticleSectorXIndex1  << "," << ParticleSectorYIndex1 << "," << ParticleSectorZIndex1 << ") S2 = " << Particles[ParticleSectorXIndex2][ParticleSectorYIndex2][ParticleSectorZIndex2].CurrentThreadIndex << "(" << ParticleSectorXIndex2  << "," << ParticleSectorYIndex2 << "," << ParticleSectorZIndex2 << ")"));

                                            CreatedMutexesCounter++;
                                        }

            LoggersManagerObject.Log(STREAM("Created Mutex Counter = " << CreatedMutexesCounter));
        }
        CATCH("getting memory for mutexes for parallel execution in particles sectors")
    }
protected:
    static inline void MoveAllAtomsInParticleAtomsListByVector(Particle &ParticleObject, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        try
        {
            for (auto& AtomObject : ParticleObject.ListOfAtoms)
            {
                AtomObject.X += VectorX;
                AtomObject.Y += VectorY;
                AtomObject.Z += VectorZ;
            }
        }
        CATCH_AND_THROW("moving all voxels in particle voxel list by vector")
    }
protected:
    static inline void MoveParticleByVector(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        try
        {
            auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
            auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
            {
                //std::lock_guard<std::mutex> LockGuard{ *InsertExtractParticleMutexObject[std::make_pair(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].CurrentThreadIndex, ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].CurrentThreadIndex)] };

                if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
                    ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
                else
                {
                    if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObject.Index) == false)
                        if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObject.Index)->second == &ParticleObject)
                        {
                            if constexpr(PrintInfoWarning == true)
                                std::cout << "Earlier correct sector for particle index: " << ParticleObject.EntityId << " " << ParticleObject.Index << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                            return;
                        }
                }
            }
        }
        CATCH_AND_THROW("moving particle by vector")
    }
protected:
    static bool CheckSectorPos(const UnsignedInt SectorPosX, UnsignedInt SectorPosY, UnsignedInt SectorPosZ)
    {
        return (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
    }
protected:
    static bool CheckDistanceOfParticlesInSector(const RealType Radius, const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles)
            if (ParticleInSectorObject.second.Index != Index)
                if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                    return false;

        return true;
    }
protected:
    static bool CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles)
            if (ParticleInSectorObject.second.Index != Index)
                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                    if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                        return false;

        return true;
    }
protected:
    static PosType GetNewPosMovedByVector(const RealType X, const RealType Y, const RealType Z, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        return { X + VectorX, Y + VectorY, Z + VectorZ };
    }
protected:
    static inline bool CheckBoundsForSector(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (!(SectorPosX >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos) && SectorPosX < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX) && SectorPosY >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos) && SectorPosY < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY) && SectorPosZ >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos) && SectorPosZ < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ)))
            return false;

        return true;
    }
    static inline bool CheckBoundsForSpace(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (!(TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX && TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY && TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ))
            return false;

        return true;
    }
protected:
    static inline bool CheckBoundsBySectorAndBySpace(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckBounds == true)
        {
            if (CompareBoundsBySectorsBounds == true)
            {
                if (CheckBoundsForSector(SectorPosX, SectorPosY, SectorPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
            if (CompareBoundsBySpaceBounds == true)
            {
                if (CheckBoundsForSpace(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
        }

        return true;
    }
protected:
    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CheckBoundsBySectorAndBySpace(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                return false;

            if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                return false;

            if (CheckDistanceOfParticlesInSector(Radius, Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                return false;
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }
protected:
    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);

                if (CheckBoundsBySectorAndBySpace(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }
protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
    }
protected:
    //CCC - w sprawdzeniu czy nowa czastka w reakcji
    static inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
    {
         return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, { 0, 0, 0 }, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, false, true);
    }
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const bool CheckOnlyParticlesCenters)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ 0, 0, 0, 0, 0, 0, 0, 0, 0}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
                MoveParticleByVector(ParticleObject, ParticlesInSector, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ)
    {
        try
        {
            bool FoundFreeSpace = false;

            const RealType VecX = NewPositionParticleObject.ListOfAtoms[0].X - ParticleObject.ListOfAtoms[0].X;
            const RealType VecY = NewPositionParticleObject.ListOfAtoms[0].Y - ParticleObject.ListOfAtoms[0].Y;
            const RealType VecZ = NewPositionParticleObject.ListOfAtoms[0].Z - ParticleObject.ListOfAtoms[0].Z;

            for (RealType PosX = VecX - AddX; PosX < VecX + AddX; PosX += 1.0)
                for (RealType PosY = VecY - AddY; PosY < VecY + AddY; PosY += 1.0)
                    for (RealType PosZ = VecZ - AddZ; PosZ < VecZ + AddZ; PosZ += 1.0)
                        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, CellEngineConfigDataObject.CheckOnlyParticlesCenters) == true)
                        {
                            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "FREE SPACE FOUND " << VecX << " " << VecY << " " << VecZ << " " << PosX << " " << PosY << " " << PosZ << terminal_colors_utils::white));
                            MoveParticleByVector(ParticleObject, ParticlesInSector, PosX, PosY, PosZ);
                            FoundFreeSpace = true;
                            goto Outside;
                        }
            Outside:

            if (FoundFreeSpace == false)
                LoggersManagerObject.Log(STREAM(terminal_colors_utils::red << "FREE SPACE NOT FOUND " << VecX << " " << VecY << " " << VecZ << " " << terminal_colors_utils::white));
        }
        CATCH("moving particle near other particles")
    }
};

#endif
