#ifndef CELL_ENGINE_SIMULATION_PARALLEL_EXECUTION_MANAGER_H
#define CELL_ENGINE_SIMULATION_PARALLEL_EXECUTION_MANAGER_H

#include <barrier>

#include "CellEngineTypes.h"
#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineSimulationSpace.h"

class ReactionStatistics;
class CellEngineSimulationSpace;

class CellEngineSimulationParallelExecutionManager : virtual public CellEngineBasicParticlesOperations
{
public:
    template <class SimulationSpaceType>
    static void CreateSimulationSpaceForParallelExecution(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);
public:
    virtual bool CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossible(const ParticleSenderStruct& ThreadsParticleSenderToInsert) = 0;
    [[nodiscard]] SignedInt GetProcessPrevNeighbor(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
    [[nodiscard]] SignedInt GetProcessNextNeighbor(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
    void CreateDataEveryMPIProcessForParallelExecution();
public:
    virtual void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;
    virtual void GenerateOneRandomReactionForSelectedSpace(RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, bool FindParticlesInProximityBool) = 0;
public:
    void ExchangeParticlesBetweenThreads(UnsignedInt StepOutside, bool StateOfSimulationSpaceDivisionForThreads, bool PrintInfo) const;
    void ExchangeParticlesBetweenThreadsParallelInsert(UnsignedInt StepOutside, bool StateOfSimulationSpaceDivisionForThreads, bool PrintInfo) const;
    void ExchangeParticlesBetweenThreadsParallelExtract(UnsignedInt StepOutside, bool StateOfSimulationSpaceDivisionForThreads, bool PrintInfo) const;
public:
    void CheckParticlesCenters(bool PrintAllParticles);
    void GatherParticlesFromThreadsToParticlesInMainThread();
    void FirstSendParticlesForThreads(bool PrintCenterOfParticleWithThreadIndex, bool PrintTime);
public:
    void GatherCancelledParticlesIndexesFromThreads();
    void GatherParticlesFromThreads();
public:
    void SaveFormerParticlesAsVectorElements();
public:
    void JoinReactionsStatisticsFromThreads(std::vector<std::map<UnsignedInt, ReactionStatistics>>& SavedReactionsMap, UnsignedInt SimulationStepNumber) const;
private:
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint);
    void GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam);
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
    void GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, bool PrintTime);
private:
    void GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex);
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam);
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
public:
    void ExchangeParticlesBetweenMPIProcessesVer2();
    void ExchangeParticlesBetweenMPIProcessesGroup1();
    void ExchangeParticlesBetweenMPIProcessesGroup2Ver2();
public:
    void ExchangeParticlesBetweenThreadsVer2ConditionalVariableTwoMutexes();
    void ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes();
    void ExchangeParticlesBetweenThreadsGroup2Ver2ConditionalVariableTwoMutexes();
public:
    void ExchangeParticlesBetweenThreadsVer2ConditionalVariableOneMutex();
    void ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex();
    void ExchangeParticlesBetweenThreadsGroup2Ver2ConditionalVariableOneMutex();
public:
    void SynchronizeWithNeighborByLocalBarrier() const;
    void ExchangeParticlesBetweenThreadsVer2LocalBarrier();
    void ExchangeParticlesBetweenThreadsVer2OneGlobalBarrier(std::barrier<>* SyncPoint);
public:
    void ExchangeParticlesBetweenThreadsGroup1Barrier();
    void ExchangeParticlesBetweenThreadsGroup2Ver2Barrier();
    void ExchangeParticlesBetweenThreadsGroup3Barrier();
private:
    void SetZeroForAllParallelExecutionVariables();
    void GatherAllParallelExecutionVariables();
private:
    ParticlesContainer<Particle> GatherParticlesToExchangeBetweenThreads(UnsignedInt TypeOfGet, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, UnsignedInt& ExchangedParticleCounter, bool StateOfSimulationSpaceDivisionForThreads, bool PrintInfo) const;
private:
    SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& SimulationSpaceDataForThreads;
public:
    CellEngineSimulationParallelExecutionManager();
};

#endif
