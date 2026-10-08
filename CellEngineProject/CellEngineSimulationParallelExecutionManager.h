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
protected:
    virtual bool CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMPIMultiProcessDiffusion(const ParticleSenderStructMPIMultiProcess& ThreadsParticleSenderToInsert) = 0;
    virtual bool CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(const ParticleSenderStructMultiThreaded& ThreadsParticleSenderToInsert, const Particle& ParticleObject) = 0;
public:
    [[nodiscard]] SignedInt GetProcessPrevNeighbor(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
    [[nodiscard]] SignedInt GetProcessNextNeighbor(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
    void CreateDataEveryMPIProcessForParallelExecution();
    void CreateDataEveryThreadForParallelExecution() const;
    void CreateDataEveryThreadForParallelExecutionADD() const;
public:
    SignedInt GetProcessNeighbor(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
public:
    virtual void GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInMPIProcesses(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, bool InBounds, UnsignedInt StartSectorXPosParam, UnsignedInt StartSectorYPosParam, UnsignedInt StartSectorZPosParam, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;
    virtual void GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInThreads(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, bool InBounds, UnsignedInt StartSectorXPosParam, UnsignedInt StartSectorYPosParam, UnsignedInt StartSectorZPosParam, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;
    virtual void GenerateOneStepOfDiffusionForSelectedSpaceForExecutionNonParallel(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, bool InBounds, UnsignedInt StartSectorXPosParam, UnsignedInt StartSectorYPosParam, UnsignedInt StartSectorZPosParam, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;
public:
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
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint);
    void GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam, const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
    void GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, bool PrintTime);
private:
    void GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex);
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam);
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
private:
    static void SychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    void ExchangeParticlesBetweenThreadsAndSychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam);
public:
    void ExchangeParticlesBetweenMPIProcessesVer2();
    void ExchangeParticlesBetweenMPIProcessesGroup1();
    void ExchangeParticlesBetweenMPIProcessesGroup2Ver2();
public:
    void ExchangeParticlesBetweenThreadsVer2ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    void ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    static void ExchangeParticlesBetweenThreadsGroup2Ver2ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
public:
    static void ExchangeParticlesBetweenThreadsVer2ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    static void ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    static void ExchangeParticlesBetweenThreadsGroup2Ver2ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
public:
    static void SynchronizeWithNeighborByLocalBarrier(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
    void ExchangeParticlesBetweenThreadsVer2LocalBarrier(const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam) const;
    void ExchangeParticlesBetweenThreadsVer2OneGlobalBarrier(std::barrier<>* SyncPoint, const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam) const;
public:
    void ExchangeParticlesBetweenThreadsGroup2Ver2Barrier(const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam) const;
    static void ExchangeParticlesBetweenThreadsGroup3Barrier(const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData);
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
