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
    virtual bool CheckInsertOfParticle(const MPIParticleSenderStruct& MPIParticleSenderToInsert) = 0;
    SignedInt GetProcessPrevNeighbour(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
    SignedInt GetProcessNextNeighbour(SignedInt ThreadXIndex, SignedInt ThreadYIndex, SignedInt ThreadZIndex) const;
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
    void JoinStatisticsFromThreads(std::vector<std::map<UnsignedInt, ReactionStatistics>>& SavedReactionsMap, UnsignedInt SimulationStepNumber) const;
private:
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint);
    void GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam) const;
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
    void GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, bool PrintTime);
private:
    void GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex);
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam);
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
    void ExchangeParticlesBetweenMPIProcessesVer1();
    void ExchangeParticlesBetweenMPIProcessesVer2();
    void ExchangeParticlesBetweenMPIProcessesGroup1();
    void ExchangeParticlesBetweenMPIProcessesGroup2Ver1();
    void ExchangeParticlesBetweenMPIProcessesGroup2Ver2();
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
