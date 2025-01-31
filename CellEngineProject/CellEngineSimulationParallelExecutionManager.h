#ifndef CELL_ENGINE_SIMULATION_PARALLEL_EXECUTION_MANAGER_H
#define CELL_ENGINE_SIMULATION_PARALLEL_EXECUTION_MANAGER_H

#include <barrier>

#include "CellEngineTypes.h"
#include "CellEngineBasicParticlesOperations.h"

class ReactionStatistics;
class CellEngineSimulationSpace;

class CellEngineSimulationParallelExecutionManager : virtual public CellEngineBasicParticlesOperations
{
public:
    template <class SimulationSpaceType>
    static void CreateSimulationSpaceForParallelExecution(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);
public:
    virtual void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;
    //virtual void GenerateOneRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, bool FindParticlesInProximityBool) = 0;
    virtual void GenerateOneRandomReactionForSelectedSpace(float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam, bool FindParticlesInProximityBool) = 0;
public:
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads);
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
public:
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
    void GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, bool PrintTime);
private:
    void SetZeroForAllParallelExecutionVariables();
    void GatherAllParallelExecutionVariables();
private:
    void GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, ThreadIdType CurrentThreadIndexParam, UnsignedInt ThreadXIndexParam, UnsignedInt ThreadYIndexParam, UnsignedInt ThreadZIndexParam) const;
private:
    vector<vector<vector<ParticlesContainerInternal<Particle>>>> GatherParticlesToExchangeBetweenThreads(UnsignedInt TypeOfGet, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, UnsignedInt& ExchangedParticleCounter, bool StateOfSimulationSpaceDivisionForThreads, bool PrintInfo) const;
private:
    std::vector<std::vector<std::vector<std::shared_ptr<CellEngineSimulationSpace>>>>& SimulationSpaceDataForThreads;
public:
    CellEngineSimulationParallelExecutionManager();
};

#endif
