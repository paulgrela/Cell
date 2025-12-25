
#include <omp.h>
#include <execution>

#include "../Common/Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartPositions()
{
    return { SelectionStartXPos, SelectionStartYPos, SelectionStartZPos };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetSizes()
{
    return { SelectionSizeX, SelectionSizeY, SelectionSizeZ };
}

bool CheckVisibility(const bool Visibility)
{
    if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile)
        return true;
    else
        return Visibility;
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix)
{
    try
    {
        // const auto start_time111 = chrono::high_resolution_clock::now();
        //
        // std::vector<GPUParticle> GPUParticles;
        // std::vector<GPUAtom> GPUAtoms;
        //
        // uint32_t atomOffset = 0;
        //
        // FOR_EACH_SECTOR_IN_XYZ_ONLY
        //     if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
        //     {
        //         bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex);
        //
        //         if (FinalVisibilityInModelWorld == true)
        //             if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
        //                 for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
        //                 {
        //                     GPUParticle GPUParticleObject;
        //                     GPUParticleObject.EntityId = ParticleObject.second.EntityId;
        //                     GPUParticleObject.ChainId = ParticleObject.second.ChainId;
        //                     GPUParticleObject.Index = ParticleObject.second.Index;
        //
        //                     GPUParticleObject.AtomOffset = atomOffset;
        //                     GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();
        //                     memset(GPUParticleObject._padding, 0, sizeof(GPUParticleObject._padding));
        //
        //                     GPUParticles.emplace_back(GPUParticleObject);
        //
        //                     for (const auto& atom : ParticleObject.second.ListOfAtoms)
        //                     {
        //                         GPUAtom gpuAtom;
        //                         gpuAtom.EntityId = atom.EntityId;
        //                         gpuAtom.X = atom.X;
        //                         gpuAtom.Y = atom.Y;
        //                         gpuAtom.Z = atom.Z;
        //                         // gpuAtom.ColorR = atom.AtomColor.X;
        //                         // gpuAtom.ColorG = atom.AtomColor.Y;
        //                         // gpuAtom.ColorB = atom.AtomColor.Z;
        //                         gpuAtom.AtomColor = atom.AtomColor;
        //
        //                         gpuAtom._padding1 = 0;
        //
        //                         // Copy strings with padding
        //                         memset(gpuAtom.Name, 0, 8);
        //                         strncpy(gpuAtom.Name, atom.Name, 6);
        //
        //                         memset(gpuAtom.ResName, 0, 8);
        //                         strncpy(gpuAtom.ResName, atom.ResName, 5);
        //
        //                         memset(gpuAtom.Chain, 0, 8);
        //                         strncpy(gpuAtom.Chain, atom.Chain, 7);
        //
        //                         GPUAtoms.push_back(gpuAtom);
        //                     }
        //
        //                     atomOffset += ParticleObject.second.ListOfAtoms.size();
        //                 }
        //     }
        //
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticleSSBO);
        // glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUParticles.size() * sizeof(GPUParticle), GPUParticles.data());
        //
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, AtomSSBO);
        // glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUAtoms.size() * sizeof(GPUAtom), GPUAtoms.data());
        //
        // // glDispatchCompute((MAX_PARTICLES + 255) / 256, 1, 1);
        // // glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        //
        // const auto stop_time111 = chrono::high_resolution_clock::now();
        //
        // ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);








        lock_guard LockGuard{ RenderMenuAndFullAtomSimulationSpaceMutexObject };

        GLuint PartOfStencilBufferIndex[3];
        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        UniformsBlocks.clear();
        //UniformsBlocks.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
        //UniformsBlocks.resize(128);
        //UniformsBlocks.resize(1);

        //std::vector<UniformsBlock> UniformsBlocksTotal{};

        const auto start_time = chrono::high_resolution_clock::now();

        //for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFile::ChosenStructureMutexObject };


            // const UnsignedInt chunksX = 4;
            // const UnsignedInt chunksY = 4;
            // const UnsignedInt chunksZ = 4;
            //
            // const UnsignedInt chunkSizeX = (CellEngineConfigDataObject.NumberOfParticlesSectorsInX + chunksX - 1) / chunksX;
            // const UnsignedInt chunkSizeY = (CellEngineConfigDataObject.NumberOfParticlesSectorsInY + chunksY - 1) / chunksY;
            // const UnsignedInt chunkSizeZ = (CellEngineConfigDataObject.NumberOfParticlesSectorsInZ + chunksZ - 1) / chunksZ;



            //omp_set_num_threads(64);
            //#pragma omp parallel for collapse(3) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool)  reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            //#pragma omp parallel for collapse(3) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool)  reduction(+:NumberOfAllRenderedAtoms) schedule(static)


            //FOR_EACH_SECTOR_IN_XYZ_ONLY
            // #pragma omp parallel num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix,  ParticlesKindsManagerObject, RenderObjectsBool, chunkSizeX, chunkSizeY, chunkSizeZ)  reduction(+:NumberOfAllRenderedAtoms)
            // {
            //
            //     int thread_id = omp_get_thread_num();
            //
            //     int chunk_x = thread_id / 16;        // 0-3
            //     int chunk_y = (thread_id / 4) % 4;   // 0-3
            //     int chunk_z = thread_id % 4;         // 0-3
            //
            //     UnsignedInt startX = chunk_x * chunkSizeX;
            //     UnsignedInt endX = std::min(startX + chunkSizeX, CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
            //
            //     UnsignedInt startY = chunk_y * chunkSizeY;
            //     UnsignedInt endY = std::min(startY + chunkSizeY, CellEngineConfigDataObject.NumberOfParticlesSectorsInY);
            //
            //     UnsignedInt startZ = chunk_z * chunkSizeZ;
            //     UnsignedInt endZ = std::min(startZ + chunkSizeZ, CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
            //
            //     for (UnsignedInt ParticleSectorXIndex = startX; ParticleSectorXIndex < endX; ParticleSectorXIndex++)
            //         for (UnsignedInt ParticleSectorYIndex = startY; ParticleSectorYIndex < endY; ParticleSectorYIndex++)
            //             for (UnsignedInt ParticleSectorZIndex = startZ; ParticleSectorZIndex < endZ; ParticleSectorZIndex++)

            //FOR_EACH_SECTOR_IN_XYZ_ONLY
            //omp_set_num_threads(40);
            // UnsignedInt ParticleSectorXIndex;
            // UnsignedInt ParticleSectorYIndex;
            // UnsignedInt ParticleSectorZIndex;

            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) private(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) num_threads(CellEngineConfigDataObject.NumberOfParticlesSectorsInX) reduction(+:NumberOfAllRenderedAtoms) schedule(static, 1)
            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) private(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) num_threads(CellEngineConfigDataObject.NumberOfParticlesSectorsInX) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            //omp_set_num_threads(40);
            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) private(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)

            //#pragma omp parallel for collapse(3) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            //#pragma omp parallel for collapse(3) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) reduction(+:NumberOfAllRenderedAtoms) schedule(static)

            //#pragma omp parallel for collapse(3) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            // for (ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++)
            //     for (ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++)
            //         for (ParticleSectorZIndex = 0; ParticleSectorZIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex++)


            //Visibility powoduje ze tylko fragment czastek jest zrobiony i watki nic nie daja
            //trzeba wiedziec gdzie jest widzialnosc i tam dzielic wiec -
            //pragma omp musi byc z sekscja krytytczna









            /*
            const auto start_time10 = chrono::high_resolution_clock::now();

            vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt>> ParticlesSectorsToBeRendered;

            // struct SectorXYZ
            // {
            //     UnsignedInt ParticleSectorXIndex;
            //     UnsignedInt ParticleSectorYIndex;
            //     UnsignedInt ParticleSectorZIndex;
            // };
            // vector<SectorXYZ> ParticlesSectorsToBeRendered;

            for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++)
                for (UnsignedInt ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++)
                    for (UnsignedInt ParticleSectorZIndex = 0; ParticleSectorZIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ; ParticleSectorZIndex++)
                        if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                        {
                            bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, 0);

                            FinalVisibilityInModelWorld = CheckVisibility(FinalVisibilityInModelWorld);

                            if (FinalVisibilityInModelWorld == true)
                                if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                                    ParticlesSectorsToBeRendered.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex);
                        }

            const auto stop_time10 = chrono::high_resolution_clock::now();
            ExecutionDurationTimeForCheckingPreparingParticles += chrono::duration(stop_time10 - start_time10);

            //UniformsBlocks.clear();
            //UniformsBlocks.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
            //UniformsBlocks.resize(128);
            //UniformsBlocks.resize(1);

            LoggersManagerObject.Log(STREAM("T=" << ParticlesSectorsToBeRendered.size()));


            const auto start_time2 = chrono::high_resolution_clock::now();

            vector<thread> Threads;
            UnsignedInt ParticlesSectorToBeRenderedIndex;
            //#pragma omp parallel for collapse(1) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered) num_threads(128) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)
            //#pragma omp parallel for collapse(1) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered) num_threads(128) reduction(+:NumberOfAllRenderedAtoms) schedule(static, 1)
            //#pragma omp parallel for collapse(1) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
            //#pragma omp parallel for collapse(1) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered) reduction(+:NumberOfAllRenderedAtoms) schedule(dynamic, 1)

            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
            for (ParticlesSectorToBeRenderedIndex = 0; ParticlesSectorToBeRenderedIndex < ParticlesSectorsToBeRendered.size(); ParticlesSectorToBeRenderedIndex++)
            {
                const UnsignedInt ParticleSectorXIndex = get<0>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
                const UnsignedInt ParticleSectorYIndex = get<1>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
                const UnsignedInt ParticleSectorZIndex = get<2>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);

                // for (const auto& ParticlesSectorsToBeRenderedObject : ParticlesSectorsToBeRendered)
                // //for (ParticlesSectorToBeRenderedIndex = 0; ParticlesSectorToBeRenderedIndex < ParticlesSectorsToBeRendered.size(); ParticlesSectorToBeRenderedIndex++)
                // {
                //     // const UnsignedInt ParticleSectorXIndex = get<0>(ParticlesSectorsToBeRenderedObject);
                //     // const UnsignedInt ParticleSectorYIndex = get<1>(ParticlesSectorsToBeRenderedObject);
                //     // const UnsignedInt ParticleSectorZIndex = get<2>(ParticlesSectorsToBeRenderedObject);
                //     const UnsignedInt ParticleSectorXIndex = ParticlesSectorsToBeRenderedObject.ParticleSectorXIndex;
                //     const UnsignedInt ParticleSectorYIndex = ParticlesSectorsToBeRenderedObject.ParticleSectorYIndex;
                //     const UnsignedInt ParticleSectorZIndex = ParticlesSectorsToBeRenderedObject.ParticleSectorZIndex;

                //LoggersManagerObject.Log(STREAM("R=" << CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size()));

                                        //UnsignedInt ParticleObjectIndex;
                                        //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) private(ParticleObjectIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                                        //for (ParticleObjectIndex = 0; ParticleObjectIndex < CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size(); ParticleObjectIndex++)
                //LoggersManagerObject.Log(STREAM("G=" << CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size()));

                // std::atomic<int> thread_counter{0};
                // thread_local int thread_index = -1;

                std::vector<ParticlesDetailedContainer<Particle>::iterator> iterators;
                iterators.reserve(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size());
                for (auto it = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin(); it != CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.end(); ++it)
                    iterators.push_back(it);

                size_t i;
                #pragma omp parallel for default(none) shared(iterators, CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered, ParticlesSectorToBeRenderedIndex, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) private(i) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                //#pragma omp parallel for
                for (i = 0; i < iterators.size(); i++)
                //for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                    // std::for_each(std::execution::par_unseq, CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin(), CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.end(),
                // [&, ViewMatrix, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex](auto& ParticleObject)
                {
                    //auto& ParticleObject = iterators[i];
                    //for (auto& i1 = iterators[i]->first; i1 != iterators[i]->second; i1++)

                    // if (thread_index == -1)
                    //     thread_index = thread_counter.fetch_add(1);
                    UnsignedInt ParticleObjectIndex = iterators[i]->first;
                    //auto& ParticleObject = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles[ParticleObjectIndex];
                    auto& ParticleObject = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles[ParticleObjectIndex];
                    //if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size() > 0)
                    {
                        //if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                        if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true)
                        {
                            //DrawBonds(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);
                            // Threads.emplace_back([&, ParticleSectorXIndex, this]()
                            // {

                            //ParticleObject.second.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };
                            ParticleObject.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };
                            //ParticleObject.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };

                            //for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                            for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                //for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                            {
                                // if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                // {
                                //     glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                //     TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                                // }

                                //RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, thread_index);
                                //RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, omp_get_thread_num());
                                RenderObject(ParticleObject.ListOfAtoms[AtomObjectIndex], ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, omp_get_thread_num());
                                //RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, 0);
                                //RenderObject(ParticleObject.ListOfAtoms[AtomObjectIndex], ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, 0);
                                //RenderObject(ParticleObject.ListOfAtoms[AtomObjectIndex], ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, 0);
                            }


                        }
                    }
                    //});
                }

                //CZY WATKI

                // UnsignedInt ParticleObjectIndex;
                // ParticlesDetailedContainer<Particle>& GetParticlesObject = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles;
                // LoggersManagerObject.Log(STREAM("G=" << GetParticlesObject.size()));
                // //if (GetParticlesObject.size() > 0)
                // {
                //     //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorsToBeRendered, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, GetParticlesObject) private(ParticleObjectIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                //     //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorToBeRenderedIndex, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex) private(GetParticlesObject, ParticleObjectIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                //     //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, GetParticlesObject) private(ParticleObjectIndex, ParticlesSectorToBeRenderedIndex) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                //     //#pragma omp parallel for collapse(1) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, ParticlesSectorToBeRenderedIndex, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, GetParticlesObject) num_threads(128) reduction(+:NumberOfAllRenderedAtoms)
                //
                //     //for (ParticleObjectIndex = 0; ParticleObjectIndex < CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size(); ParticleObjectIndex++)
                //     for (ParticleObjectIndex = 0; ParticleObjectIndex < GetParticlesObject.size(); ParticleObjectIndex++)
                //     {
                //         //auto& ParticleObject = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles[ParticleObjectIndex];
                //         auto& ParticleObject = GetParticlesObject[ParticleObjectIndex];
                //
                //          if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true)
                //          {
                //              //DrawBonds(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);
                //
                //              ParticleObject.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };
                //
                //              for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                //              {
                //                  // if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                //                  // {
                //                  //     glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                //                  //     TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                //                  // }
                //
                //                  //RenderObject(ParticleObject.ListOfAtoms[AtomObjectIndex], ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, omp_get_thread_num());
                //                  RenderObject(ParticleObject.ListOfAtoms[AtomObjectIndex], ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, 0);
                //              }
                //          }
                //     }
                // }
            }

            const auto stop_time2 = chrono::high_resolution_clock::now();
            ExecutionDurationTimeForCheckingPreparingParticles2 += chrono::duration(stop_time2 - start_time2);




















        const auto start_time111 = chrono::high_resolution_clock::now();

        std::vector<GPUParticle> GPUParticles;
        std::vector<GPUAtom> GPUAtoms;

        uint32_t atomOffset = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
            {
                bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex);

                if (FinalVisibilityInModelWorld == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                        {
                            GPUParticle GPUParticleObject;
                            GPUParticleObject.EntityId = ParticleObject.second.EntityId;
                            GPUParticleObject.ChainId = ParticleObject.second.ChainId;
                            GPUParticleObject.Index = ParticleObject.second.Index;

                            GPUParticleObject.AtomOffset = atomOffset;
                            GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();
                            memset(GPUParticleObject._padding, 0, sizeof(GPUParticleObject._padding));

                            GPUParticles.emplace_back(GPUParticleObject);

                            for (const auto& atom : ParticleObject.second.ListOfAtoms)
                            {
                                GPUAtom gpuAtom;
                                gpuAtom.EntityId = atom.EntityId;
                                gpuAtom.X = atom.X;
                                gpuAtom.Y = atom.Y;
                                gpuAtom.Z = atom.Z;
                                // gpuAtom.ColorR = atom.AtomColor.X;
                                // gpuAtom.ColorG = atom.AtomColor.Y;
                                // gpuAtom.ColorB = atom.AtomColor.Z;
                                gpuAtom.AtomColor = atom.AtomColor;

                                gpuAtom._padding1 = 0;

                                // Copy strings with padding
                                memset(gpuAtom.Name, 0, 8);
                                strncpy(gpuAtom.Name, atom.Name, 6);

                                memset(gpuAtom.ResName, 0, 8);
                                strncpy(gpuAtom.ResName, atom.ResName, 5);

                                memset(gpuAtom.Chain, 0, 8);
                                strncpy(gpuAtom.Chain, atom.Chain, 7);

                                GPUAtoms.push_back(gpuAtom);
                            }

                            atomOffset += ParticleObject.second.ListOfAtoms.size();
                        }
            }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticleSSBO);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUParticles.size() * sizeof(GPUParticle), GPUParticles.data());

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, AtomSSBO);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUAtoms.size() * sizeof(GPUAtom), GPUAtoms.data());

        // glDispatchCompute((MAX_PARTICLES + 255) / 256, 1, 1);
        // glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);















            const auto start_time3 = chrono::high_resolution_clock::now();

            size_t TotalSize = 0;
            for (const auto& UniformsBlocksSingle : UniformsBlocks)
                if (UniformsBlocksSingle.empty() == false)
                    TotalSize += UniformsBlocksSingle.size();
            UniformsBlocksTotal.reserve(TotalSize);

            for (const auto& UniformsBlocksSingle : UniformsBlocks)
                if (UniformsBlocksSingle.empty() == false)
                    //UniformsBlocksTotal.insert(UniformsBlocksTotal.end(), UniformsBlocksSingle.begin(), UniformsBlocksSingle.end());
                    UniformsBlocksTotal.insert(UniformsBlocksTotal.end(), std::make_move_iterator(UniformsBlocksSingle.begin()), std::make_move_iterator(UniformsBlocksSingle.end()));

            const auto stop_time3 = chrono::high_resolution_clock::now();
            ExecutionDurationTimeForCheckingPreparingParticles3 += chrono::duration(stop_time3 - start_time3);
            */













            //for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++)
            // {
            // FOR_EACH_SECTOR_IN_XYZ_ONLY
            // {
            //     if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
            //     {
            //         bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex);
            //
            //         //FinalVisibilityInModelWorld = CheckVisibility(FinalVisibilityInModelWorld);
            //
            //         if (FinalVisibilityInModelWorld == true)
            //             if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
            //             {
            //                 for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
            //                 {
            //                     if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
            //                     {
            //                         // Threads.emplace_back([&, ParticleSectorXIndex, this]()
            //                         // {
            //                             //DrawBonds(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);
            //
            //                             ParticleObject.second.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };
            //
            //                             for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
            //                             {
            //                                 // if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
            //                                 // {
            //                                 //     glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
            //                                 //     TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
            //                                 // }
            //
            //                                 RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool, ParticleSectorXIndex);
            //                             }
            //                     }
            //                 }
            //             }
            //     }
            // }
            // }


        }
        //}

        const auto stop_time = chrono::high_resolution_clock::now();
        ExecutionDurationTimeForTotalPreparingParticles += chrono::duration(stop_time - start_time);











        const auto start_time111 = chrono::high_resolution_clock::now();

        vector<GPUParticle> GPUParticles;
        GPUParticles.reserve(10'000'000);
        vector<GPUAtom> GPUAtoms;
        GPUAtoms.reserve(1000'000'000);

        uint32_t atomOffset = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (const bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex); FinalVisibilityInModelWorld == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                        {
                            GPUParticle GPUParticleObject;
                            GPUParticleObject.EntityId = ParticleObject.second.EntityId;
                            GPUParticleObject.ChainId = ParticleObject.second.ChainId;
                            GPUParticleObject.Index = ParticleObject.second.Index;

                            GPUParticleObject.AtomOffset = atomOffset;
                            GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();
                            //memset(GPUParticleObject._padding, 0, sizeof(GPUParticleObject._padding));
                            //memset(GPUParticleObject._padding, 0, 12);

                            GPUParticles.emplace_back(GPUParticleObject);

                            for (const auto& atom : ParticleObject.second.ListOfAtoms)
                            {
                                GPUAtom gpuAtom;
                                //gpuAtom.EntityId = atom.EntityId;
                                gpuAtom.X = atom.X;
                                gpuAtom.Y = atom.Y;
                                gpuAtom.Z = atom.Z;
                                // gpuAtom.ColorR = atom.AtomColor.X;
                                // gpuAtom.ColorG = atom.AtomColor.Y;
                                // gpuAtom.ColorB = atom.AtomColor.Z;

                                // gpuAtom.AtomColor = atom.AtomColor;
                                //ParticleObject.second.RandomParticleKindColor = GetColor<CellEngineAtom>(atom, ParticleObject.second, false);

                                // gpuAtom.ColorR = ParticleObject.second.RandomParticleKindColor.X;
                                // gpuAtom.ColorG = ParticleObject.second.RandomParticleKindColor.Y;
                                // gpuAtom.ColorB = ParticleObject.second.RandomParticleKindColor.Z;

                                gpuAtom.ColorR = static_cast<float>(ParticleObject.second.RandomParticleKindColor.X);
                                gpuAtom.ColorG = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Y);
                                gpuAtom.ColorB = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Z);

                                //gpuAtom._padding1[0] = 0.0f;
                                //gpuAtom._padding1[1] = 0.0f;
                                //memset(gpuAtom._padding, 0, 14);

                                // gpuAtom._padding1 = 0;
                                //
                                // // Copy strings with padding
                                // memset(gpuAtom.Name, 0, 8);
                                // strncpy(gpuAtom.Name, atom.Name, 6);
                                //
                                // memset(gpuAtom.ResName, 0, 8);
                                // strncpy(gpuAtom.ResName, atom.ResName, 5);
                                //
                                // memset(gpuAtom.Chain, 0, 8);
                                // strncpy(gpuAtom.Chain, atom.Chain, 7);

                                GPUAtoms.emplace_back(gpuAtom);
                            }

                            atomOffset += ParticleObject.second.ListOfAtoms.size();
                        }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);

        glUseProgram(ComputeShaderProgramPhong);

        glUniform3fv(glGetUniformLocation(ComputeShaderProgramPhong, "Center"), 1, Center);
        glUniformMatrix4fv(glGetUniformLocation(ComputeShaderProgramPhong, "ViewMatrix"), 1, GL_FALSE, ViewMatrix);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ParticleSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, AtomSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        const auto start_time112 = chrono::high_resolution_clock::now();

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticleSSBO);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUParticles.size() * sizeof(GPUParticle), GPUParticles.data());

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, AtomSSBO);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUAtoms.size() * sizeof(GPUAtom), GPUAtoms.data());

        const auto stop_time112 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory1 += chrono::duration(stop_time112 - start_time112);


        const auto start_time113 = chrono::high_resolution_clock::now();

        glDispatchCompute((GPUParticles.size() + 255) / 256, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        //glDispatchCompute((MAX_PARTICLES + 255) / 256, 1, 1);
        //glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);

        const auto stop_time113 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory2 += chrono::duration(stop_time113 - start_time113);

















        glUseProgram(ShaderProgramPhong);







        GLint ProjLoc = glGetUniformLocation(ShaderProgramPhong, "ProjectionMatrix");
        //glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(ProjectionMatrixGlobal));
        glUniformMatrix4fv(ProjLoc, 1, GL_FALSE, ProjectionMatrixGlobal);




























        //if (PressedRightMouseButton != 1)
        //    DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedAtomsList);

        const auto start_time1 = chrono::high_resolution_clock::now();

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        //glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        //glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, UniformsBlocks.size() * sizeof(UniformsBlock), UniformsBlocks.data());




        // UnsignedInt Offset = 0;
        // for (const auto& UniformsBlocksSingle : UniformsBlocks)
        // {
        //     glBufferSubData(GL_SHADER_STORAGE_BUFFER, Offset, UniformsBlocksSingle.size() * sizeof(UniformsBlock), UniformsBlocksSingle.data());
        //     Offset += UniformsBlocksSingle.size() * sizeof(UniformsBlock);
        //     LoggersManagerObject.Log(STREAM("S=" << UniformsBlocksSingle.size()));
        // }

        //glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, UniformsBlocksTotal.size() * sizeof(UniformsBlock), UniformsBlocksTotal.data());

        //glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, UniformsBlocks[0].size() * sizeof(UniformsBlock), UniformsBlocks[0].data());

        //LoggersManagerObject.Log(STREAM("S=" << UniformsBlocks.size()));
        //LoggersManagerObject.Log(STREAM("S=" << UniformsBlocksTotal.size()));
        //LoggersManagerObject.Log(STREAM("NX=" << CellEngineConfigDataObject.NumberOfParticlesSectorsInX));

        AtomGraphicsObject.RenderSubGraphicObject(0, GPUAtoms.size(), 0);

        const auto stop_time1 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory += chrono::duration(stop_time1 - start_time1);

        //AtomGraphicsObject.RenderSubGraphicObject(0, UniformsBlocks.size(), 0);
        //AtomGraphicsObject.RenderSubGraphicObject(0, UniformsBlocksTotal.size(), 0);
        //AtomGraphicsObject.RenderSubGraphicObject(0, UniformsBlocks[0].size(), 0);

        UniformsBlocks.clear();
        // for (auto& UniformBlock : UniformsBlocks)
        //     UniformBlock.clear();
    }
    CATCH("rendering full atom simulation space");
}

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>>& TemporaryRenderedAtomsList)
{
    try
    {
        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
        {
            UnsignedInt ChosenParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

            if (ChosenParticleCenterIndex > 0)
            {
                Particle ChosenParticleObject{};
                CellEngineAtom ChosenAtomObject{};
                if (ChosenParticleCenterIndex < TemporaryRenderedAtomsList.size())
                {
                    const UnsignedInt ParticleSectorXIndex = get<0>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);
                    const UnsignedInt ParticleSectorYIndex = get<1>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);
                    const UnsignedInt ParticleSectorZIndex = get<2>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);

                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.find(get<3>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.end())
                    {
                        if (get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]) > ParticleIter->second.ListOfAtoms.size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])));
                        else
                        {
                            ChosenParticleObject = ParticleIter->second;
                            ChosenAtomObject = ParticleIter->second.ListOfAtoms[get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])];
                        }
                    }
                    else
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(get<3>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])));
                }

                RenderObject(ChosenAtomObject, ChosenParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool, 1);

                PrintAtomDescriptionOnScreen(ChosenAtomObject, ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartCenterPoint()
{
    Center = CellEngineDataFileObjectPointer->GetCenterForAllParticles();

    LoggersManagerObject.Log(STREAM("CENTER OF CELL = " << Center.X() << " " << Center.Y() << " " << Center.Z()));
}

#endif