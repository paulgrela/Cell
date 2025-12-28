
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
        lock_guard LockGuard{ RenderMenuAndFullAtomSimulationSpaceMutexObject };

        GLuint PartOfStencilBufferIndex[3];
        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        //UniformsBlocks.clear();

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFile::ChosenStructureMutexObject };

            //const auto start_time = chrono::high_resolution_clock::now();
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


            //}
            //}

            //const auto stop_time = chrono::high_resolution_clock::now();
            //ExecutionDurationTimeForTotalPreparingParticles += chrono::duration(stop_time - start_time);











            const auto start_time114 = chrono::high_resolution_clock::now();









            uint32_t ParticlesOffsetTotal = 0;
            uint32_t AtomOffsetTotal = 0;

            vector<GPUParticle> GPUParticles;
            GPUParticles.reserve(10'000'000);
            vector<GPUAtom> GPUAtoms;
            GPUAtoms.reserve(1000'000'000);




            vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt>> ParticlesSectorsToBeRendered;

            //NIECH SPAMIETUJE POCZATEK I KONIEC KAZDEGO PRZEDZIALU ATOMOW DO KOPIOWANIA
            uint32_t ParticleOffsetInSectors[40][40][40][20000];
            uint32_t AtomOffsetInSectors[40][40][40][20000];
            //uint32_t AtomOffsetInSectors[40][40][40];
            FOR_EACH_SECTOR_IN_XYZ_ONLY
            {
                // ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;
                // AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;

                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                    if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex))
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        {
                            ParticlesSectorsToBeRendered.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex);

                            UnsignedInt ParticleIndex = 0;
                            for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            {
                                ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] = ParticlesOffsetTotal;
                                ParticlesOffsetTotal++;

                                //LoggersManagerObject.Log(STREAM("P=" << ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex]));

                                AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] = AtomOffsetTotal;
                                //AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = AtomOffsetTotal;
                                AtomOffsetTotal += ParticleObject.second.ListOfAtoms.size();

                                ParticleIndex++;

                                //LoggersManagerObject.Log(STREAM("A=" << AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] << " S=" << ParticleObject.second.ListOfAtoms.size()));
                            }
                        }
            }

            const auto stop_time114 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);

            //LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time114, stop_time114, "Execution B has taken time = ","Execution in threads")));


            const auto start_time111 = chrono::high_resolution_clock::now();

            GPUParticles.resize(ParticlesOffsetTotal);
            GPUAtoms.resize(AtomOffsetTotal);


            //LoggersManagerObject.Log(STREAM("A=" << ParticlesSectorsToBeRendered.size()));
            //cout << "MAX THREAD " << omp_get_max_threads() << " Thread " << omp_get_thread_num() << " at level " << omp_get_level() << endl;


            // #pragma omp parallel for collapse(3) num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, NumberOfAllRenderedAtoms) schedule(dynamic)
            // FOR_EACH_SECTOR_IN_XYZ_ONLY

            vector<thread> Threads;

            UnsignedInt ParticlesSectorToBeRenderedIndex;
            // omp_set_nested(1);
            // omp_set_max_active_levels(2);
            // omp_set_dynamic(0);
            //#pragma omp parallel for num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered,NumberOfAllRenderedAtoms) private(ParticlesSectorToBeRenderedIndex) schedule(static)

            //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered,NumberOfAllRenderedAtoms) private(ParticlesSectorToBeRenderedIndex) schedule(static)

            //#pragma omp parallel for num_threads(2) default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)

            //#pragma omp parallel for default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)
            //#pragma omp parallel for default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)
            //#pragma omp parallel for

            // size_t totalSize = ParticlesSectorsToBeRendered.size();
            // size_t halfSize = (totalSize + 1) / 2;
            // #pragma omp parallel for num_threads(2) schedule(static, halfSize)

            #pragma omp parallel for num_threads(4) schedule(static)


            for (ParticlesSectorToBeRenderedIndex = 0; ParticlesSectorToBeRenderedIndex < ParticlesSectorsToBeRendered.size(); ParticlesSectorToBeRenderedIndex++)
            {
                //Threads.emplace_back([&, this]()
                //{

                const UnsignedInt ParticleSectorXIndex = get<0>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
                const UnsignedInt ParticleSectorYIndex = get<1>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
                const UnsignedInt ParticleSectorZIndex = get<2>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);

                    //cout << "MAX THREAD S " << omp_get_max_threads() << " Thread S " << omp_get_thread_num() << " at level S " << omp_get_level() << endl;

                    // if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                    //     if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex))
                    //         if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)


                    //this_thread::sleep_for(chrono::milliseconds(1));
                //});

                        {
                            UnsignedInt ParticleIndex = 0;
                            for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            {
                                GPUParticle GPUParticleObject;

                                GPUParticleObject.EntityId = ParticleObject.second.EntityId;
                                GPUParticleObject.ChainId = ParticleObject.second.ChainId;
                                GPUParticleObject.Index = ParticleObject.second.Index;

                                //GPUParticleObject.AtomOffset = AtomOffsetTotal;
                                //TU JEDEN ATOM_OFFSET DLA CZASTKI ATOM - A POWINIEN BYC JESZCZE ZAPAMIETANY DLA KAZDEJ CZASTKI Z OSOBNA
                                //BO W JEDNYM SEKTORZE MUSI BYC ZAPAMIETANYCH WIELE PRZESUNIEC DLA ATOMOW
                                GPUParticleObject.AtomOffset = AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex];
                                GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                                //LoggersManagerObject.Log(STREAM("P=" << ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex]));

                                //GPUParticles[ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + ParticleIndex] = GPUParticleObject;
                                GPUParticles[ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex]] = GPUParticleObject;

                                UnsignedInt AtomIndex = 0;
                                for (const auto& Atom : ParticleObject.second.ListOfAtoms)
                                {
                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = Atom.X;
                                    GPUAtomObject.Y = Atom.Y;
                                    GPUAtomObject.Z = Atom.Z;

                                    GPUAtomObject.ColorR = static_cast<float>(ParticleObject.second.RandomParticleKindColor.X);
                                    GPUAtomObject.ColorG = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Y);
                                    GPUAtomObject.ColorB = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Z);

                                    //LoggersManagerObject.Log(STREAM("A=" << AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + AtomIndex << " S=" << ParticleObject.second.ListOfAtoms.size()));

                                    GPUAtoms[AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] + AtomIndex] = GPUAtomObject;
                                    //GPUAtoms[AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + AtomIndex] = GPUAtomObject;
                                    AtomIndex++;
                                }
                                ParticleIndex++;
                            }
                        }
            }


            //for (auto& t : Threads)
            //    t.join();

            const auto stop_time111 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);


            //LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time111, stop_time111, "Execution B has taken time = ","Execution in threads")));

            //const auto start_time111 = chrono::high_resolution_clock::now();

            // FOR_EACH_SECTOR_IN_XYZ_ONLY
            //     if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
            //         if (const bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex); FinalVisibilityInModelWorld == true)
            //             if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
            //                 for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
            //                 {
            //                     GPUParticle GPUParticleObject;
            //                     GPUParticleObject.EntityId = ParticleObject.second.EntityId;
            //                     GPUParticleObject.ChainId = ParticleObject.second.ChainId;
            //                     GPUParticleObject.Index = ParticleObject.second.Index;
            //
            //                     GPUParticleObject.AtomOffset = AtomOffsetTotal;
            //                     GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();
            //
            //                     GPUParticles.emplace_back(GPUParticleObject);
            //
            //                     for (const auto& atom : ParticleObject.second.ListOfAtoms)
            //                     {
            //                         GPUAtom gpuAtom;
            //
            //                         gpuAtom.X = atom.X;
            //                         gpuAtom.Y = atom.Y;
            //                         gpuAtom.Z = atom.Z;
            //
            //                         gpuAtom.ColorR = static_cast<float>(ParticleObject.second.RandomParticleKindColor.X);
            //                         gpuAtom.ColorG = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Y);
            //                         gpuAtom.ColorB = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Z);
            //
            //                         GPUAtoms.emplace_back(gpuAtom);
            //                     }
            //
            //                     AtomOffsetTotal += ParticleObject.second.ListOfAtoms.size();
            //                 }
            // const auto stop_time111 = chrono::high_resolution_clock::now();
            //
            // ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);




            /*
            const auto start_time111 = chrono::high_resolution_clock::now();

            uint32_t AtomOffsetInSectors[40][40][40];
            vector<GPUParticle> GPUParticlesInSectors[40][40][40];
            vector<GPUAtom> GPUAtomsInSectors[40][40][40];

            FOR_EACH_SECTOR_IN_XYZ_ONLY
                AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;

            omp_set_nested(1);
            omp_set_max_active_levels(2);
            omp_set_dynamic(0);

            #pragma omp parallel for collapse(3) num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, AtomOffsetInSectors, GPUAtomsInSectors, GPUParticlesInSectors, NumberOfAllRenderedAtoms) schedule(dynamic)
            FOR_EACH_SECTOR_IN_XYZ_ONLY
                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                    if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex) == true)
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                            for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            {
                                GPUParticle GPUParticleObject;

                                GPUParticleObject.EntityId = ParticleObject.second.EntityId;
                                GPUParticleObject.ChainId = ParticleObject.second.ChainId;
                                GPUParticleObject.Index = ParticleObject.second.Index;

                                GPUParticleObject.AtomOffset = AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex];
                                GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                                GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUParticleObject);

                                for (const auto& Atom : ParticleObject.second.ListOfAtoms)
                                {
                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = Atom.X;
                                    GPUAtomObject.Y = Atom.Y;
                                    GPUAtomObject.Z = Atom.Z;

                                    GPUAtomObject.ColorR = static_cast<float>(ParticleObject.second.RandomParticleKindColor.X);
                                    GPUAtomObject.ColorG = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Y);
                                    GPUAtomObject.ColorB = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Z);

                                    GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUAtomObject);
                                }

                                AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] += ParticleObject.second.ListOfAtoms.size();
                            }

            const auto stop_time111 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);


            const auto start_time114 = chrono::high_resolution_clock::now();

            FOR_EACH_SECTOR_IN_XYZ_ONLY
                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                {
                    AtomOffsetTotal += AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex];
                    for (auto& GPUParticleInSectors : GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex])
                    {
                        GPUParticleInSectors.AtomOffset += AtomOffsetTotal;
                        GPUParticles.emplace_back(GPUParticleInSectors);
                    }
                    GPUAtoms.insert(GPUAtoms.end(), make_move_iterator(GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].begin()), make_move_iterator(GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].end()));
                }

            const auto stop_time114 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);
            */




            //lock_guard LockGuard{ RenderMenuAndFullAtomSimulationSpaceMutexObject };

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

            const auto stop_time113 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory2 += chrono::duration(stop_time113 - start_time113);



            glUseProgram(ShaderProgramPhong);

            glUniformMatrix4fv(glGetUniformLocation(ShaderProgramPhong, "ProjectionMatrix"), 1, GL_FALSE, ProjectionMatrixGlobal);



            //if (PressedRightMouseButton != 1)
            //    DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedAtomsList);

            const auto start_time1 = chrono::high_resolution_clock::now();

            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

            //glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, UniformsBlocks.size() * sizeof(UniformsBlock), UniformsBlocks.data());

            AtomGraphicsObject.RenderSubGraphicObject(0, GPUAtoms.size(), 0);

            const auto stop_time1 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory += chrono::duration(stop_time1 - start_time1);
        }
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