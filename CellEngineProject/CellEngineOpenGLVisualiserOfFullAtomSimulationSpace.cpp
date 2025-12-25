
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

        UniformsBlocks.clear();

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

                                GPUParticles.emplace_back(GPUParticleObject);

                                for (const auto& atom : ParticleObject.second.ListOfAtoms)
                                {
                                    GPUAtom gpuAtom;

                                    gpuAtom.X = atom.X;
                                    gpuAtom.Y = atom.Y;
                                    gpuAtom.Z = atom.Z;

                                    gpuAtom.ColorR = static_cast<float>(ParticleObject.second.RandomParticleKindColor.X);
                                    gpuAtom.ColorG = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Y);
                                    gpuAtom.ColorB = static_cast<float>(ParticleObject.second.RandomParticleKindColor.Z);

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

            //AtomGraphicsObject.RenderSubGraphicObject(0, UniformsBlocks.size(), 0);

            UniformsBlocks.clear();
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