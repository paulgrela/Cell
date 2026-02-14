
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

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace1(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms, vector<GPUAtomLocal>& GPUAtomsLocal)
{
    try
    {
        const auto start_time111 = chrono::high_resolution_clock::now();

        uint32_t AtomOffsetTotal = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (const bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale); FinalVisibilityInModelWorld == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                            {
                                GPUParticle GPUParticleObject;

                                GPUParticleObject.AtomOffset = AtomOffsetTotal;
                                GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                                GPUParticles.emplace_back(GPUParticleObject);

                                UnsignedInt AtomIndex = 0;
                                for (const auto& AtomObject : ParticleObject.second.ListOfAtoms)
                                {
                                    if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                    {
                                        GPUAtomLocal GPUAtomLocalObject;
                                        GPUAtomLocalObject.ParticleSectorXIndex = ParticleSectorXIndex;
                                        GPUAtomLocalObject.ParticleSectorYIndex = ParticleSectorYIndex;
                                        GPUAtomLocalObject.ParticleSectorZIndex = ParticleSectorZIndex;
                                        GPUAtomLocalObject.Index = ParticleObject.second.Index;
                                        GPUAtomLocalObject.AtomOffset = AtomIndex++;
                                        GPUAtomsLocal.emplace_back(GPUAtomLocalObject);
                                    }

                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = AtomObject.X;
                                    GPUAtomObject.Y = AtomObject.Y;
                                    GPUAtomObject.Z = AtomObject.Z;

                                    const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject.second, ChosenParticleObject.Index == ParticleObject.second.Index && ChosenAtomObjectIndex == AtomIndex));

                                    GPUAtomObject.ColorR = ParticleColor.X();
                                    GPUAtomObject.ColorG = ParticleColor.Y();
                                    GPUAtomObject.ColorB = ParticleColor.Z();

                                    GPUAtoms.emplace_back(GPUAtomObject);
                                }

                                AtomOffsetTotal += ParticleObject.second.ListOfAtoms.size();
                            }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace0(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix)
{
    try
    {
        GLuint PartOfStencilBufferIndex[3];
        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFile::ChosenStructureMutexObject };

            FOR_EACH_SECTOR_IN_XYZ_ONLY
            {
                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                {
                    bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale);

                    FinalVisibilityInModelWorld = CheckVisibility(FinalVisibilityInModelWorld);

                    if (FinalVisibilityInModelWorld == true)
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                            for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                                if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                                {
                                    DrawBonds(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);

                                    ParticleObject.second.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };

                                    for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                    {
                                        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                        {
                                            glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                            TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                                        }

                                        RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                                    }
                                }
                }
            }

            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
        }

        if (PressedRightMouseButton != 1)
            DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedAtomsList);
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix)
{
    try
    {
        lock_guard LockGuard{ RenderMenuAndFullAtomSimulationSpaceMutexObject };

        //RenderSpace0(NumberOfAllRenderedAtoms, ViewMatrix);

        vector<GPUParticle> GPUParticles;
        GPUParticles.reserve(10'000'000);
        vector<GPUAtom> GPUAtoms;
        GPUAtoms.reserve(1000'000'000);
        vector<GPUAtomLocal> GPUAtomsLocal;
        GPUAtomsLocal.reserve(1000'000'000);

        RenderSpace1(NumberOfAllRenderedAtoms, ViewMatrix, GPUParticles, GPUAtoms, GPUAtomsLocal);
        //RenderSpace2(NumberOfAllRenderedAtoms, ViewMatrix, AtomOffsetTotal, GPUParticles, GPUAtoms);
        //RenderSpace3(NumberOfAllRenderedAtoms, ViewMatrix, AtomOffsetTotal, ParticlesOffsetTotal, GPUParticles, GPUAtoms);

        //LoggersManagerObject.Log(STREAM("P=" << GPUParticles.size() << " " << GPUAtoms.size()));

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






                                                                                                                        glBindFramebuffer(GL_FRAMEBUFFER, fbo);

                                                                                                                        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
                                                                                                                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                                                                                                                        GLuint ClearValue = 0xFFFFFFFF;
                                                                                                                        glClearTexImage(instanceTexture, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClearValue);

                                                                                                                        glEnable(GL_DEPTH_TEST);
                                                                                                                        glDepthFunc(GL_LESS);



        glUseProgram(ShaderProgramPhong);

        glUniformMatrix4fv(glGetUniformLocation(ShaderProgramPhong, "ProjectionMatrix"), 1, GL_FALSE, ProjectionMatrixGlobal);

        const auto start_time1 = chrono::high_resolution_clock::now();

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;


        AtomGraphicsObject.RenderSubGraphicObject(0, GPUAtoms.size(), 0);

                                                                                                                        glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
                                                                                                                        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

                                                                                                                        glReadBuffer(GL_COLOR_ATTACHMENT0);

                                                                                                                        glBlitFramebuffer(0, 0, TextureWidth, TextureHeight,0, 0, TextureWidth, TextureHeight, GL_COLOR_BUFFER_BIT,GL_LINEAR);

                                                                                                                        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        const auto stop_time1 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory += chrono::duration(stop_time1 - start_time1);

                            uint32_t ClickedObjectID;

                            glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
                            glReadBuffer(GL_COLOR_ATTACHMENT1);

                            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClickedObjectID);

                            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

                            if (PressedRightMouseButton != 1)
                                DrawChosenAtomUsingStencilBuffer1(ClickedObjectID, GPUAtomsLocal);
    }
    CATCH("rendering full atom simulation space");
}

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer1(const GLuint ChosenParticleCenterIndex, const vector<GPUAtomLocal>& GPUAtomsLocal)
{
    try
    {
        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
        {
            if (ChosenParticleCenterIndex > 0)
            {
                if (ChosenParticleCenterIndex < GPUAtomsLocal.size())
                {
                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorXIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorYIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorZIndex].Particles.find(GPUAtomsLocal[ChosenParticleCenterIndex].Index); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorXIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorYIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorZIndex].Particles.end())
                    {
                        if (GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset > ParticleIter->second.ListOfAtoms.size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(ChosenParticleCenterIndex));
                        else
                        {
                            ChosenParticleObject = ParticleIter->second;
                            ChosenAtomObject = ParticleIter->second.ListOfAtoms[GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset];
                            ChosenAtomObjectIndex = GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset + 1;
                        }
                    }
                    else
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(ChosenParticleCenterIndex));
                }

                PrintAtomDescriptionOnScreen(ChosenAtomObject, ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using buffer")
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

                RenderObject(ChosenAtomObject, ChosenParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

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