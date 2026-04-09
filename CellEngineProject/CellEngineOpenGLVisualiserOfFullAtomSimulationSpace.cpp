
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

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace1(const vmath::mat4& ViewMatrix, uint32_t& ParticlesOffsetTotal, uint32_t& AtomOffsetTotal, uint32_t& AtomLocalOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms, vector<GPUAtomLocal>& GPUAtomsLocal)
{
    try
    {
        const auto start_time111 = chrono::high_resolution_clock::now();

        uint32_t AtomTotalIndex = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles | views::values)
                            if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true)
                            {
                                UnsignedInt AtomIndex = 0;
                                UnsignedInt AtomLocalIndex = 0;

                                for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                {
                                    const auto& AtomObject = ParticleObject.ListOfAtoms[AtomObjectIndex];

                                    if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                    {
                                        GPUAtomLocal GPUAtomLocalObject;
                                        GPUAtomLocalObject.ParticleSectorXIndex = ParticleSectorXIndex;
                                        GPUAtomLocalObject.ParticleSectorYIndex = ParticleSectorYIndex;
                                        GPUAtomLocalObject.ParticleSectorZIndex = ParticleSectorZIndex;
                                        GPUAtomLocalObject.Index = ParticleObject.Index;
                                        GPUAtomLocalObject.AtomOffset = AtomLocalIndex + 1;
                                        GPUAtomsLocal[AtomLocalOffsetTotal] = std::move(GPUAtomLocalObject);
                                        AtomLocalOffsetTotal++;
                                        AtomLocalIndex++;
                                    }

                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = AtomObject.X;
                                    GPUAtomObject.Y = AtomObject.Y;
                                    GPUAtomObject.Z = AtomObject.Z;

                                    const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject, ChosenParticleObject.Index == ParticleObject.Index && ChosenAtomObjectIndex == AtomLocalIndex));

                                    GPUAtomObject.ColorR = ParticleColor.X();
                                    GPUAtomObject.ColorG = ParticleColor.Y();
                                    GPUAtomObject.ColorB = ParticleColor.Z();

                                    GPUAtoms[AtomTotalIndex] = std::move(GPUAtomObject);

                                    AtomIndex++;
                                    AtomTotalIndex++;
                                }


                                GPUParticle GPUParticleObject;

                                GPUParticleObject.AtomOffset = AtomOffsetTotal;
                                GPUParticleObject.AtomCount = AtomIndex;

                                GPUParticles[ParticlesOffsetTotal++] = std::move(GPUParticleObject);

                                AtomOffsetTotal += AtomIndex;
                            }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace2(const vmath::mat4& ViewMatrix, uint32_t& ParticlesOffsetTotal, uint32_t& AtomOffsetTotal, uint32_t& AtomLocalOffsetTotal)
{
    try
    {
        const auto start_time111 = chrono::high_resolution_clock::now();

        uint32_t AtomOffsetInSectors[40][40][40] = {};
        vector<GPUParticle> GPUParticlesInSectors[40][40][40];
        vector<GPUAtom> GPUAtomsInSectors[40][40][40];
        vector<GPUAtomLocal> GPUAtomsLocalInSectors[40][40][40];

        omp_set_nested(1);
        omp_set_max_active_levels(2);
        omp_set_dynamic(0);

        #pragma omp parallel for collapse(3) num_threads(128) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, AtomOffsetInSectors, GPUAtomsInSectors, GPUParticlesInSectors, GPUAtomsLocalInSectors) schedule(dynamic)
        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                            {
                                GPUParticle GPUParticleObject;

                                GPUParticleObject.AtomOffset = AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex];
                                GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                                GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUParticleObject);

                                UnsignedInt AtomIndex = 0;
                                UnsignedInt AtomLocalIndex = 0;

                                for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                {
                                    const auto& AtomObject = ParticleObject.second.ListOfAtoms[AtomObjectIndex];

                                    if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                    {
                                        GPUAtomLocal GPUAtomLocalObject;
                                        GPUAtomLocalObject.ParticleSectorXIndex = ParticleSectorXIndex;
                                        GPUAtomLocalObject.ParticleSectorYIndex = ParticleSectorYIndex;
                                        GPUAtomLocalObject.ParticleSectorZIndex = ParticleSectorZIndex;
                                        GPUAtomLocalObject.Index = ParticleObject.second.Index;
                                        GPUAtomLocalObject.AtomOffset = AtomLocalIndex++;
                                        GPUAtomsLocalInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUAtomLocalObject);
                                    }

                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = AtomObject.X;
                                    GPUAtomObject.Y = AtomObject.Y;
                                    GPUAtomObject.Z = AtomObject.Z;

                                    const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject.second, ChosenParticleObject.Index == ParticleObject.second.Index && ChosenAtomObjectIndex == AtomIndex));

                                    GPUAtomObject.ColorR = ParticleColor.X();
                                    GPUAtomObject.ColorG = ParticleColor.Y();
                                    GPUAtomObject.ColorB = ParticleColor.Z();

                                    GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUAtomObject);

                                    AtomIndex++;
                                }

                                AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] += AtomIndex;
                            }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);


        const auto start_time114 = chrono::high_resolution_clock::now();

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
            {
                for (auto& GPUParticleInSectors : GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex])
                {
                    GPUParticleInSectors.AtomOffset += AtomOffsetTotal;
                    GPUParticles[ParticlesOffsetTotal++] = std::move(GPUParticleInSectors);
                }

                memcpy(&GPUAtoms[AtomOffsetTotal], GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].data(), GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size() * sizeof(GPUAtom));

                if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                    memcpy(&GPUAtomsLocal[AtomLocalOffsetTotal], GPUAtomsLocalInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].data(), GPUAtomsLocalInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size() * sizeof(GPUAtomLocal));

                AtomOffsetTotal += GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size();
                AtomLocalOffsetTotal += GPUAtomsLocalInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size();
            }

        const auto stop_time114 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);
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
                    bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale);

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
                                        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                        {
                                            glStencilFunc(GL_ALWAYS, static_cast<uint8_t>((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                            TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                                        }

                                        RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, false, RenderObjectsBool);
                                        NumberOfAllRenderedAtoms++;
                                    }
                                }
                }
            }

            if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                glReadPixels(static_cast<GLint>(MousePositionLocal.s.X), static_cast<GLint>(static_cast<float>(Info.WindowHeight) - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
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
        //RenderSpace0(ViewMatrix);

        //RenderSpace1(ViewMatrix, ParticlesOffsetTotal, AtomOffsetTotal, AtomLocalOffsetTotal, GPUParticles, GPUAtoms, GPUAtomsLocal);
        RenderSpace2(ViewMatrix, ParticlesOffsetTotal, AtomOffsetTotal, AtomLocalOffsetTotal);

        // //LoggersManagerObject.Log(STREAM("P=" << GPUParticles.size() << " " << GPUAtoms.size()));
        //
        // NumberOfAllRenderedAtoms = AtomOffsetTotal;
        // NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = ParticlesOffsetTotal;
        //
        // glUseProgram(ComputeShaderProgramPhong);
        //
        // glUniform3fv(glGetUniformLocation(ComputeShaderProgramPhong, "Center"), 1, Center);
        // glUniformMatrix4fv(glGetUniformLocation(ComputeShaderProgramPhong, "ViewMatrix"), 1, GL_FALSE, ViewMatrix);
        //
        // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ParticleSSBO);
        // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, AtomSSBO);
        // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);
        //
        // const auto start_time112 = chrono::high_resolution_clock::now();
        //
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticleSSBO);
        // glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, ParticlesOffsetTotal * sizeof(GPUParticle), GPUParticles.data());
        //
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, AtomSSBO);
        // glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (AtomOffsetTotal - 1) * sizeof(GPUAtom), GPUAtoms.data());
        //
        // const auto stop_time112 = chrono::high_resolution_clock::now();
        //
        // ExecutionDurationTimeForCopyingParticlesToGraphicMemory1 += chrono::duration(stop_time112 - start_time112);
        //
        // const auto start_time113 = chrono::high_resolution_clock::now();
        //
        // glDispatchCompute((ParticlesOffsetTotal + 255) / 256, 1, 1);
        // glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        //
        // const auto stop_time113 = chrono::high_resolution_clock::now();
        //
        // ExecutionDurationTimeForCopyingParticlesToGraphicMemory2 += chrono::duration(stop_time113 - start_time113);
        //
        //                                                                                                                 glBindFramebuffer(GL_FRAMEBUFFER, FrameBufferObject);
        //
        // const vmath::vec3 BackgroundColor = CellEngineConfigDataObject.BackgroundColors[CellEngineConfigDataObject.ChosenBackgroundColor];
        // glClearColor(BackgroundColor.data[0], BackgroundColor.data[1], BackgroundColor.data[2], 0.0f);
        //
        //                                                                                                                 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //
        //                                                                                                                 constexpr GLuint ClearValue = 0xFFFFFFFF;
        //                                                                                                                 glClearTexImage(ScreenBufferInstanceTexture, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClearValue);
        //
        //                                                                                                                 glEnable(GL_DEPTH_TEST);
        //                                                                                                                 glDepthFunc(GL_LESS);
        //
        // glUseProgram(ShaderProgramPhong);
        //
        // glUniformMatrix4fv(glGetUniformLocation(ShaderProgramPhong, "ProjectionMatrix"), 1, GL_FALSE, ProjectionMatrixGlobal);
        //
        //                                                                                                                 glUniform1f(glGetUniformLocation(ShaderProgramPhong, "billboardDistance"), CellEngineConfigDataObject.Distance);
        //
        //                                                                                                                 vmath::vec2 screenSize(Info.WindowWidth, Info.WindowHeight);
        //                                                                                                                 glUniform2fv(glGetUniformLocation(ShaderProgramPhong, "screenSize"), 1, reinterpret_cast<float*>(&screenSize));
        //
        //
        // const auto start_time1 = chrono::high_resolution_clock::now();
        //
        // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);
        //
        // vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;
        //
        // if (CellEngineConfigDataObject.ViewPositionZ <= CellEngineConfigDataObject.Distance)
        //     AtomGraphicsObject.RenderSubGraphicObjectTriangles(0, AtomOffsetTotal, 0);
        // else
        //     AtomGraphicsObject.RenderSubGraphicObjectPoints(0, AtomOffsetTotal, 0);
        //
        //
        //                                                                                                                 glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        //
        //                                                                                                                 glReadBuffer(GL_COLOR_ATTACHMENT0);
        //
        //                                                                                                                 glBlitFramebuffer(0, 0, Info.WindowWidth, Info.WindowHeight, 0, 0, Info.WindowWidth, Info.WindowHeight, GL_COLOR_BUFFER_BIT,GL_LINEAR);
        //
        // const auto stop_time1 = chrono::high_resolution_clock::now();
        //
        // ExecutionDurationTimeForCopyingParticlesToGraphicMemory += chrono::duration(stop_time1 - start_time1);
        //
        //
        //
        //
        //                                                                                                                 if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        //                                                                                                                 {
        //                                                                                                                     uint32_t ClickedObjectID = 0xFFFFFFFF;
        //
        //                                                                                                                     glBindFramebuffer(GL_READ_FRAMEBUFFER, FrameBufferObject);
        //                                                                                                                     glReadBuffer(GL_COLOR_ATTACHMENT1);
        //
        //                                                                                                                     glReadPixels(static_cast<GLint>(MousePositionLocal.s.X), static_cast<GLint>(static_cast<float>(Info.WindowHeight) - MousePositionLocal.s.Y - 1), 1, 1, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClickedObjectID);
        //
        //                                                                                                                     glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
        //
        //                                                                                                                     if (PressedRightMouseButton != 1)
        //                                                                                                                         DrawChosenAtomUsingStencilBuffer1(ClickedObjectID);
        //
        //                                                                                                                     //LoggersManagerObject.Log(STREAM("C=" << CellEngineConfigDataObject.NumberOfStencilBufferLoops << " " << MousePositionLocal.s.X << " " << MousePositionLocal.s.Y << " " << Info.WindowWidth << " " << Info.WindowHeight << " " << ClickedObjectID));
        //                                                                                                                 }
    }
    CATCH("rendering full atom simulation space");
}

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer1(const GLuint ChosenParticleCenterIndex)
{
    try
    {
        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
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
                            ChosenAtomObjectIndex = GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset;
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
        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        {
            if (const UnsignedInt ChosenParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16); ChosenParticleCenterIndex > 0)
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

                RenderObject(ChosenAtomObject, ChosenParticleObject, ViewMatrix, false, false, false, true, RenderObjectsBool);

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