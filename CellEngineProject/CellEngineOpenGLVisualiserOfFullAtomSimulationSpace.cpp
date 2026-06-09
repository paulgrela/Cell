
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
                if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile || (CellEngineConfigDataObject.TypeOfFileToRead != CellEngineConfigData::TypesOfFileToRead::PDBFile && RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false) == true))
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

        #pragma omp parallel for collapse(3) num_threads(128) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, AtomOffsetInSectors, GPUAtomsInSectors, GPUParticlesInSectors, GPUAtomsLocalInSectors) schedule(dynamic)
        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile || (CellEngineConfigDataObject.TypeOfFileToRead != CellEngineConfigData::TypesOfFileToRead::PDBFile && RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false) == true))
                    for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                        if ((CellEngineConfigDataObject.ViewPositionZ <= CellEngineConfigDataObject.Distance + 700 && (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile || (CellEngineConfigDataObject.TypeOfFileToRead != CellEngineConfigData::TypesOfFileToRead::PDBFile && RenderObject(ParticleObject.second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false) == true))) || CellEngineConfigDataObject.ViewPositionZ > CellEngineConfigDataObject.Distance + 700)
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

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::FindAndDrawAllBondsBetweenAtoms(const vmath::mat4& ViewMatrix)
{
    try
    {
        if (CellEngineConfigDataObject.DrawBondsBetweenAtoms == true)
        {
            if (DrawBondsOnePragmaVersion == true)
            {
                FOR_EACH_SECTOR_IN_XYZ_ONLY
                    if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                        if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile || (CellEngineConfigDataObject.TypeOfFileToRead != CellEngineConfigData::TypesOfFileToRead::PDBFile && RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false) == true))
                            for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                                if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                                    FindAllBondsToDrawForParticle(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);
            }
            else
            {
                omp_set_max_active_levels(2);
                omp_set_dynamic(0);

                #pragma omp parallel for collapse(3) num_threads(16) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject) schedule(dynamic)
                FOR_EACH_SECTOR_IN_XYZ_ONLY
                    if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                        if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile || (CellEngineConfigDataObject.TypeOfFileToRead != CellEngineConfigData::TypesOfFileToRead::PDBFile && RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false) == true))
                            for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                                if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                                    FindAllBondsToDrawForParticle(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);
            }
        }
    }
    CATCH("drawing all bonds between atoms")
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(const vmath::mat4& ViewMatrix)
{
    try
    {
        if (CellEngineConfigDataObject.RenderCellWithParallelCPUComputing == false)
            RenderSpace1(ViewMatrix, ParticlesOffsetTotal, AtomOffsetTotal, AtomLocalOffsetTotal, GPUParticles, GPUAtoms, GPUAtomsLocal);
        else
            RenderSpace2(ViewMatrix, ParticlesOffsetTotal, AtomOffsetTotal, AtomLocalOffsetTotal);
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

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartCenterPoint()
{
    Center = CellEngineDataFileObjectPointer->GetCenterForAllParticles();

    LoggersManagerObject.Log(STREAM("CENTER OF CELL = " << Center.X() << " " << Center.Y() << " " << Center.Z()));
}

#endif