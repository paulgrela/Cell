
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, vmath::mat4& ViewMatrix)
{
    try
    {
        GLuint PartOfStencilBufferIndex[3];
        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFile::ChosenStructureMutexObject };

            FOR_EACH_PARTICLE_IN_XYZ_ONLY
            {
                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].empty() == false)
                {
                    bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale);

                    if (FinalVisibilityInModelWorld == true)
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        {
                            for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex])
                            {
                                if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters)
                                    glStencilFunc(GL_ALWAYS, uint8_t((NumberOfAllRenderedAtoms >> (8 * StencilBufferLoopCounter))), -1);

                                if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                                {
                                    NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                                    //DrawBonds(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex], BondsBetweenAtomsToDraw[ParticlesCenterObject.AtomIndex], CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);

                                    for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                    {
                                        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                            if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale)
                                            {
                                                glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                                TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                                            }

                                        RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                                    }
                                }
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

                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].find(get<3>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].end())
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
    //Center = CellEngineDataFile::GetCenter(CellEngineDataFileObjectPointer->GetParticlesCenters());
    Center = CellEngineDataFileObjectPointer->GetCenterForAllParticles();

    LoggersManagerObject.Log(STREAM("CENTER OF CELL = " << Center.X() << " " << Center.Y() << " " << Center.Z()));
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetMemoryForBondsBetweenAtomsToDraw()
{
    //BondsBetweenAtomsToDraw.resize(CellEngineDataFileObjectPointer->GetParticlesCenters().size());
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawBondsForParticlesCenters(std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix)
{
    //CellEngineOpenGLVisualiser::DrawBonds(CellEngineDataFileObjectPointer->GetParticlesCenters(), BondsToDraw, CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters, ViewMatrix);
}