
#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"

#include "CellEngineUseful.h"
#include "CellEngineConstants.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool DrawError = false;

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetStartPositions()
{
    return { SelectionStartXPos, SelectionStartYPos, SelectionStartZPos };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetSteps()
{
    return { SelectionStepX, SelectionStepY, SelectionStepZ };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetSizes()
{
    return { SelectionSizeX, SelectionSizeY, SelectionSizeZ };
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetVoxelSpaceSelection(const UnsignedInt SelectionStartXParam, const UnsignedInt SelectionStartYParam, const UnsignedInt SelectionStartZParam, const UnsignedInt SelectionStepXParam, const UnsignedInt SelectionStepYParam, const UnsignedInt SelectionStepZParam, const UnsignedInt SelectionSizeXParam, const UnsignedInt SelectionSizeYParam, const UnsignedInt SelectionSizeZParam)
{
    SelectionStartXPos = SelectionStartXParam, SelectionStartYPos = SelectionStartYParam, SelectionStartZPos = SelectionStartZParam;
    SelectionStepX = SelectionStepXParam, SelectionStepY = SelectionStepYParam, SelectionStepZ = SelectionStepZParam;
    SelectionSizeX = SelectionSizeXParam, SelectionSizeY = SelectionSizeYParam, SelectionSizeZ = SelectionSizeZParam;
}

CellEngineOpenGLVisualiserOfVoxelSimulationSpace::CellEngineOpenGLVisualiserOfVoxelSimulationSpace()
{
    SetVoxelSpaceSelection(CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, CellEngineConfigDataObject.SimulationSpaceSelectionStepX, CellEngineConfigDataObject.SimulationSpaceSelectionStepY, CellEngineConfigDataObject.SimulationSpaceSelectionStepZ, CellEngineConfigDataObject.SimulationSpaceSelectionSizeX, CellEngineConfigDataObject.SimulationSpaceSelectionSizeY, CellEngineConfigDataObject.SimulationSpaceSelectionSizeZ);
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SaveVoxelPositionChosenByMouse()
{
    SelectionStartXPos = SaveXMousePosition;
    SelectionStartYPos = SaveYMousePosition;
    SelectionStartZPos = SaveZMousePosition;
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetSaveXYZPositions(const UnsignedInt SaveXParam, const UnsignedInt SaveYParam, const UnsignedInt SaveZParam)
{
    SaveXMousePosition = SaveXParam;
    SaveYMousePosition = SaveYParam;
    SaveZMousePosition = SaveZParam;
}

static inline float CovertToGraphicsCoordinateSelected(const UnsignedInt StartParam, const UnsignedInt SpaceParam, const UnsignedInt SizeParam)
{
    return ((static_cast<float>((StartParam + SizeParam) - SpaceParam) - static_cast<float>(SizeParam) / 2) * 4);
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::ConvertAtomPosToGraphicCoordinate(CellEngineAtom& CellEngineAtomObjectParam, const UnsignedInt StartXParam, const UnsignedInt StartYParam, const UnsignedInt StartZParam, const UnsignedInt SpaceXParam, const UnsignedInt SpaceYParam, const UnsignedInt SpaceZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam) const
{
    if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
        CellEngineAtomObjectParam.SetAtomPositionsData(CovertToGraphicsCoordinateSelected(StartXParam, SpaceXParam, SizeXParam), CovertToGraphicsCoordinateSelected(StartYParam, SpaceYParam, SizeYParam), CovertToGraphicsCoordinateSelected(StartZParam, SpaceZParam, SizeZParam));
    else
    if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
        CellEngineAtomObjectParam.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceXParam), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceYParam), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceZParam));
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetParticleParametersToDraw(CellEngineAtom& TempAtomObject, Particle& ParticleObject)
{
    try
    {
        TempAtomObject.EntityId = ParticleObject.EntityId;
        auto ParticleKindObject = ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId);
        TempAtomObject.AtomColor = ParticleKindObject.AtomColor;
        // TempAtomObject.ParticleColor = ParticleKindObject.ParticleColor;
        // TempAtomObject.UniqueParticleColor = ParticleObject.UniqueColor;
        // TempAtomObject.RandomParticleKindColor = (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true ? CellEngineColorsObject.GetDNAorRNAColor(TempAtomObject.EntityId, ParticleObject.ChainId) : ParticleKindObject.RandomParticleColor);
        ParticleObject.RandomParticleKindColor = (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true ? CellEngineColorsObject.GetDNAorRNAColor(TempAtomObject.EntityId, ParticleObject.ChainId) : ParticleKindObject.RandomParticleColor);

        if (CellEngineConfigDataObject.DNAPaired == true)
            if (CellEngineUseful::IsDNA(TempAtomObject.EntityId) == true && ((CellEngineConfigDataObject.GenomeReadFromFile == false && ParticleObject.GenomeIndex == 0) || (CellEngineConfigDataObject.GenomeReadFromFile == true && (ParticleObject.Prev != nullptr && ParticleObject.Next != nullptr && ParticleObject.PairedNucleotidePtr == nullptr))))
                //TempAtomObject.RandomParticleKindColor = TempAtomObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Purple));
                //TempAtomObject.AtomColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Purple));
                ParticleObject.RandomParticleKindColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Purple));

        if (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true && ((CellEngineConfigDataObject.GenomeReadFromFile == false && ParticleObject.GenomeIndex == 0) || (CellEngineConfigDataObject.GenomeReadFromFile == true && (ParticleObject.Prev == nullptr || ParticleObject.Next == nullptr))))
            //TempAtomObject.RandomParticleKindColor = TempAtomObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Blue));
            //TempAtomObject.AtomColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Blue));
            ParticleObject.RandomParticleKindColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Blue));

        if (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true && (CellEngineConfigDataObject.GenomeReadFromFile == true))
        {
            ParticleObject.GenomeIndexPrev = (ParticleObject.Prev != nullptr ? ParticleObject.Prev->GenomeIndex : 0);
            ParticleObject.GenomeIndexNext = (ParticleObject.Next != nullptr ? ParticleObject.Next->GenomeIndex : 0);
        }
        ParticleObject.Nucleotide = ((CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.ChainId) : '0');
    }
    CATCH("setting particle parameters to draw")
};

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSelectedSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList, UnsignedInt StencilBufferLoopCounter)
{
    try
    {
        for (UnsignedInt PosX = XStartParam; PosX < XStartParam + XSizeParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam; PosY < YStartParam + YSizeParam; PosY += YStepParam)
                for (UnsignedInt PosZ = ZStartParam; PosZ < ZStartParam + ZSizeParam; PosZ += ZStepParam)
                    if (PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                    {
                        SimulationSpaceVoxel SimulationSpaceVoxelObject = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, PosZ);
                        if (SimulationSpaceVoxelObject != CellEngineParticlesVoxelsOperations::GetZeroSimulationSpaceVoxel())
                        {
                            if (auto FoundParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(SimulationSpaceVoxelObject); FoundParticleIter != CellEngineDataFileObjectPointer->GetParticleEnd())
                            {
                                Particle& ParticleObject = FoundParticleIter->second;

                                if (DrawEmptyVoxels == true || (DrawEmptyVoxels == false && SimulationSpaceVoxelObject != 0 && ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true))
                                {
                                    ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                                    if (DrawEmptyVoxels == false || (DrawEmptyVoxels == true && SimulationSpaceVoxelObject != 0))
                                        SetParticleParametersToDraw(TempAtomObject, ParticleObject);
                                    // else
                                    // if (DrawEmptyVoxels == true)
                                    //     //TempAtomObject.AtomColor = TempAtomObject.ParticleColor = TempAtomObject.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::DeepSkyBlue));
                                    //     TempAtomObject.AtomColor = ParticleObject.ParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::DeepSkyBlue));

                                    if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                    {
                                        glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedVoxelsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                        TemporaryRenderedVoxelsList.emplace_back(TemporaryRenderedVoxel{ TempAtomObject, ParticleObject, PosX, PosY, PosZ });
                                    }

                                    RenderObject(TempAtomObject, ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                                }
                            }
                            else
                            if (DrawError ==  true)
                                LoggersManagerObject.LogError(STREAM("Try to draw the particle from not existing index = " << SimulationSpaceVoxelObject));
                        }
                        else
                        if (DrawEmptyVoxels == true)
                        {
                            Particle ParticleObject;
                            TempAtomObject.AtomColor = ParticleObject.ParticleColor = ParticleObject.UniqueParticleColor = ParticleObject.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::DeepSkyBlue));

                            ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                            RenderObject(TempAtomObject, ParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                        }
                    }
    }
    CATCH("rendering selected voxel simulation space")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, vmath::mat4& ViewMatrix)
{
    try
    {
        conditional_lock_guard<recursive_mutex> LockGuardCond(CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads, &CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderMenuAndVoxelSimulationSpaceMutexObject);

        GLuint PartOfStencilBufferIndex[3];

        CellEngineAtom TempAtomObject;

        NumberOfAllRenderedAtoms = 0;

        std::vector<TemporaryRenderedVoxel> TemporaryRenderedVoxelsList;

        CellEngineConfigDataObject.LoadOfAtomsStep > 10 ? CellEngineConfigDataObject.LoadOfAtomsStep = 4 : 1;

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedVoxelsList.clear();

            if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
                for (UnsignedInt PosX = SelectionStartXPos; PosX < SelectionSizeX; PosX += SelectionStepX)
                    for (UnsignedInt PosY = SelectionStartYPos; PosY < SelectionSizeY; PosY += SelectionStepY)
                        for (UnsignedInt PosZ = SelectionStartZPos; PosZ < SelectionSizeZ; PosZ += SelectionStepZ)
                        {
                            TempAtomObject.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosX), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosY), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosZ));

                            if (RenderObject(TempAtomObject, Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                            {
                                NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;
                                RenderSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);
                            }
                        }
            else
            {
                UnsignedInt SubStartPos = 0;
                if (CellEngineConfigDataObject.SelectedSpaceStartParametersDrawTypesObject == CellEngineConfigData::SelectedSpaceStartParametersDrawTypes::DrawFromCenter)
                    SubStartPos = SelectionSizeX / 2;

                RenderSelectedSpace(SelectionStartXPos - SubStartPos, SelectionStartYPos - SubStartPos, SelectionStartZPos - SubStartPos, SelectionStepX, SelectionStepY, SelectionStepZ, SelectionSizeX, SelectionSizeY, SelectionSizeY, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);
            }

            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
        }

        if (PressedRightMouseButton != 1)
            DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedVoxelsList);
    }
    CATCH("rendering voxel simulation space");
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList)
{
    try
    {
        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
        {
            UnsignedInt ChosenVoxelIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

            if (ChosenVoxelIndex > 0)
            {
                TemporaryRenderedVoxel ChosenObject{};

                if (ChosenVoxelIndex > TemporaryRenderedVoxelsList.size())
                    throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + std::to_string(ChosenVoxelIndex) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(TemporaryRenderedVoxelsList.size()));
                else
                    ChosenObject = TemporaryRenderedVoxelsList[ChosenVoxelIndex];

                SetSaveXYZPositions(ChosenObject.X, ChosenObject.Y, ChosenObject.Z);

                RenderObject(ChosenObject.CellEngineAtomObject, Particle(), ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

                PrintAtomDescriptionOnScreen(ChosenObject.CellEngineAtomObject, ChosenObject.CellEngineParticleObject);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetStartCenterPoint()
{
    Center = { 0.0f, 0.0f, 0.0f };
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetMemoryForBondsBetweenAtomsToDraw()
{
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::DrawBondsForParticlesCenters(std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix)
{
}