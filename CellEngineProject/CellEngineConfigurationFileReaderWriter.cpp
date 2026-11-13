
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <sb7color.h>

#include "StringUtils.h"
#include "ExceptionsMacro.h"
#include "DestinationPlatform.h"

#include "CellEngineUseful.h"
#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineConfigurationFileReaderWriter.h"

using namespace std;

void CellEngineConfigurationFileReaderWriter::ReadCellConfigurationFile(const char* ConfigFileNameParameter, const UnsignedInt ExecuteCellStateId)
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        ConfigFileName = ConfigFileNameParameter;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        LoggersManagerObject.Log(STREAM("Reading xml config file finished"));

        LoggersManagerObject.Log(STREAM("ExecuteCellStateId = " << ExecuteCellStateId));

        for (const ptree::value_type& MainConfigPropertyTreeElement : MainConfigPropertyTree.get_child("Settings"))
        {
            if (MainConfigPropertyTreeElement.first == "WindowsParameters")
            {
                #ifdef WINDOWS_PLATFORM
                CellEngineConfigDataObject.XTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMainWindows");
                CellEngineConfigDataObject.YTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMainWindows");
                CellEngineConfigDataObject.WidthMainWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMainWindows");
                CellEngineConfigDataObject.HeightMainWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMainWindows");
                CellEngineConfigDataObject.XTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMenuWindows");
                CellEngineConfigDataObject.YTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMenuWindows");
                CellEngineConfigDataObject.WidthMenuWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMenuWindows");
                CellEngineConfigDataObject.HeightMenuWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMenuWindows");
                CellEngineConfigDataObject.XTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondWindows");
                CellEngineConfigDataObject.YTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("YTopSecondWindows");
                CellEngineConfigDataObject.WidthSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondWindows");
                CellEngineConfigDataObject.HeightSecondWindow = MainConfigPropertyTreeElement.second.get<int>("HeightSecondWindows");
                #endif
                #ifdef UNIX_PLATFORM
                CellEngineConfigDataObject.XTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMainLinux");
                CellEngineConfigDataObject.YTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMainLinux");
                CellEngineConfigDataObject.WidthMainWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMainLinux");
                CellEngineConfigDataObject.HeightMainWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMainLinux");
                CellEngineConfigDataObject.XTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMenuLinux");
                CellEngineConfigDataObject.YTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMenuLinux");
                CellEngineConfigDataObject.WidthMenuWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMenuLinux");
                CellEngineConfigDataObject.HeightMenuWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMenuLinux");
                CellEngineConfigDataObject.XTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondLinux");
                CellEngineConfigDataObject.YTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("YTopSecondLinux");
                CellEngineConfigDataObject.WidthSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondLinux");
                CellEngineConfigDataObject.HeightSecondWindow = MainConfigPropertyTreeElement.second.get<int>("HeightSecondLinux");
                #endif
            }
            else
            if (MainConfigPropertyTreeElement.first == "Algorithm")
            {
                CellEngineConfigDataObject.MultiThreaded = MainConfigPropertyTreeElement.second.get<bool>("MultiThreaded");
                CellEngineConfigDataObject.SetProcessPriorityHighest = MainConfigPropertyTreeElement.second.get<bool>("SetProcessPriorityHighest");
                CellEngineConfigDataObject.PrintAtomDescriptionOnScreen = MainConfigPropertyTreeElement.second.get<bool>("PrintAtomDescriptionOnScreen");
                CellEngineConfigDataObject.LogParametersOfRenderingToFile = MainConfigPropertyTreeElement.second.get<bool>("LogParametersOfRenderingToFile");
                CellEngineConfigDataObject.RandomColorEngineObject = static_cast<CellEngineConfigData::RandomColorEngineTypes>(MainConfigPropertyTreeElement.second.get<UnsignedInt>("RandomColorEngineTypes"));
            }
            else
            if (MainConfigPropertyTreeElement.first == "Logger")
            {
                CellEngineConfigDataObject.PrintLogToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintToConsole");
                CellEngineConfigDataObject.PrintLogToFiles = MainConfigPropertyTreeElement.second.get<bool>("PrintToFiles");

                CellEngineConfigDataObject.PrintLogLineNumberToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToConsole");
                CellEngineConfigDataObject.PrintLogDateTimeToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToConsole");
                CellEngineConfigDataObject.PrintLogProcessIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToConsole");
                CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToConsole");
                CellEngineConfigDataObject.PrintLogThreadIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToConsole");

                CellEngineConfigDataObject.PrintLogLineNumberToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToFile");
                CellEngineConfigDataObject.PrintLogDateTimeToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToFile");
                CellEngineConfigDataObject.PrintLogProcessIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToFile");
                CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToFile");
                CellEngineConfigDataObject.PrintLogThreadIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToFile");

                CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile = MainConfigPropertyTreeElement.second.get<UnsignedInt>("MaximalNumberOfLinesInOneFile");

                CellEngineConfigDataObject.PrintLogToCommonFileWhenPrintLogToSpecialFile = MainConfigPropertyTreeElement.second.get<bool>("PrintLogToCommonFileWhenPrintLogToSpecialFile");
            }
            else
            if (MainConfigPropertyTreeElement.first == "CellsStates")
                for (const auto& CellStatePropertyTreeElement : MainConfigPropertyTree.get_child("Settings.CellsStates") | views::values)
                    if (CellStatePropertyTreeElement.get<UnsignedInt>("<xmlattr>.id") == ExecuteCellStateId)
                    {
                        CellEngineConfigDataObject.TypeOfSpace = static_cast<CellEngineConfigData::TypesOfSpace>(CellStatePropertyTreeElement.get<UnsignedInt>("TypeOfSpace"));

                        if (CellStatePropertyTreeElement.get_child_optional("DNAPaired"))
                            CellEngineConfigDataObject.DNAPaired = CellStatePropertyTreeElement.get<bool>("DNAPaired");

                        if (CellStatePropertyTreeElement.get_child_optional("ReverseReactantsAndProductsBecauseOfFormerErrorBool"))
                            CellEngineConfigDataObject.ReverseReactantsAndProductsBecauseOfFormerErrorBool = CellStatePropertyTreeElement.get<bool>("ReverseReactantsAndProductsBecauseOfFormerErrorBool");

                        if (CellStatePropertyTreeElement.get_child_optional("Radius1ForGenerationOfParticles"))
                            CellEngineConfigDataObject.Radius1ForGenerationOfParticles = CellStatePropertyTreeElement.get<UnsignedInt>("Radius1ForGenerationOfParticles");
                        if (CellStatePropertyTreeElement.get_child_optional("Radius1SizeForGenerationOfParticles"))
                            CellEngineConfigDataObject.Radius1SizeForGenerationOfParticles = CellStatePropertyTreeElement.get<UnsignedInt>("Radius1SizeForGenerationOfParticles");
                        if (CellStatePropertyTreeElement.get_child_optional("Radius2ForGenerationOfParticles"))
                            CellEngineConfigDataObject.Radius2ForGenerationOfParticles = CellStatePropertyTreeElement.get<UnsignedInt>("Radius2ForGenerationOfParticles");
                        if (CellStatePropertyTreeElement.get_child_optional("Radius2SizeForGenerationOfParticles"))
                            CellEngineConfigDataObject.Radius2SizeForGenerationOfParticles = CellStatePropertyTreeElement.get<UnsignedInt>("Radius2SizeForGenerationOfParticles");

                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfSimulationSpaceInEachDimension"))
                            CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfSimulationSpaceInEachDimension");
                        if (CellStatePropertyTreeElement.get_child_optional("DivisionFactorForSimulationSpace"))
                            CellEngineConfigDataObject.DivisionFactorForSimulationSpace = CellStatePropertyTreeElement.get<float>("DivisionFactorForSimulationSpace");
                        if (CellStatePropertyTreeElement.get_child_optional("DivisionFactorForReadingPositionsOfParticles"))
                            CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles = CellStatePropertyTreeElement.get<float>("DivisionFactorForReadingPositionsOfParticles");
                        if (CellStatePropertyTreeElement.get_child_optional("DivisionFactorForGeneratingPositionsOfParticles"))
                            CellEngineConfigDataObject.DivisionFactorForGeneratingPositionsOfParticles = CellStatePropertyTreeElement.get<float>("DivisionFactorForGeneratingPositionsOfParticles");

                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfStepsInSimulationOutside"))
                            CellEngineConfigDataObject.NumberOfStepsInSimulationOutside = CellStatePropertyTreeElement.get<int>("NumberOfStepsInSimulationOutside");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfStepsInSimulationInside"))
                            CellEngineConfigDataObject.NumberOfStepsInSimulationInside = CellStatePropertyTreeElement.get<int>("NumberOfStepsInSimulationInside");
                        if (CellStatePropertyTreeElement.get_child_optional("TypeOfSimulation"))
                            CellEngineConfigDataObject.TypeOfSimulation = static_cast<CellEngineConfigData::TypesOfSimulation>(CellStatePropertyTreeElement.get<UnsignedInt>("TypeOfSimulation"));

                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfXThreadsInSimulation"))
                            CellEngineConfigDataObject.NumberOfXThreadsInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfXThreadsInSimulation");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfYThreadsInSimulation"))
                            CellEngineConfigDataObject.NumberOfYThreadsInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfYThreadsInSimulation");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfZThreadsInSimulation"))
                            CellEngineConfigDataObject.NumberOfZThreadsInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfZThreadsInSimulation");

                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfXSectorsInOneThreadInSimulation"))
                            CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfXSectorsInOneThreadInSimulation");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfYSectorsInOneThreadInSimulation"))
                            CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfYSectorsInOneThreadInSimulation");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfZSectorsInOneThreadInSimulation"))
                            CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfZSectorsInOneThreadInSimulation");

                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfXInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfXInOneThreadInSimulationSpace");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfYInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfYInOneThreadInSimulationSpace");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfZInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfZInOneThreadInSimulationSpace");

                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfXInOneSectorInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfXInOneSectorInOneThreadInSimulationSpace");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfYInOneSectorInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfYInOneSectorInOneThreadInSimulationSpace");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfZInOneSectorInOneThreadInSimulationSpace"))
                            CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfZInOneSectorInOneThreadInSimulationSpace");

                        if (CellStatePropertyTreeElement.get_child_optional("StepToChangeSpaceDivisionForThreads"))
                            CellEngineConfigDataObject.StepToChangeSpaceDivisionForThreads = CellStatePropertyTreeElement.get<UnsignedInt>("StepToChangeSpaceDivisionForThreads");

                        if (CellStatePropertyTreeElement.get_child_optional("TypeOfExchangeOfParticlesBetweenThreads"))
                            CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads = static_cast<CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads>(CellStatePropertyTreeElement.get<UnsignedInt>("TypeOfExchangeOfParticlesBetweenThreads"));

                        if (CellStatePropertyTreeElement.get_child_optional("TypeOfSychronizationfSimulationExecutedInThreads"))
                            CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInThreads = static_cast<CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads>(CellStatePropertyTreeElement.get<UnsignedInt>("TypeOfSychronizationfSimulationExecutedInThreads"));

                        if (CellStatePropertyTreeElement.get_child_optional("TypeOfSychronizationfSimulationExecutedInMPIProcesses"))
                            CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInMPIProcesses = static_cast<CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInMPIProcesses>(CellStatePropertyTreeElement.get<UnsignedInt>("TypeOfSychronizationfSimulationExecutedInMPIProcesses"));

                        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                        {
                            CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStartXPos");
                            CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStartYPos");
                            CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStartZPos");
                            CellEngineConfigDataObject.SimulationSpaceSelectionStepX = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStepX");
                            CellEngineConfigDataObject.SimulationSpaceSelectionStepY = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStepY");
                            CellEngineConfigDataObject.SimulationSpaceSelectionStepZ = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionStepZ");
                            CellEngineConfigDataObject.SimulationSpaceSelectionSizeX = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionSizeX");
                            CellEngineConfigDataObject.SimulationSpaceSelectionSizeY = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionSizeY");
                            CellEngineConfigDataObject.SimulationSpaceSelectionSizeZ = CellStatePropertyTreeElement.get<UnsignedInt>("SimulationSpaceSelectionSizeZ");

                            CellEngineConfigDataObject.SizeOfBigPartOfTheCellMultiplyFactor = CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfBigPartOfTheCellMultiplyFactor");

                            CellEngineConfigDataObject.RadiusOfCellForDNA = CellStatePropertyTreeElement.get<double>("RadiusOfCellForDNA");
                        }

                        CellEngineConfigDataObject.CellStateFileName = CellStatePropertyTreeElement.get<string>("CellStateFileName");
                        CellEngineConfigDataObject.CellStateFileNameBackup = CellStatePropertyTreeElement.get<string>("CellStateFileNameBackup");

                        if (CellStatePropertyTreeElement.get_child_optional("UseNewerAutinNamesToAssociateWithIllinoisNames"))
                            CellEngineConfigDataObject.UseNewerAutinNamesToAssociateWithIllinoisNames = CellStatePropertyTreeElement.get<bool>("UseNewerAutinNamesToAssociateWithIllinoisNames");

                        if (CellStatePropertyTreeElement.get_child_optional("MixedFullAtomWithVoxelSpace"))
                            CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace = CellStatePropertyTreeElement.get<bool>("MixedFullAtomWithVoxelSpace");

                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfParticlesSectorsInX"))
                            CellEngineConfigDataObject.NumberOfParticlesSectorsInX = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfParticlesSectorsInX");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfParticlesSectorsInY"))
                            CellEngineConfigDataObject.NumberOfParticlesSectorsInY = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfParticlesSectorsInY");
                        if (CellStatePropertyTreeElement.get_child_optional("NumberOfParticlesSectorsInZ"))
                            CellEngineConfigDataObject.NumberOfParticlesSectorsInZ = CellStatePropertyTreeElement.get<UnsignedInt>("NumberOfParticlesSectorsInZ");

                        if (CellStatePropertyTreeElement.get_child_optional("ShiftCenterX"))
                            CellEngineConfigDataObject.ShiftCenterX = CellStatePropertyTreeElement.get<float>("ShiftCenterX");
                        if (CellStatePropertyTreeElement.get_child_optional("ShiftCenterY"))
                            CellEngineConfigDataObject.ShiftCenterY = CellStatePropertyTreeElement.get<float>("ShiftCenterY");
                        if (CellStatePropertyTreeElement.get_child_optional("ShiftCenterZ"))
                            CellEngineConfigDataObject.ShiftCenterZ = CellStatePropertyTreeElement.get<float>("ShiftCenterZ");

                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfParticlesSectorX"))
                            CellEngineConfigDataObject.SizeOfParticlesSectorX = CellStatePropertyTreeElement.get<float>("SizeOfParticlesSectorX");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfParticlesSectorY"))
                            CellEngineConfigDataObject.SizeOfParticlesSectorY = CellStatePropertyTreeElement.get<float>("SizeOfParticlesSectorY");
                        if (CellStatePropertyTreeElement.get_child_optional("SizeOfParticlesSectorZ"))
                            CellEngineConfigDataObject.SizeOfParticlesSectorZ = CellStatePropertyTreeElement.get<float>("SizeOfParticlesSectorZ");

                        if (CellStatePropertyTreeElement.get_child_optional("CellGenomePositionsFileName"))
                            CellEngineConfigDataObject.CellGenomePositionsFileName = CellStatePropertyTreeElement.get<string>("CellGenomePositionsFileName");
                        if (CellStatePropertyTreeElement.get_child_optional("CellGenomeSequenceFileName"))
                            CellEngineConfigDataObject.CellGenomeSequenceFileName = CellStatePropertyTreeElement.get<string>("CellGenomeSequenceFileName");
                        if (CellStatePropertyTreeElement.get_child_optional("CellGenomeSequenceFastaFileName"))
                            CellEngineConfigDataObject.CellGenomeSequenceFastaFileName = CellStatePropertyTreeElement.get<string>("CellGenomeSequenceFastaFileName");
                        if (CellStatePropertyTreeElement.get_child_optional("CellGenomeSequenceFastaOneFileName"))
                            CellEngineConfigDataObject.CellGenomeSequenceFastaOneFileName = CellStatePropertyTreeElement.get<string>("CellGenomeSequenceFastaOneFileName");

                        if (CellStatePropertyTreeElement.get_child_optional("FullAtomMPIParallelProcessesExecution"))
                            CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution = CellStatePropertyTreeElement.get<bool>("FullAtomMPIParallelProcessesExecution");

                        if (CellStatePropertyTreeElement.get_child_optional("OpenGLGraphicsSwitchedOff"))
                            CellEngineConfigDataObject.OpenGLGraphicsSwitchedOff = CellStatePropertyTreeElement.get<bool>("OpenGLGraphicsSwitchedOff");

                        CellEngineConfigDataObject.ChosenStructureIndex = CellStatePropertyTreeElement.get<UnsignedInt>("ChosenStructureIndex");

                        CellEngineConfigDataObject.SpecularPower = CellStatePropertyTreeElement.get<float>("SpecularPower");
                        CellEngineConfigDataObject.SpecularAlbedo = CellStatePropertyTreeElement.get<float>("SpecularAlbedo");
                        CellEngineConfigDataObject.MakeColorsTypeObject = static_cast<CellEngineConfigData::MakeColorsType>(CellStatePropertyTreeElement.get<UnsignedInt>("MakeColorsType"));

                        CellEngineConfigDataObject.DrawBondsBetweenAtoms = CellStatePropertyTreeElement.get<bool>("DrawBondsBetweenAtoms");

                        CellEngineConfigDataObject.CutZ = CellStatePropertyTreeElement.get<float>("CutZ");
                        CellEngineConfigDataObject.Distance = CellStatePropertyTreeElement.get<float>("Distance");
                        CellEngineConfigDataObject.LoadOfAtomsStep = CellStatePropertyTreeElement.get<UnsignedInt>("LoadOfAtomsStep");

                        CellEngineConfigDataObject.XLowToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("XLowToDrawInAtomScale");
                        CellEngineConfigDataObject.XHighToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("XHighToDrawInAtomScale");
                        CellEngineConfigDataObject.YLowToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("YLowToDrawInAtomScale");
                        CellEngineConfigDataObject.YHighToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("YHighToDrawInAtomScale");
                        CellEngineConfigDataObject.ZLowToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("ZLowToDrawInAtomScale");
                        CellEngineConfigDataObject.ZHighToDrawInAtomScale = CellStatePropertyTreeElement.get<float>("ZHighToDrawInAtomScale");

                        CellEngineConfigDataObject.SizeOfAtomsDrawingTypesObject = static_cast<CellEngineConfigData::SizeOfAtomsDrawingTypes>(CellStatePropertyTreeElement.get<UnsignedInt>("SizeOfAtomsDrawingTypes"));

                        CellEngineConfigDataObject.SizeOfAtomChangeStep = CellStatePropertyTreeElement.get<float>("SizeOfAtomChangeStep");
                        CellEngineConfigDataObject.SizeOfAtomX = CellStatePropertyTreeElement.get<float>("SizeOfAtomX");
                        CellEngineConfigDataObject.SizeOfAtomY = CellStatePropertyTreeElement.get<float>("SizeOfAtomY");
                        CellEngineConfigDataObject.SizeOfAtomZ = CellStatePropertyTreeElement.get<float>("SizeOfAtomZ");
                        CellEngineConfigDataObject.MainSizeOfAtom = 1.0;

                        CellEngineConfigDataObject.RotationAngle1 = CellStatePropertyTreeElement.get<float>("RotationAngle1");
                        CellEngineConfigDataObject.RotationAngle2 = CellStatePropertyTreeElement.get<float>("RotationAngle2");
                        CellEngineConfigDataObject.RotationAngle3 = CellStatePropertyTreeElement.get<float>("RotationAngle3");
                        CellEngineConfigDataObject.ViewPositionX = CellStatePropertyTreeElement.get<float>("ViewPositionX");
                        CellEngineConfigDataObject.ViewPositionY = CellStatePropertyTreeElement.get<float>("ViewPositionY");
                        CellEngineConfigDataObject.ViewPositionZ = CellStatePropertyTreeElement.get<float>("ViewPositionZ");
                        CellEngineConfigDataObject.CameraXPosition = CellStatePropertyTreeElement.get<float>("CameraXPosition");
                        CellEngineConfigDataObject.CameraYPosition = CellStatePropertyTreeElement.get<float>("CameraYPosition");
                        CellEngineConfigDataObject.CameraZPosition = CellStatePropertyTreeElement.get<float>("CameraZPosition");
                        CellEngineConfigDataObject.CameraXMoveShortStep = CellStatePropertyTreeElement.get<float>("CameraXMoveShortStep");
                        CellEngineConfigDataObject.CameraYMoveShortStep = CellStatePropertyTreeElement.get<float>("CameraYMoveShortStep");
                        CellEngineConfigDataObject.CameraZMoveShortStep = CellStatePropertyTreeElement.get<float>("CameraZMoveShortStep");
                        CellEngineConfigDataObject.CameraXMoveLongStep = CellStatePropertyTreeElement.get<float>("CameraXMoveLongStep");
                        CellEngineConfigDataObject.CameraYMoveLongStep = CellStatePropertyTreeElement.get<float>("CameraYMoveLongStep");
                        CellEngineConfigDataObject.CameraZMoveLongStep = CellStatePropertyTreeElement.get<float>("CameraZMoveLongStep");
                        CellEngineConfigDataObject.ViewPositionX = CellStatePropertyTreeElement.get<float>("ViewPositionX");
                        CellEngineConfigDataObject.ViewPositionY = CellStatePropertyTreeElement.get<float>("ViewPositionY");
                        CellEngineConfigDataObject.ViewPositionZ = CellStatePropertyTreeElement.get<float>("ViewPositionZ");
                        CellEngineConfigDataObject.ViewXMoveShortStep = CellStatePropertyTreeElement.get<float>("ViewXMoveShortStep");
                        CellEngineConfigDataObject.ViewYMoveShortStep = CellStatePropertyTreeElement.get<float>("ViewYMoveShortStep");
                        CellEngineConfigDataObject.ViewZMoveShortStep = CellStatePropertyTreeElement.get<float>("ViewZMoveShortStep");
                        CellEngineConfigDataObject.ViewXMoveLongStep = CellStatePropertyTreeElement.get<float>("ViewXMoveLongStep");
                        CellEngineConfigDataObject.ViewYMoveLongStep = CellStatePropertyTreeElement.get<float>("ViewYMoveLongStep");
                        CellEngineConfigDataObject.ViewZMoveLongStep = CellStatePropertyTreeElement.get<float>("ViewZMoveLongStep");

                        CellEngineConfigDataObject.ViewChangeUsingLongStep = CellStatePropertyTreeElement.get<bool>("ViewChangeUsingLongStep");
                        CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = CellStatePropertyTreeElement.get<bool>("AutomaticChangeOfLoadAtomsStep");

                        CellEngineConfigDataObject.BackgroundColors[1] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.get<string>("BackgroundColor1")));
                        CellEngineConfigDataObject.BackgroundColors[2] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.get<string>("BackgroundColor2")));
                        CellEngineConfigDataObject.BackgroundColors[3] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.get<string>("BackgroundColor3")));
                        CellEngineConfigDataObject.ChosenBackgroundColor = CellStatePropertyTreeElement.get<UnsignedInt>("ChosenBackgroundColor");

                        ParticlesKindsManagerObject.AtomsKindsGraphicData.clear();
                        for (const auto& AtomKindPropertyTreeElement : CellStatePropertyTreeElement.get_child("Atoms") | views::values)
                        {
                            AtomKindGraphicData AtomKindObject;

                            LoggersManagerObject.Log(STREAM("Atom Kind Name = " << AtomKindPropertyTreeElement.get<string>("Name")));

                            AtomKindObject.Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.get<string>("Color"))));
                            AtomKindObject.ColorVmathVec3 = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.get<string>("Color")));
                            AtomKindObject.SizeX = AtomKindPropertyTreeElement.get<float>("SizeX");
                            AtomKindObject.SizeY = AtomKindPropertyTreeElement.get<float>("SizeY");
                            AtomKindObject.SizeZ = AtomKindPropertyTreeElement.get<float>("SizeZ");
                            AtomKindObject.Name = AtomKindPropertyTreeElement.get<string>("Name");

                            ParticlesKindsManagerObject.AtomsKindsGraphicData.emplace_back(AtomKindObject);
                        }

                        ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.clear();
                        for (const auto& ParticleKindPropertyTreeElement : CellStatePropertyTreeElement.get_child("Particles") | views::values)
                        {
                            ParticleKindGraphicData ParticleKindGraphicObject;

                            LoggersManagerObject.Log(STREAM("ParticleKind Name = " << ParticleKindPropertyTreeElement.get<string>("Name") << " ID = " << ParticleKindPropertyTreeElement.get<UnsignedInt>("<xmlattr>.id") << " COLOR = " << ParticleKindPropertyTreeElement.get<string>("Color")));

                            ParticleKindGraphicObject.NameFromXML = ParticleKindPropertyTreeElement.get<string>("Name");
                            ParticleKindGraphicObject.NameFromXML == "DNA" ? CellEngineConfigDataObject.DNAIdentifier = ParticleKindPropertyTreeElement.get<UnsignedInt>("<xmlattr>.id") : 0;
                            ParticleKindGraphicObject.NameFromXML == "RNA" ? CellEngineConfigDataObject.RNAIdentifier = ParticleKindPropertyTreeElement.get<UnsignedInt>("<xmlattr>.id") : 0;
                            CellEngineConfigDataObject.ProteinInBuildingProcessIdentifier = CellEngineConfigDataObject.RNAIdentifier + 1;
                            ParticleKindGraphicObject.ParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(ParticleKindPropertyTreeElement.get<string>("Color"))));
                            ParticleKindGraphicObject.Visible = ParticleKindPropertyTreeElement.get<bool>("Visible");
                            ParticleKindGraphicObject.SizeX = ParticleKindPropertyTreeElement.get<float>("SizeX");
                            ParticleKindGraphicObject.SizeY = ParticleKindPropertyTreeElement.get<float>("SizeY");
                            ParticleKindGraphicObject.SizeZ = ParticleKindPropertyTreeElement.get<float>("SizeZ");
                            ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML[ParticleKindPropertyTreeElement.get<UnsignedInt>("<xmlattr>.id")] = ParticleKindGraphicObject;
                        }
                    }
        }
    }
    CATCH("reading cell configuration xml file")
}

void CellEngineConfigurationFileReaderWriter::SaveTestStatisticsToFile(const UnsignedInt ExecuteCellStateId) const
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        for (auto& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.Tests") | views::values)
            if (ExecuteCellStateId == TestPropertyTreeElement.get<UnsignedInt>("<xmlattr>.id"))
            {
                ptree& TestPropertyTreeElementToWriteInFile = TestPropertyTreeElement;

                TestPropertyTreeElementToWriteInFile.put("CellStateFileName", CellEngineConfigDataObject.CellStateFileName);

                TestPropertyTreeElementToWriteInFile.put("ChosenStructureIndex", CellEngineConfigDataObject.ChosenStructureIndex);
            }

        std::ofstream OutputConfigFile(ConfigFileName);
        write_xml(OutputConfigFile, MainConfigPropertyTree, boost::property_tree::xml_writer_make_settings<string>('\t', 1));
    }
    CATCH("saving test statistics configuration xml file")
}
