#pragma once

#ifndef _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_
#define _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_

#include <string>
#include <vector>
#include <atomic>

#include "CellEngineTypes.h"

struct CellEngineState
{
public:
    std::uint64_t CellStateId;
    bool ExecuteStateBool;
    std::string CellStateFileName;
public:
    std::uint64_t ChosenStructureIndex;
//    std::uint64_t MaxDeepLevel;
//    bool ReverseColorOfPieces;
//public:
//    bool WriteTestStatisticsToFile;
//    bool PrintTestComparisonResult;
//public:
//    std::uint64_t NumberOfAllMoves;
//    std::uint64_t NumberOfAllWhitePiecesMoves;
//    std::uint64_t NumberOfAllBlackPiecesMoves;
//    std::uint64_t NumberOfAllWhiteMovesWhenCheckedKing;
//    std::uint64_t NumberOfAllBlackMovesWhenCheckedKing;
//    std::uint64_t NumberOfAllOppositePiecesRemovedByMoveOfWhite;
//    std::uint64_t NumberOfAllOppositePiecesRemovedByMoveOfBlack;
//    std::uint64_t NumberOfAllWhiteMovesCancelledBecauseOfDiscoveredPossibleCheckOfKing;
//    std::uint64_t NumberOfAllBlackMovesCancelledBecauseOfDiscoveredPossibleCheckOfKing;
//    std::uint64_t NumberOfAllFoundMovesPathsToMateWhite;
//    std::uint64_t NumberOfAllFoundMovesPathsToMateBlack;
//    std::uint64_t NumberOfAllFoundMovesPathsToStalemate;
//    PathScoreType TheHighestFoundCaptureMovesPathScoreForAdvantageOfWhite;
//    PathScoreType TheHighestFoundCaptureMovesPathScoreForAdvantageOfBlack;
//    PathScoreType TheHighestFoundMateMovesPathScoreForAdvantageOfWhite;
//    PathScoreType TheHighestFoundMateMovesPathScoreForAdvantageOfBlack;
};

class CellEngineConfigurationFileReaderWriter
{
public:
//    std::string ChessBoardFileName;
//
//    LengthType MaxDeepLevel;
//
//    bool ReverseColorOfPieces;
//
    bool MultiThreaded;
    bool SetProcessPriorityHighest;
//
//    bool ReadPreviousMovePathsFromFile;
//
//    bool CompareTestStatistics;
//    bool WriteTestsStatisticsToFile;
//    std::uint64_t MaxNumberOfChessBoardsForFastExecution;
//
//    std::uint8_t Color1;
//    std::uint8_t Color2;
//    std::uint8_t Color3;
//    std::uint8_t Color4;
//
//    std::string PieceString[NumberOfPieces];
//
//    PathScoreType PiecePower[NumberOfPieces];
//
    bool PrintLogToConsole = true;
    bool PrintLogToFiles = true;

    bool PrintLogLineNumberToConsole = false;
    bool PrintLogDateTimeToConsole = false;
    bool PrintLogProcessIdToConsole = false;
    bool PrintLogProcessPriorityLevelToConsole = false;
    bool PrintLogThreadIdToConsole = false;

    bool PrintLogLineNumberToFile = true;
    bool PrintLogDateTimeToFile = true;
    bool PrintLogProcessIdToFile = false;
    bool PrintLogProcessPriorityLevelToFile = false;
    bool PrintLogThreadIdToFile = false;
//
//    uint64_t MaximalNumberOfLinesInOneFile = 100000;
//
//    bool PrintMPIProcessesInformationData;
//    bool PrintActualStartPositionOfPieces;
//
//    bool CleanScreenBeforePrintingActualState;
//
//    bool PrintPossibleMoveFalseState;
//    bool StopWhenPossibleMoveFalseState;
//    bool PrintDetailedMoveDescription;
//    bool PrintGeneralInterMoveInfo;
//    bool StopWhenMove;
//    bool StopWhenMateFound;
//
//    bool PrintActualMove;
//
//    bool PrintActualChessBoardState;
//    bool PrintActualStateWhenMateFound;
//    bool PrintChessBoardStateForCheck;
//    bool PrintWhenRemovingPiece;
//    bool StopWhenRemovingPiece;
//
//    bool PrintTheHighestScoredMovesPathToCaptureWhenPrintingStatistics;
//    bool PrintTheHighestScoredMovesPathToCaptureOnChessBoardWhenPrintingStatistics;
//    bool PrintTheHighestScoredMovesPathToMateWhenPrintingStatistics;
//    bool PrintTheHighestScoredMovesPathToMateOnChessBoardWhenPrintingStatistics;
//    bool PrintEveryPathToMateWhenPrintingStatistics;
//    bool PrintEveryMoveToMateOnChessBoardWhenPrintingStatistics;
//    bool StopWhenPrintEveryMoveToMateOnChessBoardWhenPrintingStatistics;
//    bool PrintChosenMovesPathWhenPrintingStatistics;
//    bool PrintChosenMovesPathOnChessBoardWhenPrintingStatistics;
//
//    std::uint64_t PrintOnlyIfChosenMovesPathsToMateIsSmallerThanThisNumber;
//    bool PrintLookingForEqualMoveFromChosenMovesPaths;
//    bool PrintToMovesFile;
//
//    std::uint64_t ModuloToWriteMove;
//
//    bool WaitForKeyOnTheEnd;
public:
    std::vector<CellEngineState> CellEngineStates;
//public:
//    bool FirstMove;
//    std::string FirstPossibleMove;
public:
    std::string ConfigFileName;
public:
    CellEngineConfigurationFileReaderWriter() = default;
public:
    void ReadChessConfigurationFile(const char* ConfigFileNameParameter, const std::unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer);
public:
    void SaveTestStatisticsToFile(const CellEngineState& CellEngineStateObject) const;
};

inline CellEngineConfigurationFileReaderWriter CellEngineConfigurationFileReaderWriterObject;

#endif