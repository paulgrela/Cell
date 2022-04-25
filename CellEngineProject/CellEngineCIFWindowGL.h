#pragma once

#ifndef CELL_ENGINE_CIF_WINDOW_GL_H
#define CELL_ENGINE_CIF_WINDOW_GL_H

/*
#include <memory>

#include "VectorType.h"

#include "CellEngineWindowGL.h"
#include "CellEngineCIFDataFile.h"

class CIFWindowGL : public WindowGL
{
private:
    void DrawActors() final;
private:
    void LightSources() final;
    static void MilkyBulb(float Brightness);
public:
    explicit CIFWindowGL(const std::string_view FileName);
    ~CIFWindowGL() = default;
public:
    LRESULT WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam) final;
private:
    void StartTimerEvent();
    void FilmTimerEvent();
private:
    void ChangeAtomsSize();
    void ChangeCameraShift();
    void ChangeShowOfBonds();
    void ChangeShowOfSpheres();
    void ShowNextStructure();
    void ShowPrevStructure();
private:
    void KeyboardKeyDownPressedEvent(WPARAM wParam);
    void MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam);
private:
    void FullDrawStage();
private:
    std::string ChosenAtomDescription;
    std::unique_ptr<CIFDataFile> CIFDataFileObjectPointer;
    bool OpenCIFFile(const std::string_view FileName);
    void DrawAtoms(const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors) const;
    void DrawBonds(const double LengthUnit, bool MakeColors) const;
    void ChooseAtomColor(const std::string_view Name, const float Alpha) const;
    [[nodiscard]] UnsignedIntType CreateListOfDrawing(const double LengthUnit);
private:
    bool ShowSpheres;
    bool ShowBonds;
    float ShowCIFSize;
    bool RefreshListOfDrawing;
    IntType ChooseAtom(POINT MouseCursorPosition);
    void DrawChosenAtom(const double LengthUnit, const double AtomSizeLengthUnit, IntType LocalChosenAtomIndex, bool UseGrid) const;
    void DrawChosenAtomDescription(IntType LocalChosenAtomIndex, IntType BitmapFont);
    IntType ChosenAtomIndex;
};
*/

#endif