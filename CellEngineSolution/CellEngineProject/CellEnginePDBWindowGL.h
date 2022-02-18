#pragma once

#ifndef CELL_ENGINE_PDB_WINDOW_GL_H
#define CELL_ENGINE_PDB_WINDOW_GL_H

#include <memory>

#include "VectorType.h"

#include "CellEngineWindowGL.h"
#include "CellEnginePDBDataFile.h"

class PDBWindowGL : public WindowGL
{
private:
	void DrawActors() final;
private:
	void LightSources() final;
	static void MilkyBulb(float Brightness);
public:
	PDBWindowGL();
	~PDBWindowGL() = default;
public:
    LRESULT WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam) final;
private:
    void StartTimerEvent();
    void FilmTimerEvent();
private:
    void ChangeElementsSize();
    void ChangeShowOfBonds();
    void ShowNextStructure();
    void ShowPrevStructure();
private:
    void KeyboardKeyDownPressedEvent(WPARAM wParam);
    void MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam);
private:
    void FullDrawStage();
private:
    std::string ChosenElementDescription;
	std::unique_ptr<PDBDataFile> PDBDataFileObjectPointer;
	bool OpenPDBFile(const char* FileName);
	void DrawElements(const double LengthUnit, const double ElementSizeLengthUnit, bool MakeColors) const;
	void DrawBonds(const double LengthUnit, bool MakeColors) const;
	void ChooseElementColor(const std::string_view Name, const float Alpha) const;
    [[nodiscard]] UnsignedIntType CreateListOfDrawing(const double LengthUnit);
private:
	bool ShowBonds;
	float ShowPDBSize;
	bool RefreshListOfDrawing;
	IntType ChooseElement(POINT MouseCursorPosition);
	void DrawChosenElement(const double LengthUnit, const double ElementSizeLengthUnit, IntType LocalChosenElementIndex, bool UseGrid) const;
	void DrawChosenElementDescription(IntType LocalChosenElementIndex, IntType BitmapFont);
	IntType ChosenElementIndex;
};

#endif