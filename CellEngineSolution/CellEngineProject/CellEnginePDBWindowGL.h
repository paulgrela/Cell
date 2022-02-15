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
    std::string ChosenAtomDescription;
	std::unique_ptr<PDBDataFile> PDBDataFileObjectPointer;
	bool OpenPDBFile(const char* FileName);
	void DrawAtoms(const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors) const;
	void DrawBonds(const double LengthUnit, bool MakeColors) const;
	void ChooseAtomColor(const std::string_view Name, const float Alpha) const;
    [[nodiscard]] UnsignedIntType CreateListOfDrawing(const double LengthUnit);
private:
	bool ShowBonds;
	float ShowPDBSize;
	bool RefreshListOfDrawing;
	bool ProjectionType;
	IntType ChooseAtom(POINT MouseCursorPosition);
	void DrawChosenAtom(const double LengthUnit, const double AtomSizeLengthUnit, IntType LocalChosenAtomIndex, bool UseGrid) const;
	void DrawChosenAtomDescription(IntType LocalChosenAtomIndex, IntType BitmapFont);
	IntType ChosenAtomIndex;
};

#endif