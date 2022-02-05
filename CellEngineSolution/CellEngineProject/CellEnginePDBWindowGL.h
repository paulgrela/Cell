#pragma once

#ifndef CELL_ENGINE_PDB_WINDOW_GL_H
#define CELL_ENGINE_PDB_WINDOW_GL_H

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
	~PDBWindowGL();
private:
	PDBDataFile* PDBDataFileObjectPointer;
	bool OpenPDBFile(const char* FileName);
	static void DrawAtoms(PDBDataFile* PDBDataFileObject, const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors);
	static void DrawBonds(PDBDataFile* PDBDataFileObject, const double LengthUnit, bool MakeColors);
	static void ChooseAtomColor(const char* AtomSymbol, const float Alpha);
    [[nodiscard]] UnsignedIntType CreateListOfDrawing(const double LengthUnit) const;
	LRESULT WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) final;
private:
	bool ShowBonds;
	float ShowPDBSize;
	bool RefreshListOfDrawing;
	bool ProjectionType;
	IntType ChooseAtom(POINT MouseCursorPosition);
	static void DrawChosenAtom(PDBDataFile* PDBDataFileObject, const double LengthUnit, const double AtomSizeLengthUnit, IntType LocalChosenAtomIndex, bool UseGrid);
	void DrawChosenAtomDescription(IntType LocalChosenAtomIndex, IntType BitmapFont);
	IntType ChosenAtomIndex;
} 
PDBWindowGLObject;

WindowGL* WindowGLPointer = &PDBWindowGLObject;

#endif