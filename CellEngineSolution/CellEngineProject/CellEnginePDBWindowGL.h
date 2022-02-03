#pragma once

#ifndef CELL_ENGINE_PDB_WINDOW_GL_H
#define CELL_ENGINE_PDB_WINDOW_GL_H

#include "VectorType.h"

#include "CellEngineWindowGL.h"
#include "CellEnginePDBDataFile.h"


class PDBWindowGL : public WindowGL
{
private:
	void DrawActors();
private:
	void LightSources();
	void MilkyBulb(float Brightness);
	void YellowGreenMilkyBulbs();
	void Reflector(float FlareBrightness = 1.0f, float BrightnessDisperesed = 0.3f);
public:
	PDBWindowGL();
	~PDBWindowGL();
private:
	PDBDataFile* PDBDataFileObjectPointer;
	bool OpenPDBFile(const char* FileName);
	void DrawAtoms(PDBDataFile* PDBDataFileObject, const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors) const;
	void DrawBonds(PDBDataFile* PDBDataFileObject, const double LengthUnit, bool MakeColors) const;
	void ChooseAtomColor(const char* AtomSymbol, const float Alpha) const;
	UnsignedIntType CreateListOfDrawing(const double LengthUnit) const;
	LRESULT WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
private:
	bool ShowBonds; //wlaczanie/wylaczanie wiezow i sfer
	float ShowPDBSize;
	bool RefreshListOfDrawing;
	bool ProjectionType; //false=perspektywa,true=izometryczne
	IntType ChooseAtom(POINT MouseCursorPosition);
	void DrawChosenAtom(PDBDataFile* PDBDataFileObject, const double LengthUnit, const double AtomSizeLengthUnit, IntType ChosenAtomIndex, bool UseGrid) const;
	void DrawChosenAtomDescription(IntType ChosenAtomIndex, IntType BitmapFont);
	IntType ChosenAtomIndex;
} 
PDBWindowGLObject;

WindowGL* WindowGLPointer = &PDBWindowGLObject;

#endif