
#include <memory>
#include <exception>

#include <string>

#include "ExceptionsMacro.h"
#include "CellEnginePDBWindowGL.h"

using namespace std;

#pragma region DrawPDB

PDBWindowGL::PDBWindowGL() : WindowGL(), PDBDataFileObjectPointer(nullptr), ShowBonds(false), ShowPDBSize(1.00f), RefreshListOfDrawing(false), ProjectionType(false), ChosenAtomIndex(-1)
{
    try
    {
        char* FileName;

        if (__argc > 1)
            FileName = __argv[1];
        else
        {
            MessageBox(nullptr, "Lack of file name in program parameters", "Cell Engine View PDB", MB_OK | MB_ICONWARNING);
            PostQuitMessage(0);
            return;
        }

        if (OpenPDBFile(FileName) == false)
            PostQuitMessage(0);
	}
    CATCH("initiation PDBWindowGL")
}

bool PDBWindowGL::OpenPDBFile(const char* FileName)
{
	try
	{
		PDBDataFileObjectPointer = make_unique<PDBDataFile>(FileName);
	}
	CATCH("Error with loading PDB object");

	return true;
}

void PDBWindowGL::DrawAtoms(const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors) const
{
    try
    {
        glPushMatrix();

        DoubleVectorType MassCenter = LengthUnit * PDBDataFileObjectPointer->MassCenter();
        glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

        GLUquadricObj* Quadriga = gluNewQuadric();
        gluQuadricDrawStyle(Quadriga, GLU_FILL);
        glShadeModel(GL_SMOOTH);

        const IntType HowManyPointsInEachDimension = 10;

        glInitNames();
        glPushName(-1);

        for (const Atom& AtomObjectPtr : PDBDataFileObjectPointer->GetAtoms())
        {
            DoubleVectorType AtomPosition = LengthUnit * AtomObjectPtr.Position();

            if (MakeColors == true)
                ChooseAtomColor(AtomObjectPtr.Name, 1.0f);

            glPushMatrix();
            glTranslated(AtomPosition.X, AtomPosition.Y, AtomPosition.Z);
            glLoadName(AtomObjectPtr.AtomIndex);
            gluSphere(Quadriga, LengthUnit * AtomSizeLengthUnit, HowManyPointsInEachDimension, HowManyPointsInEachDimension);
            glPopMatrix();
        }

        glPopName();

        gluDeleteQuadric(Quadriga);
        glPopMatrix();
	}
    CATCH("drawing atoms")
}

void PDBWindowGL::DrawChosenAtom(const double LengthUnit, const double AtomSizeLengthUnit, IntType LocalChosenAtomIndex, bool UseGrid) const
{
    try
    {
        if (LocalChosenAtomIndex < 0)
            return;

        glPushMatrix();

        GLUquadricObj* Quadric = gluNewQuadric();
        if (UseGrid)
            gluQuadricDrawStyle(Quadric, GLU_LINE);
        glLineWidth(1.0f);
        glShadeModel(GL_SMOOTH);
        DoubleVectorType ChosenAtomPosition = LengthUnit * (PDBDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Position() - PDBDataFileObjectPointer->MassCenter());
        glTranslated(ChosenAtomPosition.X, ChosenAtomPosition.Y, ChosenAtomPosition.Z);
        gluSphere(Quadric, LengthUnit * AtomSizeLengthUnit, 15, 15);
        gluDeleteQuadric(Quadric);

        glPopMatrix();
	}
    CATCH("drawing chosen atom")
}

void PrintTextInGL(const char* Text, IntType CharactersNumber, UnsignedIntType Font, IntType FirstCharCode)
{
    try
    {
        if (Text == nullptr)
            return;

        glPushAttrib(GL_LIST_BIT);
        glListBase(Font - FirstCharCode);
        glCallLists(CharactersNumber, GL_UNSIGNED_BYTE, Text);
        glPopAttrib();
    }
    CATCH("printing text in gl")
}

void PDBWindowGL::DrawChosenAtomDescription(IntType LocalChosenAtomIndex, IntType BitmapFont)
{
    try
    {
        ChosenAtomDescription = to_string(PDBDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Serial) + "." + PDBDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Name + "(" + PDBDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).ResName + ")";

        glPushMatrix();
        glLoadIdentity();

        float Coordinates = UserAreaHeight / (float)UserAreaWidth;
        float TextPosition = (!ProjectionType) ? -0.095f : -2.85f;

        glRasterPos3f(TextPosition, TextPosition * Coordinates, -0.3f);
        PrintTextInGL(ChosenAtomDescription.c_str(), 256, BitmapFont, 32);

        glPopMatrix();
	}
    CATCH("drawing chosen atom description")
}

void PDBWindowGL::ChooseAtomColor(const string_view Name, const float Alpha = 1.0f) const
{
    try
    {
//        if (Name == "CA")
//            glColor4f(0.25f, 0.75f, 0.75f, Alpha);
//        else
//        if (Name == "CB")
//            glColor4f(1.00f, 0.00f, 0.00f, Alpha);

        switch(Name[0])
        {
            case 'C': glColor4f(0.25f, 0.75f, 0.75f, Alpha); break;
            case 'O': glColor4f(1.00f, 0.00f, 0.00f, Alpha); break;
            case 'H': glColor4f(1.00f, 1.00f, 1.00f, Alpha); break;
            case 'N': glColor4f(0.00f, 0.00f, 1.00f, Alpha); break;
            case 'P': glColor4f(0.50f, 0.50f, 0.20f, Alpha); break;
            default: glColor4f(0.50f, 0.50f, 0.50f, Alpha); break;
        }
    }
    CATCH("chossing atom color")
}

void PDBWindowGL::DrawBonds(const double LengthUnit, bool MakeColors) const
{
    try
    {
        glPushMatrix();

        DoubleVectorType MassCenter = LengthUnit * PDBDataFileObjectPointer->MassCenter();;
        glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

        glColor3f(0, 0, 0);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

        glShadeModel(GL_SMOOTH);
        glLineWidth(2);
        glBegin(GL_LINES);

        for (IntType AtomIndex1 = 0; AtomIndex1 < PDBDataFileObjectPointer->GetNumberOfAtoms(); AtomIndex1++)
            for (IntType AtomIndex2 = AtomIndex1 + 1; AtomIndex2 < PDBDataFileObjectPointer->GetNumberOfAtoms(); AtomIndex2++)
            {
                DoubleVectorType Atom1Position = LengthUnit * PDBDataFileObjectPointer->GetAtom(AtomIndex1).Position();
                DoubleVectorType Atom2Position = LengthUnit * PDBDataFileObjectPointer->GetAtom(AtomIndex2).Position();
                DoubleVectorType Position12 = Atom2Position - Atom1Position;

                if (Position12.Length() < 0.17)
                {
                    if (MakeColors)
                        ChooseAtomColor(PDBDataFileObjectPointer->GetAtom(AtomIndex1).Name.c_str(), 1.0f);

                    glVertex3d(Atom1Position.X, Atom1Position.Y, Atom1Position.Z);

                    if (MakeColors)
                        ChooseAtomColor(PDBDataFileObjectPointer->GetAtom(AtomIndex2).Name.c_str(), 1.0f);

                    glVertex3d(Atom2Position.X, Atom2Position.Y, Atom2Position.Z);
                }
            }

        glEnd();

        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        glPopMatrix();
	}
    CATCH("drawing bonds")
}

UnsignedIntType PDBWindowGL::CreateListOfDrawing(const double LengthUnit)
{
    GLuint DrawList;

    try
    {
        DrawList = glGenLists(1);

        glNewList(DrawList, GL_COMPILE);
        glColor3f(1, 1, 1);

        if (ShowBonds)
            DrawBonds(LengthUnit, true);

        if (ShowPDBSize != 0)
            DrawAtoms(LengthUnit, ShowPDBSize, true);

        glEndList();
    }
    CATCH("creating list of drawing")

	return DrawList;
}

LRESULT PDBWindowGL::WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	UnsignedIntType Result;

	try
    {
        Result = WindowGL::WndProc(hWnd, Message, wParam, lParam);

        switch (Message)
        {
            case WM_DESTROY: KillTimer(HandleWindow, 1); break;

            case WM_TIMER:
            {
                if (PDBDataFileObjectPointer->ChosenStructureIndex < PDBDataFileObjectPointer->GetNumberOfStructures() - 1)
                {
                    PDBDataFileObjectPointer->ChosenStructureIndex++;
                    RefreshListOfDrawing = true;
                    DrawStage();
                }
                else
                    KillTimer(HandleWindow, 1);
            }

            case WM_KEYDOWN:
            {
                bool CallDrawStage = false;
                switch (wParam)
                {
                    case VK_F1: MessageBox(HandleWindow, "Shortcut Keys: \nZ - show / hide bonds \nX - size of atoms \nC - projection mode", "PDB VIEWER", MB_OK); break;
                    case VK_F2: MessageBox(HandleWindow, ChosenAtomDescription.c_str(), "PDB VIEWER", MB_OK); break;
                    case 'Z': ShowBonds = !ShowBonds; CallDrawStage = true; break;
                    case 'X':
                    {
                        float NewSize = -1.0f;

                        if (ShowPDBSize == 0.0f)
                            NewSize = 0.25f;
                        if (ShowPDBSize == 0.25f)
                            NewSize = 0.5f;
                        if (ShowPDBSize == 0.5f)
                            NewSize = 1.0f;
                        if (ShowPDBSize == 1.0f)
                            NewSize = 0.0f;
                        if (NewSize == -1.0f)
                            NewSize = 1.0f;

                        ShowPDBSize = NewSize;
                    }
                    CallDrawStage = true;
                    break;

                    case 'C':
                    {
                        ProjectionType = !ProjectionType;
                        SetStage(ProjectionType);
                        CallDrawStage = true;
                    }
                    break;

                    case 'N':
                    {
                        if (PDBDataFileObjectPointer->ChosenStructureIndex < PDBDataFileObjectPointer->GetNumberOfStructures() - 1)
                            PDBDataFileObjectPointer->ChosenStructureIndex++;
                        CallDrawStage = true;
                    }
                    break;
                    case 'M':
                    {
                        if (PDBDataFileObjectPointer->ChosenStructureIndex > 0)
                            PDBDataFileObjectPointer->ChosenStructureIndex--;
                        CallDrawStage = true;
                    }
                    break;
                    case 'F':
                    {
                        PDBDataFileObjectPointer->ChosenStructureIndex = 0;

                        if (SetTimer(hWnd, 1, 50, NULL) == 0 )
                            MessageBox(hWnd, "Nie udało się ustawić timera", "", MB_OK | MB_ICONERROR);
                    }
                    break;

                    default: break;
                }

                if (CallDrawStage)
                {
                    RefreshListOfDrawing = true;
                    DrawStage();
                }

                break;
            }

            case WM_LBUTTONDOWN:

                if (ShowPDBSize == 0)
                    break;

                POINT CursorMouseIndex;
                CursorMouseIndex.x = LOWORD(lParam);
                CursorMouseIndex.y = HIWORD(lParam);

                IntType PrevIndex = ChosenAtomIndex;
                IntType NextIndex = ChooseAtom(CursorMouseIndex);

                if (NextIndex != -1)
                    ChosenAtomIndex = NextIndex;

                if (PrevIndex != -1 || ChosenAtomIndex != -1)
                    DrawStage();

                break;
        }
    }
    CATCH("execution of PDBWindowGL::WndProc")

	return Result;
}

IntType PDBWindowGL::ChooseAtom(POINT MouseCursorPosition)
{
    try
    {
        const IntType MarkingBufferSize = 1024;
        unsigned int MarkingBuffer[MarkingBufferSize];
        ZeroMemory(MarkingBuffer, MarkingBufferSize);
        glSelectBuffer(MarkingBufferSize, MarkingBuffer);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();

        int Viewport[4];
        glGetIntegerv(GL_VIEWPORT, Viewport);
        gluPickMatrix(MouseCursorPosition.x, UserAreaHeight - MouseCursorPosition.y, 1, 1, Viewport);
        float wsp = static_cast<float>(UserAreaHeight) / static_cast<float>(UserAreaWidth);
        if (!ProjectionType)
            glFrustum(-0.1, 0.1, wsp * -0.1, wsp * 0.1, 0.3, 100.0);
        else
            glOrtho(-3, 3, wsp * -3, wsp * 3, 0.3, 100.0);

        glMatrixMode(GL_MODELVIEW);

        glRenderMode(GL_SELECT);
        DrawStage();

        IntType HitsNumber = glRenderMode(GL_RENDER);

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();

        glMatrixMode(GL_MODELVIEW);

        if (HitsNumber > 0)
        {
            UnsignedIntType ClosestAtomIndex = MarkingBuffer[3];
            UnsignedIntType ClosestAtomDistance = MarkingBuffer[1];
            IntType CurrentIndex = 0;
            for (IntType HitIndex = 0; HitIndex < HitsNumber; HitIndex++)
            {
                if (MarkingBuffer[CurrentIndex + 1] < ClosestAtomDistance)
                {
                    ClosestAtomDistance = MarkingBuffer[CurrentIndex + 1];
                    if (MarkingBuffer[CurrentIndex] > 0)
                        ClosestAtomIndex = MarkingBuffer[CurrentIndex + 3];
                }
                CurrentIndex += 3 + MarkingBuffer[CurrentIndex];
            }
            return ClosestAtomIndex;
        }
        else
            return -1;
	}
    CATCH("choosing atom")

    return -1;
}

#pragma endregion

UnsignedIntType CreateFontGeneral(HWND HandleWindow, const char* FontName, IntType HeightInPixels, bool Bold, bool Italics, IntType FirstCharCode, IntType LastCharCode)
{
    UnsignedIntType FirstListIndex;

    try
    {
        FirstListIndex = glGenLists(LastCharCode + 1 - FirstCharCode);

        HFONT FontHandle = CreateFont(HeightInPixels, 0, 0, 0, Bold ? FW_BOLD : FALSE, Italics ? TRUE : FALSE, FALSE, FALSE, ANSI_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY, FF_DONTCARE | DEFAULT_PITCH, FontName);

        HDC HandleDC = GetDC(HandleWindow);
        auto DCFontHandle = static_cast<HFONT>(SelectObject(HandleDC, FontHandle));

        wglUseFontBitmaps(HandleDC, FirstCharCode, LastCharCode + 1 - FirstCharCode, FirstListIndex);

        SelectObject(HandleDC, DCFontHandle);
        DeleteObject(FontHandle);
    }
    CATCH("creating font general")

	return FirstListIndex;
}

void PDBWindowGL::DrawActors()
{
    try
    {
        ShowRenderingFrequency();

        const double LengthUnit = 0.1;

        static GLuint DrawingList = 0;
        if (glIsList(DrawingList) == false || RefreshListOfDrawing == true)
        {
            if (RefreshListOfDrawing == true)
                glDeleteLists(DrawingList, 1);
            DrawingList = CreateListOfDrawing(LengthUnit);
            RefreshListOfDrawing = false;
        }

        static UnsignedIntType BitmapFont = 0;
        if (BitmapFont == 0)
            BitmapFont = CreateFontGeneral(HandleWindow, "Calibri", 20, true, false, 32, 255);

        glCallList(DrawingList);
        int RenderingMode = 0;
        glGetIntegerv(GL_RENDER_MODE, &RenderingMode);
        if (RenderingMode == GL_RENDER)
        {
            if (ChosenAtomIndex >= 0 && ShowPDBSize > 0)
            {
                glColor3f(1, 1, 0);
                DrawChosenAtom(LengthUnit, 1.05 * ShowPDBSize, ChosenAtomIndex, false);
                glColor3f(5, 0, 5);
                DrawChosenAtomDescription(ChosenAtomIndex, BitmapFont);
            }
        }
	}
    CATCH("drawing actors")
}

#pragma region SourceOfLights

void PDBWindowGL::LightSources()
{
	BackgroundLightIntensity = 0.5f;
	MilkyBulb(0.5f);
}

void PDBWindowGL::MilkyBulb(float jasnosc)
{
    try
    {
        const float Color[4] = { jasnosc, jasnosc, jasnosc, 1.0f };
        const float Position[4] = { 5.0f, 0.0f, 5.0f, 1.0f };
        glLightfv(GL_LIGHT1, GL_POSITION, Position);
        glLightfv(GL_LIGHT1, GL_DIFFUSE, Color);
        glEnable(GL_LIGHT1);
	}
    CATCH("setting milkybulb")
}

#pragma endregion