/*
#include <memory>
#include <exception>

#include <string>

#include "ExceptionsMacro.h"
#include "DateTimeUtils.h"
#include "CellEnginePDBWindowGL.h"

using namespace std;

#pragma region DrawPDB

PDBWindowGL::PDBWindowGL(const string_view FileName) : WindowGL(), PDBDataFileObjectPointer(nullptr), ShowBonds(false), ShowPDBSize(1.00f), RefreshListOfDrawing(false), ChosenElementIndex(-1)
{
    try
    {
        CameraShift = 3;

        if (OpenPDBFile(FileName) == false)
            PostQuitMessage(0);
	}
    CATCH("initiation PDBWindowGL")
}

bool PDBWindowGL::OpenPDBFile(const string_view FileName)
{
	try
	{
        PDBDataFileObjectPointer = make_unique<CellEnginePDBDataFile>(FileName);
	}
	CATCH("Error with loading PDB object");

	return true;
}

void PDBWindowGL::DrawElements(const double LengthUnit, const double ElementSizeLengthUnit, bool MakeColors) const
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

        for (const Element& ElementObject : PDBDataFileObjectPointer->GetElements())
        {
            DoubleVectorType ElementPosition = LengthUnit * ElementObject.Position();

            if (MakeColors == true)
                ChooseElementColor(ElementObject.Name, 1.0f);

            glPushMatrix();
            glTranslated(ElementPosition.X, ElementPosition.Y, ElementPosition.Z);
            glLoadName(ElementObject.ElementIndex);
            gluSphere(Quadriga, LengthUnit * ElementSizeLengthUnit, HowManyPointsInEachDimension, HowManyPointsInEachDimension);
            glPopMatrix();
        }

        glPopName();

        gluDeleteQuadric(Quadriga);
        glPopMatrix();
    }
    CATCH("drawing elements")
}

void PDBWindowGL::DrawChosenElement(const double LengthUnit, const double ElementSizeLengthUnit, IntType LocalChosenElementIndex, bool UseGrid) const
{
    try
    {
        if (LocalChosenElementIndex < 0)
            return;

        glPushMatrix();

        GLUquadricObj* Quadric = gluNewQuadric();
        if (UseGrid)
            gluQuadricDrawStyle(Quadric, GLU_LINE);
        glLineWidth(1.0f);
        glShadeModel(GL_SMOOTH);
        DoubleVectorType ChosenElementPosition = LengthUnit * (PDBDataFileObjectPointer->GetElement(LocalChosenElementIndex).Position() - PDBDataFileObjectPointer->MassCenter());
        glTranslated(ChosenElementPosition.X, ChosenElementPosition.Y, ChosenElementPosition.Z);
        gluSphere(Quadric, LengthUnit * ElementSizeLengthUnit, 15, 15);
        gluDeleteQuadric(Quadric);

        glPopMatrix();
	}
    CATCH("drawing chosen element")
}

void PDBWindowGL::DrawChosenElementDescription(IntType LocalChosenElementIndex, IntType BitmapFont)
{
    try
    {
        ChosenElementDescription = to_string(PDBDataFileObjectPointer->GetElement(LocalChosenElementIndex).Serial) + "." + PDBDataFileObjectPointer->GetElement(LocalChosenElementIndex).Name + "(" + PDBDataFileObjectPointer->GetElement(LocalChosenElementIndex).ResName + ")";

        glPushMatrix();

        glLoadIdentity();

        float Coordinates = UserAreaHeight / (float)UserAreaWidth;
        float TextPosition = (IsometricProjection == false) ? -0.095f : -2.85f;

        glRasterPos3f(TextPosition, TextPosition * Coordinates, -0.3f);

        glPushAttrib(GL_LIST_BIT);
        glListBase(BitmapFont - 32);
        glCallLists(ChosenElementDescription.length(), GL_UNSIGNED_BYTE, ChosenElementDescription.c_str());
        glPopAttrib();

        glPopMatrix();
	}
    CATCH("drawing chosen Element description")
}

void PDBWindowGL::ChooseElementColor(const string_view Name, const float Alpha = 1.0f) const
{
    try
    {
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
    CATCH("chossing element color")
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

        for (IntType ElementIndex1 = 0; ElementIndex1 < PDBDataFileObjectPointer->GetNumberOfElements(); ElementIndex1++)
            for (IntType ElementIndex2 = ElementIndex1 + 1; ElementIndex2 < PDBDataFileObjectPointer->GetNumberOfElements(); ElementIndex2++)
            {
                DoubleVectorType Element1Position = LengthUnit * PDBDataFileObjectPointer->GetElement(ElementIndex1).Position();
                DoubleVectorType Element2Position = LengthUnit * PDBDataFileObjectPointer->GetElement(ElementIndex2).Position();
                DoubleVectorType Position12 = Element2Position - Element1Position;

                if (Position12.Length() < 0.17)
                {
                    if (MakeColors)
                        ChooseElementColor(PDBDataFileObjectPointer->GetElement(ElementIndex1).Name.c_str(), 1.0f);

                    glVertex3d(Element1Position.X, Element1Position.Y, Element1Position.Z);

                    if (MakeColors)
                        ChooseElementColor(PDBDataFileObjectPointer->GetElement(ElementIndex2).Name.c_str(), 1.0f);

                    glVertex3d(Element2Position.X, Element2Position.Y, Element2Position.Z);
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

        if (ShowBonds == true)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            DrawBonds(LengthUnit, true);

            const auto stop_time = chrono::high_resolution_clock::now();

            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of drawing bonds has taken time: ","executing printing duration_time")));
        }

        if (ShowPDBSize != 0)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            DrawElements(LengthUnit, ShowPDBSize, true);

            const auto stop_time = chrono::high_resolution_clock::now();

            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of drawing elements has taken time: ","executing printing duration_time")));
        }

        glEndList();
    }
    CATCH("creating list of drawing")

	return DrawList;
}

void PDBWindowGL::FullDrawStage()
{
    try
    {
        RefreshListOfDrawing = true;
        DrawStage();
    }
    CATCH("full draw stage")
}

void PDBWindowGL::StartTimerEvent()
{
    try
    {
        PDBDataFileObjectPointer->ChosenStructureIndex = 0;

        if (SetTimer(HandleWindow, 1, 50, nullptr) == 0 )
            MessageBox(HandleWindow, "Nie udało się ustawić timera", "", MB_OK | MB_ICONERROR);
    }
    CATCH("start timer event")
}

void PDBWindowGL::FilmTimerEvent()
{
    try
    {
        if (PDBDataFileObjectPointer->ChosenStructureIndex < PDBDataFileObjectPointer->GetNumberOfStructures() - 1)
        {
            PDBDataFileObjectPointer->ChosenStructureIndex++;
            FullDrawStage();
        }
        else
            KillTimer(HandleWindow, 1);
    }
    CATCH("timer event")
}

void PDBWindowGL::ChangeElementsSize()
{
    try
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

        FullDrawStage();
    }
    CATCH("changing elements size")
}

void PDBWindowGL::ShowNextStructure()
{
    try
    {
        if (PDBDataFileObjectPointer->ChosenStructureIndex < PDBDataFileObjectPointer->GetNumberOfStructures() - 1)
            PDBDataFileObjectPointer->ChosenStructureIndex++;
        FullDrawStage();
    }
    CATCH("showing next structure")
}

void PDBWindowGL::ShowPrevStructure()
{
    try
    {
        if (PDBDataFileObjectPointer->ChosenStructureIndex > 0)
            PDBDataFileObjectPointer->ChosenStructureIndex--;
        FullDrawStage();
    }
    CATCH("showing previous structure")
}

void PDBWindowGL::ChangeShowOfBonds()
{
    try
    {
        ShowBonds = !ShowBonds;
        FullDrawStage();
    }
    CATCH("changing showing of bonds")
}

void PDBWindowGL::KeyboardKeyDownPressedEvent(WPARAM wParam)
{
    try
    {
        switch (wParam)
        {
            case VK_F1: MessageBox(HandleWindow, "Shortcut Keys: \nZ - show / hide bonds \nX - size of Elements \nC - projection mode", "PDB VIEWER", MB_OK); break;
            case VK_F2: MessageBox(HandleWindow, ChosenElementDescription.c_str(), "PDB VIEWER", MB_OK); break;
            case 'Z': ChangeShowOfBonds(); break;
            case 'X': ChangeElementsSize(); break;
            case 'N': ShowNextStructure(); break;
            case 'M': ShowPrevStructure(); break;
            case 'F': StartTimerEvent(); break;

            default: break;
        }

        DrawStage();
    }
    CATCH("PDBWindow key down pressed event")
}

void PDBWindowGL::MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        if (ShowPDBSize == 0)
            return;

        POINT CursorMouseIndex;
        CursorMouseIndex.x = LOWORD(lParam);
        CursorMouseIndex.y = HIWORD(lParam);

        IntType PrevIndex = ChosenElementIndex;
        IntType NextIndex = ChooseElement(CursorMouseIndex);

        if (NextIndex != -1)
            ChosenElementIndex = NextIndex;

        if (PrevIndex != -1 || ChosenElementIndex != -1)
            DrawStage();
    }
    CATCH("mouse move event")
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
            case WM_TIMER: FilmTimerEvent(); break;
            case WM_KEYDOWN: KeyboardKeyDownPressedEvent(wParam); break;
            case WM_LBUTTONDOWN: MouseLeftButtonDownEvent(wParam, lParam); break;

            default: break;
        }
    }
    CATCH("execution of PDBWindowGL::WndProc")

	return Result;
}

IntType PDBWindowGL::ChooseElement(POINT MouseCursorPosition)
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
        float Factor = static_cast<float>(UserAreaHeight) / static_cast<float>(UserAreaWidth);
        if (IsometricProjection == false)
            glFrustum(-0.1, 0.1, Factor * -0.1, Factor * 0.1, 0.3, 100.0);
        else
            glOrtho(-3, 3, Factor * -3, Factor * 3, 0.3, 100.0);

        glMatrixMode(GL_MODELVIEW);

        glRenderMode(GL_SELECT);
        DrawStage();

        IntType HitsNumber = glRenderMode(GL_RENDER);

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();

        glMatrixMode(GL_MODELVIEW);

        if (HitsNumber > 0)
        {
            UnsignedIntType ClosestElementIndex = MarkingBuffer[3];
            UnsignedIntType ClosestElementDistance = MarkingBuffer[1];
            IntType CurrentIndex = 0;
            for (IntType HitIndex = 0; HitIndex < HitsNumber; HitIndex++)
            {
                if (MarkingBuffer[CurrentIndex + 1] < ClosestElementDistance)
                {
                    ClosestElementDistance = MarkingBuffer[CurrentIndex + 1];
                    if (MarkingBuffer[CurrentIndex] > 0)
                        ClosestElementIndex = MarkingBuffer[CurrentIndex + 3];
                }
                CurrentIndex += 3 + MarkingBuffer[CurrentIndex];
            }
            return ClosestElementIndex;
        }
        else
            return -1;
	}
    CATCH("choosing element")

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

        const auto start_time1 = chrono::high_resolution_clock::now();

        const double LengthUnit = 0.1;

        static GLuint DrawingList = 0;
        if (glIsList(DrawingList) == false || RefreshListOfDrawing == true)
        {
            if (RefreshListOfDrawing == true)
                glDeleteLists(DrawingList, 1);
            DrawingList = CreateListOfDrawing(LengthUnit);
            RefreshListOfDrawing = false;
        }

        const auto stop_time1 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of X1 has taken time: ","executing printing duration_time")));

        static UnsignedIntType BitmapFont = 0;
        if (BitmapFont == 0)
            BitmapFont = CreateFontGeneral(HandleWindow, "Calibri", 20, true, false, 32, 255);

        const auto start_time2 = chrono::high_resolution_clock::now();

        glCallList(DrawingList);

        const auto stop_time2 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time2, stop_time2, "Execution of X2 has taken time: ","executing printing duration_time")));


        const auto start_time3 = chrono::high_resolution_clock::now();

        int RenderingMode = 0;
        glGetIntegerv(GL_RENDER_MODE, &RenderingMode);
        if (RenderingMode == GL_RENDER)
            if (ChosenElementIndex >= 0 && ShowPDBSize > 0)
            {
                glColor3f(1, 1, 0);
                DrawChosenElement(LengthUnit, 1.05 * ShowPDBSize, ChosenElementIndex, false);
                glColor3f(5, 0, 5);
                DrawChosenElementDescription(ChosenElementIndex, BitmapFont);
            }

        const auto stop_time3 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time3, stop_time3, "Execution of X3 has taken time: " ,"executing printing duration_time")));
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
 */