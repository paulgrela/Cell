
#include <memory>
#include <exception>

#include <string>

#include "ExceptionsMacro.h"
#include "DateTimeUtils.h"
#include "CellEngineCIFWindowGL.h"

using namespace std;

#pragma region DrawCIF

CIFWindowGL::CIFWindowGL(const string_view FileName) : WindowGL(), CIFDataFileObjectPointer(nullptr), ShowBonds(false), ShowSpheres(false), ShowCIFSize(1.00f), RefreshListOfDrawing(false), ChosenAtomIndex(-1)
{
    try
    {
        if (OpenCIFFile(FileName) == false)
            PostQuitMessage(0);
    }
    CATCH("initiation CIFWindowGL")
}

bool CIFWindowGL::OpenCIFFile(const string_view FileName)
{
    try
    {
        CIFDataFileObjectPointer = make_unique<CIFDataFile>(FileName);
    }
    CATCH("Error with loading CIF object");

    return true;
}

void CIFWindowGL::DrawAtoms(const double LengthUnit, const double AtomsizeLengthUnit, bool MakeColors) const
{
    try
    {
        glPushMatrix();

        DoubleVectorType MassCenter = LengthUnit * CIFDataFileObjectPointer->MassCenter();
        glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

        GLUquadricObj* Quadriga = gluNewQuadric();

        if (ShowSpheres == true)
            gluQuadricDrawStyle(Quadriga, GLU_FILL);
        else
            gluQuadricDrawStyle(Quadriga, GLU_POINT);

        glShadeModel(GL_SMOOTH);

        const IntType HowManyPointsInEachDimension = 10;
        //const IntType HowManyPointsInEachDimension = 2;

        glInitNames();
        glPushName(-1);


        //TUTAJ DAĆ NIE WSZYSTKIE ATOMY ALE TE KTÓRE RYSUJĘ?
        //int i = 0;
        for (const Atom& AtomObjectPtr : CIFDataFileObjectPointer->GetAtoms())
        //if(i < 10000)
        {
            //i++;
            DoubleVectorType AtomPosition = LengthUnit * AtomObjectPtr.Position();

            if (MakeColors == true)
                ChooseAtomColor(AtomObjectPtr.Chain, 1.0f);

            glPushMatrix();
            glTranslated(AtomPosition.X, AtomPosition.Y, AtomPosition.Z);
            glLoadName(AtomObjectPtr.AtomIndex);
            gluSphere(Quadriga, LengthUnit * AtomsizeLengthUnit, HowManyPointsInEachDimension, HowManyPointsInEachDimension);
            glPopMatrix();
        }

        glPopName();

        gluDeleteQuadric(Quadriga);
        glPopMatrix();
    }
    CATCH("drawing Atoms")
}

void CIFWindowGL::DrawChosenAtom(const double LengthUnit, const double AtomsizeLengthUnit, IntType LocalChosenAtomIndex, bool UseGrid) const
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
        DoubleVectorType ChosenAtomPosition = LengthUnit * (CIFDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Position() - CIFDataFileObjectPointer->MassCenter());
        glTranslated(ChosenAtomPosition.X, ChosenAtomPosition.Y, ChosenAtomPosition.Z);
        gluSphere(Quadric, LengthUnit * AtomsizeLengthUnit, 15, 15);
        gluDeleteQuadric(Quadric);

        glPopMatrix();
    }
    CATCH("drawing chosen Atom")
}

void CIFWindowGL::DrawChosenAtomDescription(IntType LocalChosenAtomIndex, IntType BitmapFont)
{
    try
    {
        ChosenAtomDescription = to_string(CIFDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Serial) + "." + CIFDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).Name + "(" + CIFDataFileObjectPointer->GetAtom(LocalChosenAtomIndex).ResName + ")";
        ChosenAtomDescription = "A";

        glPushMatrix();

        glLoadIdentity();

        float Coordinates = UserAreaHeight / (float)UserAreaWidth;
        float TextPosition = (IsometricProjection == false) ? -0.095f : -2.85f;

        glRasterPos3f(TextPosition, TextPosition * Coordinates, -0.3f);

        glPushAttrib(GL_LIST_BIT);
        glListBase(BitmapFont - 32);
        glCallLists(ChosenAtomDescription.length(), GL_UNSIGNED_BYTE, ChosenAtomDescription.c_str());
        glPopAttrib();

        glPopMatrix();
    }
    CATCH("drawing chosen Atom description")
}

void CIFWindowGL::ChooseAtomColor(const string_view Chain, const float Alpha = 1.0f) const
{
    try
    {
        /*
        switch(Name[0])
        {
            case 'C': glColor4f(0.25f, 0.75f, 0.75f, Alpha); break;
            case 'O': glColor4f(1.00f, 0.00f, 0.00f, Alpha); break;
            case 'H': glColor4f(1.00f, 1.00f, 1.00f, Alpha); break;
            case 'N': glColor4f(0.00f, 0.00f, 1.00f, Alpha); break;
            case 'P': glColor4f(0.50f, 0.50f, 0.20f, Alpha); break;
            default: glColor4f(0.50f, 0.50f, 0.50f, Alpha); break;
        }
        */

        if(Chain.substr(0, 3) == "BAF")
            glColor4f(0.25f, 0.75f, 0.75f, Alpha);
        else
        if(Chain.substr(0, 3) == "BAE")
            glColor4f(1.00f, 0.00f, 0.00f, Alpha);
        else
        if(Chain.substr(0, 3) == "ATP")
            glColor4f(1.00f, 1.00f, 1.00f, Alpha);
        else
        if(Chain.substr(0, 3) == "BAR")
            glColor4f(0.00f, 0.00f, 1.00f, Alpha);
        else
        if(Chain.substr(0, 1) == "A")
            glColor4f(0.50f, 0.50f, 0.20f, Alpha);
        else
            glColor4f(0.50f, 0.50f, 0.50f, Alpha);
    }
    CATCH("chossing Atom color")
}

void CIFWindowGL::DrawBonds(const double LengthUnit, bool MakeColors) const
{
    try
    {
        glPushMatrix();

        DoubleVectorType MassCenter = LengthUnit * CIFDataFileObjectPointer->MassCenter();;
        glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

        glColor3f(0, 0, 0);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

        glShadeModel(GL_SMOOTH);
        glLineWidth(2);

        glBegin(GL_LINES);

        for (IntType AtomIndex1 = 0; AtomIndex1 < CIFDataFileObjectPointer->GetNumberOfAtoms(); AtomIndex1++)
            for (IntType AtomIndex2 = AtomIndex1 + 1; AtomIndex2 < CIFDataFileObjectPointer->GetNumberOfAtoms(); AtomIndex2++)
            {
                DoubleVectorType Atom1Position = LengthUnit * CIFDataFileObjectPointer->GetAtom(AtomIndex1).Position();
                DoubleVectorType Atom2Position = LengthUnit * CIFDataFileObjectPointer->GetAtom(AtomIndex2).Position();
                DoubleVectorType Position12 = Atom2Position - Atom1Position;

                if (Position12.Length() < 0.17)
                {
                    if (MakeColors)
                        ChooseAtomColor(CIFDataFileObjectPointer->GetAtom(AtomIndex1).Name, 1.0f);

                    glVertex3d(Atom1Position.X, Atom1Position.Y, Atom1Position.Z);

                    if (MakeColors)
                        ChooseAtomColor(CIFDataFileObjectPointer->GetAtom(AtomIndex2).Name, 1.0f);

                    glVertex3d(Atom2Position.X, Atom2Position.Y, Atom2Position.Z);
                }
            }

        glEnd();

        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        glPopMatrix();
    }
    CATCH("drawing bonds")
}

UnsignedIntType CIFWindowGL::CreateListOfDrawing(const double LengthUnit)
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

        if (ShowCIFSize != 0)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            DrawAtoms(LengthUnit, ShowCIFSize, true);

            const auto stop_time = chrono::high_resolution_clock::now();

            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of drawing Atoms has taken time: ","executing printing duration_time")));
        }

        glEndList();
    }
    CATCH("creating list of drawing")

    return DrawList;
}

void CIFWindowGL::FullDrawStage()
{
    try
    {
        RefreshListOfDrawing = true;
        DrawStage();
    }
    CATCH("full draw stage")
}

void CIFWindowGL::StartTimerEvent()
{
    try
    {
        CIFDataFileObjectPointer->ChosenStructureIndex = 0;

        if (SetTimer(HandleWindow, 1, 50, nullptr) == 0 )
            MessageBox(HandleWindow, "Nie udało się ustawić timera", "", MB_OK | MB_ICONERROR);
    }
    CATCH("start timer event")
}

void CIFWindowGL::FilmTimerEvent()
{
    try
    {
        if (CIFDataFileObjectPointer->ChosenStructureIndex < CIFDataFileObjectPointer->GetNumberOfStructures() - 1)
        {
            CIFDataFileObjectPointer->ChosenStructureIndex++;
            FullDrawStage();
        }
        else
            KillTimer(HandleWindow, 1);
    }
    CATCH("timer event")
}

void CIFWindowGL::ChangeAtomsSize()
{
    try
    {
        float NewSize = -1.0f;

        if (ShowCIFSize == 0.0f)
            NewSize = 0.25f;
        if (ShowCIFSize == 0.25f)
            NewSize = 0.5f;
        if (ShowCIFSize == 0.5f)
            NewSize = 1.0f;
        if (ShowCIFSize == 1.0f)
            NewSize = 0.0f;
        if (NewSize == -1.0f)
            NewSize = 1.0f;

        ShowCIFSize = NewSize;

        FullDrawStage();
    }
    CATCH("changing Atoms size")
}

void CIFWindowGL::ShowNextStructure()
{
    try
    {
        if (CIFDataFileObjectPointer->ChosenStructureIndex < CIFDataFileObjectPointer->GetNumberOfStructures() - 1)
            CIFDataFileObjectPointer->ChosenStructureIndex++;
        FullDrawStage();
    }
    CATCH("showing next structure")
}

void CIFWindowGL::ShowPrevStructure()
{
    try
    {
        if (CIFDataFileObjectPointer->ChosenStructureIndex > 0)
            CIFDataFileObjectPointer->ChosenStructureIndex--;
        FullDrawStage();
    }
    CATCH("showing previous structure")
}

void CIFWindowGL::ChangeShowOfBonds()
{
    try
    {
        ShowBonds = !ShowBonds;
        FullDrawStage();
    }
    CATCH("changing showing of bonds")
}

void CIFWindowGL::ChangeShowOfSpheres()
{
    try
    {
        ShowSpheres = !ShowSpheres;
        FullDrawStage();
    }
    CATCH("changing showing of bonds")
}

void CIFWindowGL::KeyboardKeyDownPressedEvent(WPARAM wParam)
{
    try
    {
        switch (wParam)
        {
            case VK_F1: MessageBox(HandleWindow, "Shortcut Keys: \nZ - show / hide bonds \nX - size of Atoms \nC - projection mode", "CIF VIEWER", MB_OK); break;
            case VK_F2: MessageBox(HandleWindow, ChosenAtomDescription.c_str(), "CIF VIEWER", MB_OK); break;
            case 'Z': ChangeShowOfBonds(); break;
            case 'J': ChangeShowOfSpheres(); break;
            case 'X': ChangeAtomsSize(); break;
            case 'N': ShowNextStructure(); break;
            case 'M': ShowPrevStructure(); break;
            case 'F': StartTimerEvent(); break;
            default: break;
        }

        DrawStage();
    }
    CATCH("CIFWindow key down pressed event")
}

void CIFWindowGL::MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        if (ShowCIFSize == 0)
            return;

        POINT CursorMouseIndex;
        CursorMouseIndex.x = LOWORD(lParam);
        CursorMouseIndex.y = HIWORD(lParam);

        IntType PrevIndex = ChosenAtomIndex;
        IntType NextIndex = ChooseAtom(CursorMouseIndex);

        if (NextIndex != -1)
            ChosenAtomIndex = NextIndex;

        if (PrevIndex != -1 || ChosenAtomIndex != -1)
            DrawStage();
    }
    CATCH("mouse move event")
}

LRESULT CIFWindowGL::WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam)
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
    CATCH("execution of CIFWindowGL::WndProc")

    return Result;
}

IntType CIFWindowGL::ChooseAtom(POINT MouseCursorPosition)
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
        if (IsometricProjection == false)
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
    CATCH("choosing Atom")

    return -1;
}

#pragma endregion

UnsignedIntType CreateFontGeneral1(HWND HandleWindow, const char* FontName, IntType HeightInPixels, bool Bold, bool Italics, IntType FirstCharCode, IntType LastCharCode)
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

void CIFWindowGL::DrawActors()
{
    try
    {
        ShowRenderingFrequency();

        const auto start_time = chrono::high_resolution_clock::now();

        const double LengthUnit = 0.1;

        static GLuint DrawingList = 0;
        if (glIsList(DrawingList) == false || RefreshListOfDrawing == true)
        {
            if (RefreshListOfDrawing == true)
                glDeleteLists(DrawingList, 1);
            DrawingList = CreateListOfDrawing(LengthUnit);
            RefreshListOfDrawing = false;
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of X1 has taken time: ","executing printing duration_time")));

        static UnsignedIntType BitmapFont = 0;
        if (BitmapFont == 0)
            BitmapFont = CreateFontGeneral1(HandleWindow, "Calibri", 20, true, false, 32, 255);
        const auto start_time2 = chrono::high_resolution_clock::now();

        glCallList(DrawingList);

        const auto stop_time2 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time2, stop_time2, "Execution of X2 has taken time: ","executing printing duration_time")));


        const auto start_time3 = chrono::high_resolution_clock::now();

        int RenderingMode = 0;
        glGetIntegerv(GL_RENDER_MODE, &RenderingMode);
        if (RenderingMode == GL_RENDER)
        {
            if (ChosenAtomIndex >= 0 && ShowCIFSize > 0)
            {
                glColor3f(1, 1, 0);
                DrawChosenAtom(LengthUnit, 1.05 * ShowCIFSize, ChosenAtomIndex, false);
                glColor3f(5, 0, 5);
                DrawChosenAtomDescription(ChosenAtomIndex, BitmapFont);
            }
        }

        const auto stop_time3 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time3, stop_time3, "Execution of X3 has taken time: " ,"executing printing duration_time")));
    }
    CATCH("drawing actors")
}

#pragma region SourceOfLights

void CIFWindowGL::LightSources()
{
    BackgroundLightIntensity = 0.5f;
    MilkyBulb(0.5f);
}

void CIFWindowGL::MilkyBulb(float jasnosc)
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