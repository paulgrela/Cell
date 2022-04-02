
#include <cstdint>

#include "ExceptionsMacro.h"

#include "CellEngineWindowGL.h"
#include "AdditionalFunctions.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"

#include "CellConstants.h"
#include "CellEnginePDBWindowGL.h"
#include "CellEngineCIFWindowGL.h"

using namespace std;

void InitializeLoggerManagerParameters()
{
    try
    {
        using namespace string_utils;

        LoggersManagerObject.InitializeFilesNames({ "AllMessages" });
        LoggersManagerObject.InitializeSelectiveWordsFunctions({ [](const string& s) { return true; } });
        LoggersManagerObject.InitializePrintingParameters(true, true, false, false, false, false, false, true, true, false, false, false, 10000);
        LoggersManagerObject.InitializeLoggerManagerDataForTask("CELL_RESULTS", ".\\", string("Logs." + GetActualDateTimeStandardCPP(".", ".", ".", ".", ".")), true, 0, function<void(const uint64_t& CurrentThreadId, const uint64_t FileNumber, const string& MessageStr)>());
    }
    CATCH("initializing logger manager parameters")
}

#pragma region Functions WinMain and WndProc

unique_ptr<WindowGL> WindowGLPointer;

bool check_end_str(const string_view FileName, const string_view FileExtension)
{
    return std::equal(FileExtension.rbegin(), FileExtension.rend(), string(FileName).rbegin());
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    try
    {
        InitializeLoggerManagerParameters();
        LoggersManagerObject.Log(STREAM("START CELL"));

                string FileName;
                if (__argc > 1)
                    FileName = __argv[1];
                else
                {
                    MessageBox(nullptr, "Lack of file name in program parameters", "Cell Engine View PDB", MB_OK | MB_ICONWARNING);
                    PostQuitMessage(0);
                }

        if (check_end_str(FileName, ".pdb") == true)
            WindowGLPointer = make_unique<PDBWindowGL>(FileName);
        else
            WindowGLPointer = make_unique<CIFWindowGL>(FileName);

        POINT WindowPosition = { 480,100 };
        POINT WindowSize = { 1600,1200 };

        if (WindowGLPointer->Init(hInstance, WindowPosition, WindowSize) == false)
        {
            MessageBox(nullptr, "Window initiation failed!", "OpenGL PDB Viewer Application", MB_OK | MB_ICONERROR);
            return EXIT_FAILURE;
        }
        else
            return WindowGLPointer->Run();
	}
    CATCH("execution WinMain")

    return 0;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	return WindowGLPointer->WndProc(hWnd, message, wParam, lParam);
}

#pragma endregion

#pragma region Class CWindow

bool Window::Init(HINSTANCE ApplicationHandle, POINT WindowPosition, POINT WindowSize)
{
    try
    {
        char WindowName[] = "OpenGL PDB Viewer Application";

        WNDCLASSEX WindowClassObject;
        WindowClassObject.cbSize = sizeof(WindowClassObject);
        WindowClassObject.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
        WindowClassObject.lpfnWndProc = (WNDPROC)::WndProc;
        WindowClassObject.cbClsExtra = 0;
        WindowClassObject.cbWndExtra = 0;
        WindowClassObject.hInstance = ApplicationHandle;
        WindowClassObject.hIcon = nullptr;
        WindowClassObject.hIconSm = nullptr;
        WindowClassObject.hCursor = LoadCursor(nullptr, IDC_ARROW);
        WindowClassObject.hbrBackground = nullptr;
        WindowClassObject.lpszMenuName = nullptr;
        WindowClassObject.lpszClassName = WindowName;

        if (RegisterClassEx(&WindowClassObject) == 0)
            return false;

        DWORD WindowStyle = WS_OVERLAPPEDWINDOW;

        HandleWindow = CreateWindow(WindowName, WindowName, WindowStyle, WindowPosition.x, WindowPosition.y, WindowSize.x, WindowSize.y, nullptr, nullptr, ApplicationHandle, nullptr );

        if (HandleWindow == nullptr)
            return false;

        ShowWindow(HandleWindow, SW_SHOW);
        UpdateWindow(HandleWindow);
	}
    CATCH("initiation of window")

	return true;
};

WPARAM Window::Run()
{
    MSG msg;

    try
    {
        while (GetMessage(&msg, nullptr, 0, 0))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
	}
    CATCH("running of window")

	return msg.wParam;
}

void Window::ResizeWindow()
{
    try
    {
        RECT rect;
        GetClientRect(HandleWindow, &rect);
        UserAreaWidth = rect.right - rect.left;
        UserAreaHeight = rect.bottom - rect.top;
    }
    CATCH("resizing window")
}

LRESULT Window::WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
    try
    {
        switch (Message)
        {
            case WM_DESTROY: PostQuitMessage(0); break;
            case WM_SIZE: ResizeWindow(); break;
            default: return DefWindowProc(hWnd, Message, wParam, lParam);
        }
	}
    CATCH("executing Window::WndProc")

	return 0L;
}

#pragma endregion

#pragma region Class CWindowGL

void WindowGL::ChangeProjectionType()
{
    try
    {
        IsometricProjection = !IsometricProjection;
        SetStage();
    }
    CATCH("changing projection type")
}

void WindowGL::ShiftCameraCloser()
{
    try
    {
        CameraZ += cosDegf(CameraCelPhi) * CameraShift;
        CameraX -= sinDegf(CameraCelPhi) * CameraShift;
    }
    CATCH("shifting camera up")
}

void WindowGL::ShiftCameraFarer()
{
    try
    {
        CameraZ -= cosDegf(CameraCelPhi) * CameraShift;
        CameraX += sinDegf(CameraCelPhi) * CameraShift;
    }
    CATCH("shifting camera down")
}

void WindowGL::ShiftCameraRight()
{
    try
    {
        CameraZ -= sinDegf(CameraCelPhi) * CameraShift;
        CameraX -= cosDegf(CameraCelPhi) * CameraShift;
    }
    CATCH("shifting camera right")
}

void WindowGL::ShiftCameraLeft()
{
    try
    {
        CameraZ += sinDegf(CameraCelPhi) * CameraShift;
        CameraX += cosDegf(CameraCelPhi) * CameraShift;
    }
    CATCH("shifting camera left")
}

void WindowGL::KeyboardKeyDownPressedEvent(WPARAM wParam)
{
    try
    {
        switch (wParam)
        {
            case VK_ESCAPE: SendMessage(HandleWindow, WM_DESTROY, 0, 0); break;
            case 'C': ChangeProjectionType(); break;
            case 'W':
            case VK_UP: ShiftCameraCloser(); break;
            case 'S':
            case VK_DOWN: ShiftCameraFarer(); break;
            case 'D':
            case VK_RIGHT: ShiftCameraRight(); break;
            case 'A':
            case VK_LEFT: ShiftCameraLeft(); break;

            default: break;
        }

        if (wParam >= '0' && wParam <= '7')
        {
            GLenum Light = GL_LIGHT0;
            if (wParam == 1)
                Light = GL_LIGHT1;

            if (glIsEnabled(Light))
                glDisable(Light);
            else
                glEnable(Light);
        }

        DrawStage();
    }
    CATCH("WindowGL key down pressed event")
}

UnsignedIntType WindowGL::MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        MousePt.s.X = (GLfloat)LOWORD(lParam);
        MousePt.s.Y = (GLfloat)HIWORD(lParam);
        LastRot = ThisRot;
        ArcBall->click(&MousePt);
    }
    CATCH("mouse left button down event")

    return 0;
}

UnsignedIntType WindowGL::MouseRightButtonDownEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        MouseCursorInitialPosition.x = LOWORD(lParam);
        MouseCursorInitialPosition.y = HIWORD(lParam);
    }
    CATCH("mouse left button down event")

    return 0;
}

UnsignedIntType WindowGL::MouseMoveEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        if (wParam & (MK_LBUTTON | MK_RBUTTON | MK_MBUTTON))
        {
            POINT MouseCursorCurrentPosition = { LOWORD(lParam),HIWORD(lParam) };
            POINT MouseCursorShift = { MouseCursorCurrentPosition.x - MouseCursorInitialPosition.x, MouseCursorCurrentPosition.y - MouseCursorInitialPosition.y };

            if (MouseCursorShift.x == 0 && MouseCursorShift.y == 0)
                return 0;

            if (wParam & MK_LBUTTON)
            {
                MousePt.s.X = (GLfloat)MouseCursorCurrentPosition.x;
                MousePt.s.Y = (GLfloat)MouseCursorCurrentPosition.y;
                Quat4fT ThisQuat;
                ArcBall->drag(&MousePt, &ThisQuat);
                Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);
                Matrix3fMulMatrix3f(&ThisRot, &LastRot);
                Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);
            }
            if (wParam & MK_RBUTTON)
            {
                const float MouseSensitivity = 75.0f;
                CameraX += MouseCursorShift.x / MouseSensitivity;
                CameraY -= MouseCursorShift.y / MouseSensitivity;
            }
            if (wParam & MK_MBUTTON)
            {
                const float MouseSensitivity = 5.0f;
                CameraCelPhi += MouseCursorShift.x / MouseSensitivity;
                float ChangeAimTheta = MouseCursorShift.y / MouseSensitivity;
                if (fabs(CameraCelTheta + ChangeAimTheta) < 90)
                    CameraCelTheta += ChangeAimTheta;
            }

            MouseCursorInitialPosition.x = LOWORD(lParam);
            MouseCursorInitialPosition.y = HIWORD(lParam);

            DrawStage();
        }
    }
    CATCH("mouse move event")

    return 0;
}

UnsignedIntType WindowGL::MouseWheelEvent(WPARAM wParam, LPARAM lParam)
{
    try
    {
        const float MouseSensitivity = 10.0f;
        auto RollingPositionChange = static_cast<short>(HIWORD(wParam));
        CameraR *= 1 + static_cast<float>(RollingPositionChange) / abs(static_cast<float>(RollingPositionChange)) / MouseSensitivity;

        DrawStage();
    }
    CATCH("mouse wheel event")

    return 0;
}

void WindowGL::CreateWindowEvent(HWND hWnd)
{
    try
    {
        InitWGL(hWnd);
        SetStage();
    }
    CATCH("create window event")
}

void WindowGL::PaintEvent(HWND hWnd)
{
    try
    {
        DrawStage();
        ValidateRect(hWnd, nullptr);
    }
    CATCH("paint event")
}

LRESULT WindowGL::WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
    IntType Result;

    try
    {
        Result = Window::WndProc(hWnd, Message, wParam, lParam);

        switch (Message)
        {
            case WM_CREATE: CreateWindowEvent(hWnd); break;
            case WM_DESTROY: DeleteWGL(); break;
            case WM_SIZE: SetStage(); break;
            case WM_PAINT: PaintEvent(hWnd); break;
            case WM_KEYDOWN: KeyboardKeyDownPressedEvent(wParam); break;
            case WM_LBUTTONDOWN: Result = MouseLeftButtonDownEvent(wParam, lParam); break;
            case WM_RBUTTONDOWN: Result = MouseRightButtonDownEvent(wParam, lParam); break;
            case WM_MOUSEMOVE: Result = MouseMoveEvent(wParam, lParam); break;
            case WM_MOUSEWHEEL: Result = MouseWheelEvent(wParam, lParam); break;

            default: break;
        }
    }
    CATCH("executing WindowGL::WndProc")

	return Result;
}

bool WindowGL::InitWGL(HWND HandleWindow)
{
    try
    {
        HandleDC = ::GetDC(HandleWindow);

        if (!SetPixelFormatWindow(HandleDC))
            return false;

        HandleRC = wglCreateContext(HandleDC);

        if (HandleRC == nullptr)
            return false;

        if (!wglMakeCurrent(HandleDC, HandleRC))
            return false;
    }
    CATCH("initating WGL")

	return true;
}

void WindowGL::DeleteWGL()
{
    try
    {
        wglMakeCurrent(nullptr, nullptr);
        wglDeleteContext(HandleRC);
        ::ReleaseDC(HandleWindow, HandleDC);
    }
    CATCH("deleting WGL")
}

bool WindowGL::SetPixelFormatWindow(HDC HandleDC)
{
    try
    {
        PIXELFORMATDESCRIPTOR PixelFormatDescription;
        ZeroMemory(&PixelFormatDescription, sizeof(PixelFormatDescription));
        PixelFormatDescription.nVersion = 1;
        PixelFormatDescription.dwFlags = PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW | PFD_DOUBLEBUFFER;
        PixelFormatDescription.iPixelType = PFD_TYPE_RGBA;
        PixelFormatDescription.cColorBits = 32;
        PixelFormatDescription.cDepthBits = 16;
        PixelFormatDescription.cStencilBits = 1;
        PixelFormatDescription.iLayerType = PFD_MAIN_PLANE;
        int PixelFormat = ChoosePixelFormat(HandleDC, &PixelFormatDescription);

        if (PixelFormat == 0)
            return false;

        if (!SetPixelFormat(HandleDC, PixelFormat, &PixelFormatDescription))
            return false;
    }
    CATCH("setting pixel fromat window")

	return true;
}

void WindowGL::SetStage()
{
    try
    {
        glViewport(0, 0, static_cast<GLsizei>(UserAreaWidth), static_cast<GLsizei>(UserAreaHeight));
        glClearColor(0.0, 0.0, 0.0, 1.0);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float Coordinates = static_cast<float>(UserAreaHeight) / static_cast<float>(UserAreaWidth);
        if (IsometricProjection == false)
            glFrustum(-0.1, 0.1, Coordinates * -0.1, Coordinates * 0.1, 0.3, 100.0);
        else
            glOrtho(-3, 3, Coordinates * -3, Coordinates * 3, 0.3, 100.0);
        glMatrixMode(GL_MODELVIEW);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        ArcBall->setBounds(static_cast<float>(UserAreaWidth), static_cast<float>(UserAreaHeight));

        Lighting();
    }
    CATCH("setting stage")
}

void WindowGL::SetCamera()
{
    try
    {
        glLoadIdentity();
        glRotatef(CameraCelPhi, 0, 1, 0);
        glRotatef(CameraCelTheta, cosDegf(CameraCelPhi), 0, sinDegf(CameraCelPhi));
        glTranslatef(CameraX, CameraY, CameraZ);

        glTranslatef(0, 0, -CameraR);

        glMultMatrixf(Transform.M);
    }
    CATCH("setting camera")
}

void WindowGL::DrawStage()
{
    try
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        glLoadIdentity();

        ShowRenderingFrequency();

        //const auto start_time1 = chrono::high_resolution_clock::now();
        SetCamera();
        //const auto stop_time1 = chrono::high_resolution_clock::now();
        //LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of setting camera has taken time: ","executing printing duration_time")));

        const auto start_time2 = chrono::high_resolution_clock::now();

        DrawActors();

        const auto stop_time2 = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time2, stop_time2, "Execution of drawing actors has taken time: ","executing printing duration_time")));

        SwapBuffers(HandleDC);
    }
    CATCH("drawing stage")
}

void WindowGL::ShowRenderingFrequency()
{
    try
    {
        static uint64_t OldTickNumber = GetTickCount();
        uint64_t NewTickNumber = GetTickCount();

        if (NewTickNumber == OldTickNumber)
            return;

        double f = 1E3 / (static_cast<double>(NewTickNumber) - static_cast<double>(OldTickNumber));
        f = floor(10.0 * f) / 10.0;
        OldTickNumber = NewTickNumber;

        SetWindowText(HandleWindow, string("OpenGL PDB Viewer Application - OpenGL Version = " + string((char*)glGetString(GL_VERSION)) + " GLU = " + string((char*)gluGetString(GLU_VERSION)) + " [ Rendering Frequency = " + to_string(f) + " Hz]").c_str());
    }
    CATCH("showing rendering frequency")
}

void WindowGL::Lighting()
{
    try
    {
        glEnable(GL_LIGHTING);

        const float BackgroundColor[] = { BackgroundLightIntensity, BackgroundLightIntensity, BackgroundLightIntensity };
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, BackgroundColor);

        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        glPushMatrix();
        glLoadIdentity();
        LightSources();
        glPopMatrix();

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
    CATCH("lighting")
}

void WindowGL::InitArcBall()
{
    try
    {
        Matrix3fSetIdentity(&LastRot);
        Matrix3fSetIdentity(&ThisRot);

        for (IntType Index1 = 0; Index1 < 4; Index1++)
            for (IntType Index2 = 0; Index2 < 4; Index2++)
                Transform.M[Index1 + 4 * Index2] = (Index1 == Index2) ? 1.0f : 0.0f;

        ArcBall = make_unique<ArcBallT>(ArcBallT(640.0f, 480.0f));
	}
    CATCH("initating arc ball")
}

#pragma endregion