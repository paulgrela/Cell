#include "CellEngineWindowGL.h"

#define _USE_MATH_DEFINES
#include <math.h>

#define IDI_GLICON 101

#include <cstdint>

using namespace std;


							double DegToRad(double deg)
							{
								return M_PI * deg / 180.0f;
							}

							double sinDeg(double deg)
							{
								return sin(DegToRad(deg));
							}

							double cosDeg(double deg)
							{
								return cos(DegToRad(deg));
							}

							float sinDegf(float deg)
							{
								return (float)sinDeg(deg);
							}

							float cosDegf(float deg)
							{
								return (float)cosDeg(deg);
							}

							
#pragma region Functions WinMain and WndProc
extern WindowGL* WindowGLPointer;

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	POINT WindowPosition = { 480,100 };
	POINT WindowSize = { 1600,1200 };

	if (!WindowGLPointer->Init(hInstance, WindowPosition, WindowSize))
	{
		MessageBox(NULL, "Window initation failed!", "OpenGL PDB Viewer Application", MB_OK | MB_ICONERROR);
		return EXIT_FAILURE;
	}
	else 
		return WindowGLPointer->Run();
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	return WindowGLPointer->WndProc(hWnd, message, wParam, lParam);
}

#pragma endregion

#pragma region Class CWindow
bool Window::Init(HINSTANCE ApplicationHandle, POINT WindowPosition, POINT WindowSize)
{
	char WindowName[] = "OpenGL PDB Viewer Application";

	WNDCLASSEX WindowClassObject;
	WindowClassObject.cbSize = sizeof(WindowClassObject);
	WindowClassObject.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	WindowClassObject.lpfnWndProc = (WNDPROC)::WndProc; 
	WindowClassObject.cbClsExtra = 0; 
	WindowClassObject.cbWndExtra = 0; 
	WindowClassObject.hInstance = ApplicationHandle; 
	WindowClassObject.hIcon = LoadIcon(ApplicationHandle, MAKEINTRESOURCE(IDI_GLICON));
	WindowClassObject.hIconSm = LoadIcon(ApplicationHandle, MAKEINTRESOURCE(IDI_GLICON));
	WindowClassObject.hCursor = LoadCursor(NULL, IDC_ARROW);
	WindowClassObject.hbrBackground = NULL;
	WindowClassObject.lpszMenuName = NULL;
	WindowClassObject.lpszClassName = WindowName;

	if (RegisterClassEx(&WindowClassObject) == 0) 
		return false;

	bool FullScreenMode = false;

	DWORD WindowStyle = WS_OVERLAPPEDWINDOW;

	if (FullScreenMode)
	{
		WindowPosition.x = 0;
		WindowPosition.y = 0;
		RECT ScreenSize;
		GetWindowRect(GetDesktopWindow(), &ScreenSize);
		WindowSize.x = ScreenSize.right - ScreenSize.left;
		WindowSize.y = ScreenSize.bottom - ScreenSize.top;
		WindowStyle = WS_POPUP;

		if (!ChangeResolution(WindowSize.x, WindowSize.y)) 
			return false;
	}

	HandleWindow = CreateWindow(WindowName, WindowName, WindowStyle, WindowPosition.x, WindowPosition.y, WindowSize.x, WindowSize.y, NULL, NULL, ApplicationHandle, NULL );

	if (HandleWindow == NULL) 
		return false;

	ShowWindow(HandleWindow, SW_SHOW);
	UpdateWindow(HandleWindow);

	return true;
};

WPARAM Window::Run()
{
	MSG msg;
	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return msg.wParam;
}

LRESULT Window::WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
		case WM_DESTROY:
			PostQuitMessage(0);
			break;
		case WM_SIZE:
			RECT rect;
			GetClientRect(hWnd, &rect);
			UserAreaWidth = rect.right - rect.left;
			UserAreaHeight = rect.bottom - rect.top;
			break;
		default: 
			return (DefWindowProc(hWnd, message, wParam, lParam));
	}

	return 0L;
}

bool Window::ChangeResolution(long Width, long Height, long ColorsDepth) const
{
	DEVMODE dmScreenSettings;
	memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));
	dmScreenSettings.dmSize = sizeof(dmScreenSettings);
	dmScreenSettings.dmPelsWidth = Width;
	dmScreenSettings.dmPelsHeight = Height;
	dmScreenSettings.dmBitsPerPel = ColorsDepth;
	dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;
	return ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) == DISP_CHANGE_SUCCESSFUL;
}
#pragma endregion

#pragma region Class CWindowGL
LRESULT WindowGL::WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	const bool FreeRotationOfCamera = false;
	const IntType CameraRotationIdtTimer = 2;
	const IntType TimerInterval = 50;

	long Result = Window::WndProc(hWnd, message, wParam, lParam);

	switch (message)
	{
		case WM_CREATE: 
			InitWGL(hWnd);
			SetStage();
			{
				char Title[1024] = "OpenGL ";
				strcat_s(Title, (char*)glGetString(GL_VERSION));
				strcat_s(Title, ", GLU ");
				strcat_s(Title, (char*)gluGetString(GLU_VERSION));
				SetWindowText(hWnd, Title);
			}
			if (FreeRotationOfCamera)
				if (SetTimer(hWnd, CameraRotationIdtTimer, 50, NULL) == 0)
					MessageBox(hWnd, "Setting timer failed", "", MB_OK | MB_ICONERROR);
			break;

		case WM_TIMER:
			switch (wParam)
			{
				case CameraRotationIdtTimer:
					if (FreeRotationAcitve == true)
					{
						Matrix3fSetRotationFromQuat4f(&ThisRot, &FreeRotationQuaternionOfRotation);
						Matrix3fMulMatrix3f(&ThisRot, &LastRot); //powtorzenie ostatniego obrotu
						Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);
						LastRot = ThisRot;

						//Blanking swobodnych obrotow = zmniejszanie kata kolejnego obrotu
						if (FreeRotationBlanking)
						{
							const float ExtinctionRate = 0.97f;
							FreeRotationQuaternionOfRotation.s.W /= ExtinctionRate;
							FreeRotationQuaternionOfRotation.s.X *= ExtinctionRate;
							FreeRotationQuaternionOfRotation.s.Y *= ExtinctionRate;
							FreeRotationQuaternionOfRotation.s.Z *= ExtinctionRate;
							if (fabs(FreeRotationQuaternionOfRotation.s.X) < 1E-3 && fabs(FreeRotationQuaternionOfRotation.s.Y) < 1E-3 && fabs(FreeRotationQuaternionOfRotation.s.Z) < 1E-3)
							{
								FreeRotationAcitve = false;
							}
						}

						DrawStage();
					}
					break;
			}
			Result = 0;
			break;

		case WM_DESTROY: 
			DeleteWGL();
			KillTimer(HandleWindow, CameraRotationIdtTimer);
			break;

		case WM_SIZE: 
			SetStage();
			break;

		case WM_PAINT: 
			DrawStage();
			ValidateRect(hWnd, NULL);
			break;

		case WM_KEYDOWN:
			switch (wParam)
			{
				case VK_ESCAPE:
					SendMessage(HandleWindow, WM_DESTROY, 0, 0);
					break;

				case VK_OEM_MINUS:
					BackgroundLightIntensity -= 0.01f;
					if (BackgroundLightIntensity < 0) BackgroundLightIntensity = 0;
					Lighting();
					break;

				case VK_OEM_PLUS:
				case '=':
					BackgroundLightIntensity += 0.01f;
					if (BackgroundLightIntensity > 1) BackgroundLightIntensity = 1;
					Lighting();
					break;

				case 'C':
					static bool IsometricProjection = false;
					IsometricProjection = !IsometricProjection;
					SetStage(IsometricProjection);
					break;
			}

			if (wParam >= '0' && wParam <= '7')
			{
				GLenum Light = GL_LIGHT0;
				switch (wParam)
				{
					case '1': Light = GL_LIGHT1; break;
					case '2': Light = GL_LIGHT2; break;
					case '3': Light = GL_LIGHT3; break;
					case '4': Light = GL_LIGHT4; break;
					case '5': Light = GL_LIGHT5; break;
					case '6': Light = GL_LIGHT6; break;
					case '7': Light = GL_LIGHT7; break;
					default: Light = GL_LIGHT0;
				}

				if (glIsEnabled(Light)) 
					glDisable(Light);
				else 
					glEnable(Light);
			}

			DrawStage();
			break;
	}

	//kamera
	if (ControlCameraByUser)
	{
		switch (message)
		{
			case WM_LBUTTONDOWN:

				MousePt.s.X = (GLfloat)LOWORD(lParam);
				MousePt.s.Y = (GLfloat)HIWORD(lParam);
				LastRot = ThisRot;
				ArcBall->click(&MousePt);

				FreeRotationAcitve = false;

			case WM_RBUTTONDOWN:
			case WM_MBUTTONDOWN:
				MouseCursorInitialPosition.x = LOWORD(lParam);
				MouseCursorInitialPosition.y = HIWORD(lParam);
				Result = 0;
				break;

			case WM_MOUSEMOVE:
				if (wParam & (MK_LBUTTON | MK_RBUTTON | MK_MBUTTON))
				{
					POINT MouseCursorCurrentPosition = { LOWORD(lParam),HIWORD(lParam) };
					POINT MouseCursorShift = { MouseCursorCurrentPosition.x - MouseCursorInitialPosition.x, MouseCursorCurrentPosition.y - MouseCursorInitialPosition.y };

					if (MouseCursorShift.x == 0 && MouseCursorShift.y == 0)
						break;

					if (wParam & MK_LBUTTON)
					{
						MousePt.s.X = (GLfloat)MouseCursorCurrentPosition.x;
						MousePt.s.Y = (GLfloat)MouseCursorCurrentPosition.y;
						Quat4fT ThisQuat;
						ArcBall->drag(&MousePt, &ThisQuat);
						Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);
						Matrix3fMulMatrix3f(&ThisRot, &LastRot);
						Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);
						//swobodne obroty
						FreeRotationAcitve = true;
						FreeRotationQuaternionOfRotation = ThisQuat;
					}
					if (wParam & MK_RBUTTON)
					{
						const float czuloscMyszy = 75.0f;
						CameraX += MouseCursorShift.x / czuloscMyszy;
						CameraY -= MouseCursorShift.y / czuloscMyszy;
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
					//MouseCursorPositionInWindow(HandleWindow,&pozycjaPoczatkowaKursoraMyszy);
					DrawStage();
				}
				Result = 0;
				break;

			case WM_MOUSEWHEEL:
			{
				const float MouseSensitivity = 10.0f;
				short RollingPositionChange = (short)HIWORD(wParam);
				//zmiana odleglosci kamery od pocz. ukl. wsp.
				CameraR *= 1 + RollingPositionChange / abs(RollingPositionChange) / MouseSensitivity;
				DrawStage();
				Result = 0;
				break;
			}

			case WM_KEYDOWN:
				const float Shift = 0.1f;
				switch (wParam)
				{
					case 'W':
					case VK_UP:
						CameraZ += cosDegf(CameraCelPhi) * Shift;
						CameraX -= sinDegf(CameraCelPhi) * Shift;
						break;
					case 'S':
					case VK_DOWN:
						CameraZ -= cosDegf(CameraCelPhi) * Shift;
						CameraX += sinDegf(CameraCelPhi) * Shift;
						break;
					case 'A':
					case VK_LEFT:
						CameraZ += sinDegf(CameraCelPhi) * Shift;
						CameraX += cosDegf(CameraCelPhi) * Shift;
						break;
					case 'D':
					case VK_RIGHT:
						CameraZ -= sinDegf(CameraCelPhi) * Shift;
						CameraX -= cosDegf(CameraCelPhi) * Shift;
						break;
				}

				DrawStage();
				break;
		}
	}

	return Result;
}

bool WindowGL::InitWGL(HWND HandleWindow)
{
	HandleDC = ::GetDC(HandleWindow);

	if (!UstalFormatPikseli(HandleDC)) 
		return false; //Utworzenie kontekstu renderowania Index uczynienie go aktywnym

	HandleRC = wglCreateContext(HandleDC);

	if (HandleRC == NULL) 
		return false;

	if (!wglMakeCurrent(HandleDC, HandleRC)) 
		return false;

	return true;
}

void WindowGL::DeleteWGL()
{
	wglMakeCurrent(NULL, NULL);
	wglDeleteContext(HandleRC);
	::ReleaseDC(HandleWindow, HandleDC);
}

bool WindowGL::UstalFormatPikseli(HDC HandleDC) const
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
	IntType formatPikseli = ChoosePixelFormat(HandleDC, &PixelFormatDescription);

	if (formatPikseli == 0) 
		return false;

	if (!SetPixelFormat(HandleDC, formatPikseli, &PixelFormatDescription)) 
		return false;

	return true;
}

// ----------------- OpenGL -----------------

void WindowGL::SetStage(bool IsometricProjection)
{
	glViewport(0, 0, UserAreaWidth, UserAreaHeight); //okno OpenGL = wnetrze formy (domyslnie) 
	glClearColor(0.0, 0.0, 0.0, 1.0); //czarne tlo	

	//ustawienie punktu projekcji 
	glMatrixMode(GL_PROJECTION); //prze³¹czenie na macierz projekcji
	glLoadIdentity();
	//left,right,bottom,top,znear,zfar (clipping) 
	float wsp = UserAreaHeight / (float)UserAreaWidth;
	if (!IsometricProjection)
		glFrustum(-0.1, 0.1, wsp * -0.1, wsp * 0.1, 0.3, 100.0); //mnozenie macierzy rzutowania przez macierz perspektywy - ustalanie frustum 	
		//gluPerspective(RadToDeg(2*atan(wsp*0.1/0.3)),1/wsp,0.3,100.0);
	else
		glOrtho(-3, 3, wsp * -3, wsp * 3, 0.3, 100.0); //rzutowanie rownolegle
	//glScalef(1,-1,1); //do góry-nogami
	glMatrixMode(GL_MODELVIEW); //powrót do macierzy widoku modelu 
	glEnable(GL_DEPTH_TEST); //Z-buffer aktywny = ukrywanie niewidocznych powierzchni 	
	glDepthFunc(GL_LEQUAL);

	ArcBall->setBounds((float)UserAreaWidth, (float)UserAreaHeight);

	Lighting();
}

void WindowGL::SetCamera()
{
	glLoadIdentity();
	glRotatef(CameraCelPhi, 0, 1, 0);
	glRotatef(CameraCelTheta, cosDegf(CameraCelPhi), 0, sinDegf(CameraCelPhi));
	glTranslatef(CameraX, CameraY, CameraZ);

	glTranslatef(0, 0, -CameraR);

	//gluLookAt(-CameraX,-CameraY,-CameraZ+CameraR, 
	//		  -CameraX+CameraR*sinDegf(CameraCelPhi),-CameraY-CameraR*sinDegf(CameraCelTheta),-CameraZ,
	//	      0,1,0);

	//wykonanie transformacji
	glMultMatrixf(Transform.M);

	//obliczenie polozenia kamery
	CameraPosition[0] = -CameraX - Transform.s.XZ * CameraR;
	CameraPosition[1] = -CameraY - Transform.s.YZ * CameraR;
	CameraPosition[2] = -CameraZ - Transform.s.ZZ * CameraR;

	
	//glLoadIdentity();
	//gluLookAt(-CameraPosition[0],-CameraPosition[1],-CameraPosition[2], 0,0,0, 0,1,0);
}

float* WindowGL::SetCameraPosition(float* Buffer) const
{
	for (IntType Index = 0; Index < 3; Index++) 
		Buffer[Index] = CameraPosition[Index];
	return Buffer;
}

void WindowGL::DrawStage()
{
	//Przygotowanie bufora 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); //czysci bufory 
	glLoadIdentity(); //macierz model-widok = macierz jednostkowa 

	ShowRenderingFrequency();

	SetCamera();
	DrawActors();

	//Z bufora na ekran 
	SwapBuffers(HandleDC);
}

void WindowGL::ShowRenderingFrequency()
{
	static uint64_t OldTickNumber = GetTickCount();
	uint64_t NewTickNumber = GetTickCount();

	if (NewTickNumber == OldTickNumber) 
		return;
	double f = 1E3 / (NewTickNumber - OldTickNumber);
	f = floor(10.0 * f) / 10.0;
	OldTickNumber = NewTickNumber;

	char Buffer[256];
	SetWindowText(HandleWindow, strcat(_gcvt(f, 10, Buffer), "Hz"));
}


//oswietlenie
void WindowGL::Lighting()
{
	glEnable(GL_LIGHTING); //wlaczenie systemu oswietlania

	//Light tla			
	const float kolor_tla[] = { BackgroundLightIntensity, BackgroundLightIntensity, BackgroundLightIntensity };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, kolor_tla);

	//material
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

	//zrodla swiatla
	glPushMatrix();
	glLoadIdentity();
	LightSources();
	glPopMatrix();

	//mieszanie kolorow
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void WindowGL::InitArcBall()
{
	Matrix3fSetIdentity(&LastRot);
	Matrix3fSetIdentity(&ThisRot);
	//Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);
	for (IntType Index1 = 0; Index1 < 4; Index1++)
		for (IntType Index2 = 0; Index2 < 4; Index2++)
			Transform.M[Index1 + 4 * Index2] = (Index1 == Index2) ? 1.0f : 0.0f;
	ArcBall = new ArcBallT(640.0f, 480.0f);
}

WindowGL::~WindowGL()
{
	delete ArcBall;
}

#pragma endregion