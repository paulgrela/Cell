#pragma once

#ifndef CELL_ENGINE_WINDOW_GL_H
#define CELL_ENGINE_WINDOW_GL_H

#include <windows.h>

#include <gl\gl.h>
#include <gl\glu.h>

#include <memory>

#include "CellEngineTypes.h"

class Window
{
protected:
	IntType UserAreaWidth;
	IntType UserAreaHeight;
	HWND HandleWindow;
public:
	Window() : HandleWindow(nullptr) {};
	bool Init(HINSTANCE ApplicationHandle, POINT WindowPosition, POINT WindowSize);
	static WPARAM Run();
	virtual LRESULT WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam);
    void ResizeWindow();
};

#include "ArcBall.h"

class WindowGL : public Window
{
private:
	HGLRC HandleRC; 
	HDC HandleDC; 
	static bool SetPixelFormatWindow(HDC HandleDC);
	bool InitWGL(HWND HandleWindow);
	void DeleteWGL();
protected:
	void SetCamera();
	void SetStage();
	void DrawStage();
	virtual void DrawActors() = 0;
	void ShowRenderingFrequency();
public:
	WindowGL() : Window(), HandleRC(nullptr), HandleDC(nullptr), IsometricProjection(false), MouseCursorInitialPosition(POINT()), CameraR(10), CameraCelPhi(0), CameraCelTheta(0), CameraX(0), CameraY(0), CameraZ(0), ArcBall(nullptr), BackgroundLightIntensity(0.5f)
	{
		InitArcBall();
	};
	LRESULT WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam) override;
private:
    void KeyboardKeyDownPressedEvent(WPARAM wParam);
    UnsignedIntType MouseLeftButtonDownEvent(WPARAM wParam, LPARAM lParam);
    UnsignedIntType MouseRightButtonDownEvent(WPARAM wParam, LPARAM lParam);
    UnsignedIntType MouseMoveEvent(WPARAM wParam, LPARAM lParam);
    UnsignedIntType MouseWheelEvent(WPARAM wParam, LPARAM lParam);
private:
    void ChangeProjectionType();
    void ShiftCameraCloser();
    void ShiftCameraFarer();
    void ShiftCameraRight();
    void ShiftCameraLeft();
private:
    void CreateWindowEvent(HWND hWnd);
    void PaintEvent(HWND hWnd);
protected:
    bool IsometricProjection;
private:
	POINT MouseCursorInitialPosition;
	float CameraR;
	float CameraCelPhi, CameraCelTheta;
	float CameraX, CameraY, CameraZ;
private:
	Matrix4fT Transform;
	Matrix3fT LastRot;
	Matrix3fT ThisRot;
    std::unique_ptr<ArcBallT> ArcBall;
	Point2fT MousePt;
public:
	void InitArcBall();
	~WindowGL() = default;
public:
	float BackgroundLightIntensity;
private:
	void Lighting();
protected:
	virtual void LightSources() = 0;
public:
};

#endif
