#pragma once

#ifndef CELL_ENGINE_WINDOW_GL_H
#define CELL_ENGINE_WINDOW_GL_H

#include <windows.h>

#include <gl\gl.h>
#include <gl\glu.h>

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
	void SetStage(bool IsometricProjection = false);
	void DrawStage();
	virtual void DrawActors() = 0;
	void ShowRenderingFrequency();
public:
	WindowGL() : Window(), HandleRC(nullptr), HandleDC(nullptr), ControlCameraByUser(true), MouseCursorInitialPosition(POINT()), CameraR(10), CameraCelPhi(0), CameraCelTheta(0), CameraX(0), CameraY(0), CameraZ(0), ArcBall(nullptr), BackgroundLightIntensity(0.5f)
	{
		InitArcBall();
	};
	LRESULT WndProc(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam) override;
protected:
	bool ControlCameraByUser;
private:
	POINT MouseCursorInitialPosition;
	float CameraR;
	float CameraCelPhi, CameraCelTheta;
	float CameraX, CameraY, CameraZ;
	float CameraPosition[3];
private:
	Matrix4fT Transform;
	Matrix3fT LastRot;
	Matrix3fT ThisRot;
	ArcBallT* ArcBall;
	Point2fT MousePt;
public:
	void InitArcBall();
	~WindowGL();
public:
	bool FreeRotationActive;
	bool FreeRotationBlanking;
	Quat4fT FreeRotationQuaternionOfRotation;
public:
	float BackgroundLightIntensity;
private:
	void Lighting();
protected:
	virtual void LightSources() = 0;
public:
};

#endif
