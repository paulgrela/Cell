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
	long UserAreaWidth;
	long UserAreaHeight;
	HWND HandleWindow;
public:
	Window() : HandleWindow(NULL) {};
	bool Init(HINSTANCE ApplicationHandle, POINT WindowPosition, POINT WindowSize);
	WPARAM Run();
	virtual LRESULT WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
private:
	bool ChangeResolution(long Width, long Height, long ColorsDepth = 32) const;
};

#include "ArcBall.h"

class WindowGL : public Window
{
private:
	HGLRC HandleRC; 
	HDC HandleDC; 
	bool UstalFormatPikseli(HDC HandleDC) const;
	bool InitWGL(HWND HandleWindow);
	void DeleteWGL();
protected:
	void SetCamera();
	void SetStage(bool IsometricProjection = false);
	void DrawStage();
	virtual void DrawActors() = 0;
	void ShowRenderingFrequency();
public:
	WindowGL() : Window(), HandleRC(NULL), HandleDC(NULL), ControlCameraByUser(true), MouseCursorInitialPosition(POINT()), CameraR(10), CameraCelPhi(0), CameraCelTheta(0), CameraX(0), CameraY(0), CameraZ(0), ArcBall(NULL), BackgroundLightIntensity(0.5f)
	{
		InitArcBall();
	};
	virtual LRESULT WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);

protected:
	bool ControlCameraByUser;
private:
	POINT MouseCursorInitialPosition;
	float CameraR;
	float CameraCelPhi, CameraCelTheta;
	float CameraX, CameraY, CameraZ;
	float CameraPosition[3];
public:
	float* SetCameraPosition(float* Buffer) const;

private:
	Matrix4fT Transform;
	Matrix3fT LastRot;
	Matrix3fT ThisRot;
	ArcBallT* ArcBall;
	Point2fT MousePt;
public:
	void InitArcBall();
	~WindowGL(); //usuwanie obiektu ArcBall

//swobodne obroty
	bool FreeRotationAcitve;
	bool FreeRotationBlanking;
	Quat4fT FreeRotationQuaternionOfRotation;

	//oswietlenie
public:
	float BackgroundLightIntensity;
private:
	void Lighting();
protected:
	virtual void LightSources() = 0;
public:
};

#endif
