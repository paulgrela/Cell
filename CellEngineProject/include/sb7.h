
#ifndef __SB7_H__
#define __SB7_H__

#ifdef WIN32
    #pragma once
    #define _CRT_SECURE_NO_WARNINGS 1

    #define WIN32_LEAN_AND_MEAN 1
    #include <Windows.h>
#else
    #include <unistd.h>
    #define Sleep(t) sleep(t)
#endif

#include "GL/gl3w.h"

#define GLFW_NO_GLU 1
#define GLFW_INCLUDE_GLCOREARB 1

#include "GLFW/glfw3.h"

#include "sb7ext.h"

#include <cstdio>
#include <cstring>
#include <cmath>

namespace sb7
{
    class OpenGLApplication
    {
    private:
        static void APIENTRY debug_callback(GLenum Source, GLenum Type, GLuint Id, GLenum Severity, GLsizei Length, const GLchar* Message, GLvoid* UserParam);

    public:
        OpenGLApplication() = default;

        virtual ~OpenGLApplication() = default;

        virtual void Run(sb7::OpenGLApplication* OpenGLApplicationObjectParameter)
        {
            InitExternalData();

            bool Running = true;
            OpenGLApplicationObject = OpenGLApplicationObjectParameter;

            if (!glfwInit())
            {
                fprintf(stderr, "Failed to initialize GLFW\n");
                return;
            }

            Init();

            glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, Info.MajorVersion);
            glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, Info.MinorVersion);

            #ifndef _DEBUG
            if (Info.Flags.Debug)
            #endif
                glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

            if (Info.Flags.Robust)
                glfwWindowHint(GLFW_CONTEXT_ROBUSTNESS, GLFW_LOSE_CONTEXT_ON_RESET);

            glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
            glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
            glfwWindowHint(GLFW_SAMPLES, Info.Samples);
            glfwWindowHint(GLFW_STEREO, Info.Flags.Stereo ? GL_TRUE : GL_FALSE);

    //        if (Info.Flags.FullScreen)
    //        {
    //            if (Info.WindowWidth == 0 || Info.WindowHeight == 0)
    //            {
    //                GLFWvidmode mode;
    //                glfwGetDesktopMode(&mode);
    //                Info.WindowWidth = mode.Width;
    //                Info.WindowHeight = mode.Height;
    //            }
    //
    //            glfwOpenWindow(Info.WindowWidth, Info.WindowHeight, 8, 8, 8, 0, 32, 0, GLFW_FULLSCREEN);
    //            glfwSwapInterval((int)Info.flags.vsync);
    //        }
    //        else
            {
                Window = glfwCreateWindow(Info.WindowWidth, Info.WindowHeight, Info.Title, Info.Flags.FullScreen ? glfwGetPrimaryMonitor() : NULL, NULL);
                if (!Window)
                {
                    fprintf(stderr, "Failed to open window\n");
                    return;
                }
                glfwSetWindowPos(Window, 480, 100);
            }

            glfwMakeContextCurrent(Window);

            glfwSetWindowSizeCallback(Window, glfw_onResize);
            glfwSetKeyCallback(Window, glfw_onKey);
            glfwSetMouseButtonCallback(Window, glfw_onMouseButton);
            glfwSetCursorPosCallback(Window, glfw_onMouseMove);
            glfwSetScrollCallback(Window, glfw_onMouseWheel);
            if (!Info.Flags.Cursor)
                glfwSetInputMode(Window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

            // Info.flags.Stereo = (glfwGetWindowParam(GLFW_STEREO) ? 1 : 0);

            gl3wInit();

            #ifdef _DEBUG
            fprintf(stderr, "VENDOR: %s\n", (char *)glGetString(GL_VENDOR));
            fprintf(stderr, "VERSION: %s\n", (char *)glGetString(GL_VERSION));
            fprintf(stderr, "RENDERER: %s\n", (char *)glGetString(GL_RENDERER));
            #endif

            if (Info.Flags.Debug)
            {
                if (gl3wIsSupported(4, 3))
                {
                    glDebugMessageCallback((GLDEBUGPROC)debug_callback, this);
                    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
                }
                else if (sb6IsExtensionSupported("GL_ARB_debug_output"))
                {
                    glDebugMessageCallbackARB((GLDEBUGPROC)debug_callback, this);
                    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
                }
            }

            StartUp();

            do
            {
                Render(glfwGetTime());

                glfwSwapBuffers(Window);
                glfwPollEvents();

                Running &= (glfwGetKey(Window, GLFW_KEY_ESCAPE) == GLFW_RELEASE);
                Running &= (glfwWindowShouldClose(Window) != GL_TRUE);
            }
            while (Running);

            ShutDown();

            glfwDestroyWindow(Window);
            glfwTerminate();
        }

        virtual void Init()
        {
            strcpy(Info.Title, "OpenGL Window");
            Info.WindowWidth = 1600;
            Info.WindowHeight = 1200;
            #ifdef __APPLE__
            Info.MajorVersion = 3;
            Info.MinorVersion = 2;
            #else
            Info.MajorVersion = 4;
            Info.MinorVersion = 3;
            #endif
            Info.Samples = 0;
            Info.Flags.All = 0;
            Info.Flags.Cursor = 1;
            #ifdef _DEBUG
            Info.flags.Debug = 1;
            #endif
        }

        virtual void InitExternalData()
        {
        }

        virtual void StartUp()
        {
        }

        virtual void Render(double CurrentTime)
        {
        }

        virtual void ShutDown()
        {
        }

        void SetWindowTitle(const char* Title)
        {
            glfwSetWindowTitle(Window, Title);
        }

        virtual void OnResize(int Width, int Height)
        {
            Info.WindowWidth = Width;
            Info.WindowHeight = Height;
        }

        virtual void OnKey(int Key, int Action)
        {
        }

        virtual void OnMouseButton(int Button, int Action)
        {
        }

        virtual void OnMouseMove(int X, int Y)
        {
        }

        virtual void OnMouseWheel(int Pos)
        {
        }

        virtual void OnDebugMessage(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message)
        {
            #ifdef _WIN32
            OutputDebugStringA(message);
            OutputDebugStringA("\n");
            #endif
        }

        void GetMousePosition(int& X, int& Y)
        {
            double DX, DY;
            glfwGetCursorPos(Window, &DX, &DY);

            X = static_cast<int>(floor(DX));
            Y = static_cast<int>(floor(DY));
        }

    public:
        struct APPINFO
        {
            char Title[128];
            int WindowWidth;
            int WindowHeight;
            int MajorVersion;
            int MinorVersion;
            int Samples;
            union
            {
                struct
                {
                    unsigned int FullScreen : 1;
                    unsigned int VSync : 1;
                    unsigned int Cursor : 1;
                    unsigned int Stereo : 1;
                    unsigned int Debug : 1;
                    unsigned int Robust : 1;
                };
                unsigned int All;
            }
            Flags;
        };

    protected:
        APPINFO Info;
        static sb7::OpenGLApplication* OpenGLApplicationObject;
        GLFWwindow* Window;

        static void glfw_onResize(GLFWwindow* Window, int w, int h)
        {
            OpenGLApplicationObject->OnResize(w, h);
        }

        static void glfw_onKey(GLFWwindow* Window, int Key, int ScanCode, int Action, int Mods)
        {
            OpenGLApplicationObject->OnKey(Key, Action);
        }

        static void glfw_onMouseButton(GLFWwindow* Window, int Button, int Action, int Mods)
        {
            OpenGLApplicationObject->OnMouseButton(Button, Action);
        }

        static void glfw_onMouseMove(GLFWwindow* Window, double X, double Y)
        {
            OpenGLApplicationObject->OnMouseMove(static_cast<int>(X), static_cast<int>(Y));
        }

        static void glfw_onMouseWheel(GLFWwindow* Window, double XOffset, double YOffset)
        {
            OpenGLApplicationObject->OnMouseWheel(static_cast<int>(YOffset));
        }

        void setVsync(bool Enable)
        {
            Info.Flags.VSync = Enable ? 1 : 0;
            glfwSwapInterval((int)Info.Flags.VSync);
        }
    };
};

#if defined _WIN32
#define DECLARE_MAIN(MainApplicationClass)                                                                              \
int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)                       \
{                                                                                                                       \
    auto MainApplicationObject = new MainApplicationClass;                                                              \
    MainApplicationObject->Run(MainApplicationObject);                                                                  \
    delete MainApplicationObject;                                                                                       \
    return 0;                                                                                                           \
}
#elif defined _LINUX || defined __APPLE__
#define DECLARE_MAIN(MainApplicationClass)                                                                              \
int main(int argc, const char ** argv)                                                                                  \
{                                                                                                                       \
    auto MainApplicationObject = new MainApplicationClass;                                                              \
    MainApplicationObject->Run(MainApplicationObject);                                                                  \
    delete MainApplicationObject;                                                                                       \
    return 0;                                                                                                           \
}
#else
#error Undefined platform!
#endif

#endif

