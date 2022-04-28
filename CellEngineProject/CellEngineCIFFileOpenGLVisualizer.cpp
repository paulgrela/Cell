
#include <iostream>
#include <numeric>

#include <sb7.h>
#include <vmath.h>

#include <object.h>
#include <sb7ktx.h>
#include <shader.h>

#include <string>
#include <memory>
#include "ArcBall.h"
#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"
#include "AdditionalFunctions.h"

#include "CellEnginePDBDataFile.h"

using namespace std;

class phonglighting_app : public sb7::application
{
public:
    phonglighting_app() : per_fragment_program(0), per_vertex_program(0), per_vertex(false)
    {
    }

protected:
    void init() override
    {
        static const char title[] = "Cell Engine Visualizer - OpenGL";

        sb7::application::init();

        memcpy(info.title, title, sizeof(title));
    }

    void InitArcBall();

    void startup() override;
    void load_shaders();
    void render(double currentTime) override;
    void onKey(int key, int action) override;
    void onMouseWheel(int pos) override;
    void onMouseButton(int button, int action) override;
    void onMouseMove(int x, int y) override;
    void onResize(int w, int h) override;


    GLuint per_fragment_program;
    GLuint per_vertex_program;

    struct
    {
        GLuint color;
        GLuint normals;
    } 
    textures;

    struct uniforms_block
    {
        vmath::mat4 mv_matrix;
        vmath::mat4 view_matrix;
        vmath::mat4 proj_matrix;
        vmath::vec3 color;
    };

    GLuint uniforms_buffer;

    struct
    {
        GLint diffuse_albedo;
        GLint specular_albedo;
        GLint specular_power;
    } 
    uniforms[2];

    sb7::object object;

    bool per_vertex;
private:
    Matrix3fT ArcBallPrevRotationMatrix;
    Matrix3fT ArcBallActualRotationMatrix;
    std::unique_ptr<ArcBallT> ArcBall;
    Point2fT MousePosition;
private:
    std::unique_ptr<PDBDataFile> PDBDataFileObjectPointer;
private:
    float CameraXPosition = 0.0;
    float CameraYPosition = 0.0;
    float CameraZPosition = 0.0;

    float ViewX = 0.0;
    float ViewY = 0.0;
    float ViewZ = 50.0;

    vmath::mat4 RotationMatrix;

    float RotationAngle1 = 0.0;
    float RotationAngle2 = 0.0;
    float RotationAngle3 = 0.0;
private:
    uint64_t PressedRightMouseButton = 0;
};

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

void phonglighting_app::load_shaders()
{
    try
    {
        GLuint vs;
        GLuint fs;

        vs = sb7::shader::load("..\\shaders\\per-fragment-phong.vs.glsl", GL_VERTEX_SHADER);
        fs = sb7::shader::load("..\\shaders\\per-fragment-phong.fs.glsl", GL_FRAGMENT_SHADER);

        if (per_fragment_program)
            glDeleteProgram(per_fragment_program);

        per_fragment_program = glCreateProgram();
        glAttachShader(per_fragment_program, vs);
        glAttachShader(per_fragment_program, fs);
        glLinkProgram(per_fragment_program);

        uniforms[0].diffuse_albedo = glGetUniformLocation(per_fragment_program, "diffuse_albedo");
        uniforms[0].specular_albedo = glGetUniformLocation(per_fragment_program, "specular_albedo");
        uniforms[0].specular_power = glGetUniformLocation(per_fragment_program, "specular_power");

        vs = sb7::shader::load("..\\shaders\\per-vertex-phong.vs.glsl", GL_VERTEX_SHADER);
        fs = sb7::shader::load("..\\shaders\\per-vertex-phong.fs.glsl", GL_FRAGMENT_SHADER);

        if (per_vertex_program)
            glDeleteProgram(per_vertex_program);

        per_vertex_program = glCreateProgram();
        glAttachShader(per_vertex_program, vs);
        glAttachShader(per_vertex_program, fs);
        glLinkProgram(per_vertex_program);

        uniforms[1].diffuse_albedo = glGetUniformLocation(per_vertex_program, "diffuse_albedo");
        uniforms[1].specular_albedo = glGetUniformLocation(per_vertex_program, "specular_albedo");
        uniforms[1].specular_power = glGetUniformLocation(per_vertex_program, "specular_power");
    }
    CATCH("loading shaders for cell visualization")
}

void phonglighting_app::startup()
{
    try
    {
        load_shaders();

        glGenBuffers(1, &uniforms_buffer);
        glBindBuffer(GL_UNIFORM_BUFFER, uniforms_buffer);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(uniforms_block), nullptr, GL_DYNAMIC_DRAW);

        object.load("sphere.sbm");

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        InitArcBall();

        InitializeLoggerManagerParameters();

        string FileName = R"(c:\Projects\Programs\C++\Cell\PDBFiles\struct0\117E.pdb)";

        PDBDataFileObjectPointer = make_unique<PDBDataFile>(FileName);
    }
    CATCH("initation of data for cell visualization")
}

void phonglighting_app::render(double currentTime)
{
    try
    {
        const float LengthUnit = 0.3;
        FloatVectorType MassCenter = PDBDataFileObjectPointer->MassCenter();

        static const GLfloat zeros[] = { 0.0f, 0.0f, 0.0f, 0.0f };
        static const GLfloat gray[] = { 0.1f, 0.1f, 0.1f, 0.0f };
        static const GLfloat ones[] = { 1.0f };

        glUseProgram(per_vertex ? per_vertex_program : per_fragment_program);
        glViewport(0, 0, info.windowWidth, info.windowHeight);

        glClearBufferfv(GL_COLOR, 0, gray);
        glClearBufferfv(GL_DEPTH, 0, ones);

        vmath::vec3 view_position = vmath::vec3(ViewX, ViewY, ViewZ);
        vmath::mat4 view_matrix = vmath::lookat(view_position, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix;

        int j = 1;
        for (const Element& ElementObject : PDBDataFileObjectPointer->GetElements())
        {
            glBindBufferBase(GL_UNIFORM_BUFFER, 0, uniforms_buffer);
            auto block = (uniforms_block*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(uniforms_block), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

            FloatVectorType ElementPosition = LengthUnit * ElementObject.Position();

            vmath::mat4 model_matrix = vmath::translate(ElementPosition.X * 3.0f - CameraXPosition - MassCenter.X, ElementPosition.Y * 3.0f + CameraYPosition - MassCenter.Y, ElementPosition.Z * 3.0f + CameraZPosition - MassCenter.Z);

            block->mv_matrix = view_matrix * model_matrix;
            block->view_matrix = view_matrix;
            block->proj_matrix = vmath::perspective(50.0f, (float)info.windowWidth / (float)info.windowHeight, 0.1f, 2000.0f);

            switch(ElementObject.Name[0])
            {
                case 'C': block->color = vmath::vec3(0.25f, 0.75f, 0.75f); break;
                case 'O': block->color = vmath::vec3(1.00f, 0.00f, 0.00f); break;
                case 'H': block->color = vmath::vec3(1.00f, 1.00f, 1.00f); break;
                case 'N': block->color = vmath::vec3(0.00f, 0.00f, 1.00f); break;
                case 'P': block->color = vmath::vec3(0.50f, 0.50f, 0.20f); break;
                default: block->color = vmath::vec3(0.50f, 0.50f, 0.50f); break;
            }
            //block->color = vmath::vec3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);;

            glUnmapBuffer(GL_UNIFORM_BUFFER);

            glUniform1f(uniforms[per_vertex ? 1 : 0].specular_power, powf(2.0f, (float)j + 2.0f));
            glUniform3fv(uniforms[per_vertex ? 1 : 0].specular_albedo, 1, vmath::vec3((float)j / 9.0f + 1.0f / 9.0f));
            j++;

            object.render();
        }



        glBindBufferBase(GL_UNIFORM_BUFFER, 0, uniforms_buffer);
        auto block = (uniforms_block*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(uniforms_block), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

        vmath::mat4 model_matrix = vmath::translate(0.0f, 0.0f, 0.0f);

        block->mv_matrix = model_matrix * view_matrix;
        block->view_matrix = view_matrix;
        block->proj_matrix = vmath::perspective(50.0f, (float)info.windowWidth / (float)info.windowHeight, 0.1f, 2000.0f);

        block->color = vmath::vec3(0.7, 0.2, 0.9);

        glUnmapBuffer(GL_UNIFORM_BUFFER);

        glUniform1f(uniforms[per_vertex ? 1 : 0].specular_power, powf(2.0f, (float)j + 2.0f));
        glUniform3fv(uniforms[per_vertex ? 1 : 0].specular_albedo, 1, vmath::vec3((float)j / 9.0f + 1.0f / 9.0f));

        object.render();

    }
    CATCH("rendering cell visualization")
}

void phonglighting_app::onKey(int key, int action)
{
    try
    {
        if (action)
            switch (key)
            {
                case 'M': load_shaders(); break;
                case 'V': per_vertex = !per_vertex; break;

                case '1': CameraZPosition += 1; break;
                case '2': CameraZPosition -= 1; break;
                case '3': CameraXPosition += 1; break;
                case '4': CameraXPosition -= 1; break;
                case '5': CameraYPosition += 1; break;
                case '6': CameraYPosition -= 1; break;

                case 'Q': ViewX += 1; break;
                case 'W': ViewX -= 1; break;
                case 'E': ViewY += 1; break;
                case 'R': ViewY -= 1; break;
                case 'T': ViewZ += 1; break;
                case 'Y': ViewZ -= 1; break;

                case 'A': RotationAngle1 += 1; break;
                case 'S': RotationAngle1 -= 1; break;
                case 'D': RotationAngle2 += 1; break;
                case 'F': RotationAngle2 -= 1; break;
                case 'G': RotationAngle3 += 1; break;
                case 'H': RotationAngle3 -= 1; break;

                default: break;
            }
    }
    CATCH("executing on key event for cell visualisation")
}

void phonglighting_app::onMouseWheel(int pos)
{
    try
    {
        if (pos > 0)
            ViewZ += 1;
        else
            ViewZ -= 1;
    }
    CATCH("executing on mouse wheel event")
}

void phonglighting_app::onMouseButton(int button, int action)
{
    try
    {
        if (button == 0)
        {
            PressedRightMouseButton++;
            if (PressedRightMouseButton == 1)
            {
                ArcBallPrevRotationMatrix = ArcBallActualRotationMatrix;
                ArcBall->click(&MousePosition);
            }
            else
            if (PressedRightMouseButton == 2)
                PressedRightMouseButton = 0;
        }
    }
    CATCH("executing on mouse button event for cell visualisation")
}

void phonglighting_app::onMouseMove(int x, int y)
{
    try
    {
        MousePosition.s.X = x;
        MousePosition.s.Y = y;

        if (PressedRightMouseButton == true)
        {
            Quat4fT ThisQuat;
            ArcBall->drag(&MousePosition, &ThisQuat);
            Matrix3fSetRotationFromQuat4f(&ArcBallActualRotationMatrix, &ThisQuat);
            Matrix3fMulMatrix3f(&ArcBallActualRotationMatrix, &ArcBallPrevRotationMatrix);

            RotationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);

            RotationMatrix.data[0][0] = ArcBallActualRotationMatrix.s.M00;
            RotationMatrix.data[1][0] = ArcBallActualRotationMatrix.s.M01;
            RotationMatrix.data[2][0] = ArcBallActualRotationMatrix.s.M02;
            RotationMatrix.data[0][1] = ArcBallActualRotationMatrix.s.M10;
            RotationMatrix.data[1][1] = ArcBallActualRotationMatrix.s.M11;
            RotationMatrix.data[2][1] = ArcBallActualRotationMatrix.s.M12;
            RotationMatrix.data[0][2] = ArcBallActualRotationMatrix.s.M20;
            RotationMatrix.data[1][2] = ArcBallActualRotationMatrix.s.M21;
            RotationMatrix.data[2][2] = ArcBallActualRotationMatrix.s.M22;
        }
    }
    CATCH("executing on mouse event for cell visualisation")
}

void phonglighting_app::InitArcBall()
{
    try
    {
        Matrix3fSetIdentity(&ArcBallPrevRotationMatrix);
        Matrix3fSetIdentity(&ArcBallActualRotationMatrix);

        ArcBall = make_unique<ArcBallT>(ArcBallT(640.0f, 480.0f));
        ArcBall->setBounds(static_cast<float>(info.windowWidth), static_cast<float>(info.windowHeight));

        RotationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);
    }
    CATCH("initiation of arc ball counting data for cell visualisation")
}

void phonglighting_app::onResize(int w, int h)
{
    try
    {
        info.windowWidth = w;
        info.windowHeight = h;
        ArcBall->setBounds(static_cast<float>(info.windowWidth), static_cast<float>(info.windowHeight));
    }
    CATCH("executing window resize event - setting bounds of arc ball counting data for cell visualisation")
}

DECLARE_MAIN(phonglighting_app)