
#include <sb7.h>
#include <vmath.h>
#include <object.h>
#include <shader.h>

#include <iostream>
#include <string>
#include <memory>
#include "ArcBall.h"
#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"

#include "CellEnginePDBDataFile.h"
#include "CellEngineCIFDataFile.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace std;

class CellEngineOpenGLVisualiser : public sb7::application
{
private:
    GLuint LineVAO;
    GLuint LineDataBuffer[2];
private:
    GLuint per_fragment_program;
    GLuint per_vertex_program;

    struct
    {
        GLuint color;
        GLuint normals;
    }
    textures{};

    struct uniforms_block
    {
        vmath::mat4 mv_matrix;
        vmath::mat4 view_matrix;
        vmath::mat4 proj_matrix;
        vmath::vec3 color;
    };

    GLuint uniforms_buffer{};

    struct
    {
        GLint diffuse_albedo;
        GLint specular_albedo;
        GLint specular_power;
    }
    uniforms[2]{};

    sb7::object AtomGraphicsObject;

    bool per_vertex;
private:
    Matrix3fT ArcBallPrevRotationMatrix{};
    Matrix3fT ArcBallActualRotationMatrix{};
    std::unique_ptr<ArcBallT> ArcBall;
    Point2fT MousePosition{};
private:
    float LengthUnit = 1;

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
private:
    vector<pair<uint64_t, uint64_t>> AtomBondsToDraw;
    bool FindBondsToDrawFirstTime = true;
private:
    std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;
public:
    CellEngineOpenGLVisualiser() : per_fragment_program(0), per_vertex_program(0), per_vertex(false)
    {
    }
protected:
    void init() override
    {
        static const char title[] = "Cell Engine Visualizer - OpenGL";

        sb7::application::init();

        memcpy(info.title, title, sizeof(title));
    }
protected:
    void InitArcBall();
protected:
    void InitExternalData() override;
protected:
    void InitLineVertexes();
    void DeleteLineVertexes();
    void FindBondsToDraw(const vector<CellEngineAtom>& Atoms);
    void DrawBonds(const vector<CellEngineAtom>& Atoms, const vmath::mat4& ViewMatrix, const vmath::vec3& Center);
    void DrawBond(float x1, float y1, float z1, float x2, float y2, float z2);
protected:
    void load_shaders();
    void startup() override;
    void render(double currentTime) override;
    void onKey(int key, int action) override;
    void onMouseWheel(int pos) override;
    void onMouseButton(int button, int action) override;
    void onMouseMove(int x, int y) override;
    void onResize(int w, int h) override;
protected:
    static inline void DrawCenterPoint(uniforms_block*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix);
    inline vmath::vec3 GetFinalModelPosition(const vmath::vec3& AtomPosition, uniforms_block*  MatrixUniformBlockForVertexShaderPointer, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawOutsideBorder) const;
    inline vmath::vec3 RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, uint64_t& NumberOfAllRenderedAtoms);
protected:
    [[nodiscard]] inline bool CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew) const;
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

void CellEngineOpenGLVisualiser::InitExternalData()
{
    try
    {
        InitializeLoggerManagerParameters();
        LoggersManagerObject.Log(STREAM("START CELL"));

        string FileName;
        if (__argc > 1)
            FileName = __argv[1];
        else
            LoggersManagerObject.Log(STREAM("Lack of file name in program parameters"));

        if (string_utils::check_end_str(FileName, ".pdb") == true)
            CellEngineDataFileObjectPointer = make_unique<CellEnginePDBDataFile>(FileName);
        else
            CellEngineDataFileObjectPointer = make_unique<CellEngineCIFDataFile>(FileName);
    }
    CATCH("reading of data file")
}

void CellEngineOpenGLVisualiser::load_shaders()
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

void CellEngineOpenGLVisualiser::DeleteLineVertexes()
{
    try
    {
        glDeleteVertexArrays(1, &LineVAO);
        glDeleteBuffers(2, LineDataBuffer);
    }
    CATCH("deleting line vertexes")
}

void CellEngineOpenGLVisualiser::InitLineVertexes()
{
    try
    {
        DeleteLineVertexes();

        const float LineVertexes[] = { 0.0, 0.0, 0.0,   1.0, 0.0, 0.0 };
        const float LineNormals[] = { 1.0, 1.0 };

        glGenVertexArrays(1, &LineVAO);
        glBindVertexArray(LineVAO);

        glGenBuffers(2, &LineDataBuffer[0]);

        glBindBuffer(GL_ARRAY_BUFFER, LineDataBuffer[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(LineVertexes), LineVertexes, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, LineDataBuffer[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(LineNormals), LineNormals, GL_STATIC_DRAW);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(1);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    CATCH("initiation of line vertexes")
}

void CellEngineOpenGLVisualiser::DrawBond(float x1, float y1, float z1, float x2, float y2, float z2)
{
    try
    {
        glBindVertexArray(LineVAO);
        const float LineVertexes[] = { x1, y1, z1, x2, y2, z2 };
        glBindBuffer(GL_ARRAY_BUFFER, LineDataBuffer[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(LineVertexes), LineVertexes, GL_STATIC_DRAW);
        glDrawArrays(GL_LINES, 0, 2);
    }
    CATCH("drawing bond")
}

void CellEngineOpenGLVisualiser::FindBondsToDraw(const vector<CellEngineAtom>& Atoms)
{
    try
    {
        if (CellEngineDataFileObjectPointer->DrawBonds == true)
            if (FindBondsToDrawFirstTime == true)
            {
                FindBondsToDrawFirstTime = false;
                AtomBondsToDraw.clear();
                for (uint64_t AtomObjectIndex1 = 0; AtomObjectIndex1 < Atoms.size(); AtomObjectIndex1++)
                    for (uint64_t AtomObjectIndex2 = 0; AtomObjectIndex2 < Atoms.size(); AtomObjectIndex2++)
                        if(AtomObjectIndex1 != AtomObjectIndex2)
                        {
                            const auto& ParticlesCenterObject1 = Atoms[AtomObjectIndex1];
                            const auto& ParticlesCenterObject2 = Atoms[AtomObjectIndex2];

                            float DiffX = ParticlesCenterObject2.X - ParticlesCenterObject1.X;
                            float DiffY = ParticlesCenterObject2.Y - ParticlesCenterObject1.Y;
                            float DiffZ = ParticlesCenterObject2.Z - ParticlesCenterObject1.Z;
                            float VectorLength = sqrt(DiffX * DiffX + DiffY * DiffY + DiffZ * DiffZ);
                            if (VectorLength < 1.5)
                                AtomBondsToDraw.emplace_back(make_pair(AtomObjectIndex1, AtomObjectIndex2));
                        }
            }
    }
    CATCH("finding bonds")
}

void CellEngineOpenGLVisualiser::DrawBonds(const vector<CellEngineAtom>& Atoms, const vmath::mat4& ViewMatrix, const vmath::vec3& Center)
{
    try
    {
        if (CellEngineDataFileObjectPointer->DrawBonds == true)
            for (const auto& AtomBondToDrawObject : AtomBondsToDraw)
            {
                const auto& ParticlesCenterObject1 = Atoms[AtomBondToDrawObject.first];
                const auto& ParticlesCenterObject2 = Atoms[AtomBondToDrawObject.second];

                glBindBufferBase(GL_UNIFORM_BUFFER, 0, uniforms_buffer);
                auto MatrixUniformBlockForVertexShaderPointer = (uniforms_block*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(uniforms_block), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

                MatrixUniformBlockForVertexShaderPointer->view_matrix = ViewMatrix;
                MatrixUniformBlockForVertexShaderPointer->proj_matrix = vmath::perspective(50.0f, (float)info.windowWidth / (float)info.windowHeight, 0.1f, 95000.0f);
                MatrixUniformBlockForVertexShaderPointer->color = vmath::vec3(0.25f, 0.75f, 0.75f);

                MatrixUniformBlockForVertexShaderPointer->mv_matrix = ViewMatrix;

                glUnmapBuffer(GL_UNIFORM_BUFFER);

                DrawBond(ParticlesCenterObject1.X - CameraXPosition - Center.X(), ParticlesCenterObject1.Y - CameraYPosition - Center.Y(), ParticlesCenterObject1.Z - CameraZPosition - Center.Z(), ParticlesCenterObject2.X - CameraXPosition - Center.X(), ParticlesCenterObject2.Y - CameraYPosition - Center.Y(), ParticlesCenterObject2.Z - CameraZPosition - Center.Z());
            }
    }
    CATCH("drawing bonds")
}

void CellEngineOpenGLVisualiser::startup()
{
    try
    {
        load_shaders();

        glGenBuffers(1, &uniforms_buffer);
        glBindBuffer(GL_UNIFORM_BUFFER, uniforms_buffer);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(uniforms_block), nullptr, GL_DYNAMIC_DRAW);

        AtomGraphicsObject.load("..//objects//sphere.sbm");
        InitLineVertexes();

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        InitArcBall();
    }
    CATCH("initiation of data for cell visualization")
}

inline bool CellEngineOpenGLVisualiser::CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew) const
{
    if (CellEngineDataFileObjectPointer->CheckAtomVisibility == true)
    {
        if (ViewZ > CellEngineDataFileObjectPointer->Distance)
            return sqrt((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew)) > CellEngineDataFileObjectPointer->Distance && ZNew > CellEngineDataFileObjectPointer->CutZ;
        else
            return (ZNew > ViewZ - 300 && ZNew < ViewZ && XNew >- 200 && XNew < 200 && YNew > -200 && YNew < 200);
    }
    else
        return false;
}

inline void CellEngineOpenGLVisualiser::DrawCenterPoint(uniforms_block*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix)
{
    try
    {
        ModelMatrix = vmath::translate(0.0f, 0.0f, 0.0f) * vmath::scale(vmath::vec3(0.5, 0.5, 0.5));
        MatrixUniformBlockForVertexShaderPointer->color = vmath::vec3(0.7, 0.2, 0.9);
    }
    CATCH("drawing center point for data for cell visualization")
}

inline vmath::vec3 CellEngineOpenGLVisualiser::GetFinalModelPosition(const vmath::vec3& AtomPosition, uniforms_block*  MatrixUniformBlockForVertexShaderPointer, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawOutsideBorder) const
{
    try
    {
        if (CountNewPosition == true)
        {
            float XNew = MatrixUniformBlockForVertexShaderPointer->mv_matrix[0][0] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[1][0] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[2][0] * (AtomPosition.Z() + CameraZPosition - Center.Z());
            float YNew = MatrixUniformBlockForVertexShaderPointer->mv_matrix[0][1] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[1][1] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[2][1] * (AtomPosition.Z() + CameraZPosition - Center.Z());
            float ZNew = MatrixUniformBlockForVertexShaderPointer->mv_matrix[0][2] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[1][2] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->mv_matrix[2][2] * (AtomPosition.Z() + CameraZPosition - Center.Z());

            if (DrawOutsideBorder == true)
                if (CheckDistanceToDrawDetailsInAtomScale(XNew, YNew, ZNew) == true)
                    MatrixUniformBlockForVertexShaderPointer->color = vmath::vec3(0.7, 0.2, 0.9);

            return vmath::vec3(XNew, YNew, ZNew);
        }
    }
    CATCH("getting final model position for data for cell visualization")

    return vmath::vec3(0, 0, 0);
}

inline vmath::vec3 CellEngineOpenGLVisualiser::RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, uint64_t& NumberOfAllRenderedAtoms)
{
    vmath::vec3 FinalModelPosition;

    try
    {
        NumberOfAllRenderedAtoms++;

        glBindBufferBase(GL_UNIFORM_BUFFER, 0, uniforms_buffer);
        auto MatrixUniformBlockForVertexShaderPointer = (uniforms_block*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(uniforms_block), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

        vmath::vec3 AtomPosition = LengthUnit * AtomObject.Position();
        vmath::mat4 ModelMatrix = vmath::translate(AtomPosition.X() - CameraXPosition - Center.X(), AtomPosition.Y() + CameraYPosition - Center.Y(), AtomPosition.Z() + CameraZPosition - Center.Z()) * vmath::scale(vmath::vec3(CellEngineDataFileObjectPointer->SizeX, CellEngineDataFileObjectPointer->SizeY, CellEngineDataFileObjectPointer->SizeZ));

        MatrixUniformBlockForVertexShaderPointer->view_matrix = ViewMatrix;
        MatrixUniformBlockForVertexShaderPointer->proj_matrix = vmath::perspective(50.0f, (float)info.windowWidth / (float)info.windowHeight, 0.1f, 95000.0f);
        MatrixUniformBlockForVertexShaderPointer->color = AtomObject.Color;

        FinalModelPosition = GetFinalModelPosition(AtomPosition, MatrixUniformBlockForVertexShaderPointer, Center, CountNewPosition, DrawOutsideBorder);
        if (DrawCenter == true)
            DrawCenterPoint(MatrixUniformBlockForVertexShaderPointer, ModelMatrix);

        MatrixUniformBlockForVertexShaderPointer->mv_matrix = ViewMatrix * ModelMatrix;

        glUnmapBuffer(GL_UNIFORM_BUFFER);

        AtomGraphicsObject.render();
    }
    CATCH("rendering object for data for cell visualization")

    return FinalModelPosition;
}

void CellEngineOpenGLVisualiser::render(double currentTime)
{
    try
    {
        static const GLfloat gray[] = { 0.1f, 0.1f, 0.1f, 0.0f };
        static const GLfloat ones[] = { 1.0f };

        glUseProgram(per_vertex ? per_vertex_program : per_fragment_program);
        glViewport(0, 0, info.windowWidth, info.windowHeight);

        glClearBufferfv(GL_COLOR, 0, gray);
        glClearBufferfv(GL_DEPTH, 0, ones);

        vmath::vec3 ViewPositionVector = vmath::vec3(ViewX, ViewY, ViewZ);
        vmath::mat4 ViewMatrix = vmath::lookat(ViewPositionVector, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix;

        vmath::vec3 Center = CellEngineDataFileObjectPointer->GetCenter(CellEngineDataFileObjectPointer->GetParticlesCenters());

        glUniform1f(uniforms[per_vertex ? 1 : 0].specular_power, powf(2.0f, 3.0f));
        glUniform3fv(uniforms[per_vertex ? 1 : 0].specular_albedo, 1, vmath::vec3(1.0f / 9.0f + 1.0f / 9.0f));

        uint64_t NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
        uint64_t NumberOfAllRenderedAtoms = 0;

        FindBondsToDraw(CellEngineDataFileObjectPointer->GetParticlesCenters());
        DrawBonds(CellEngineDataFileObjectPointer->GetParticlesCenters(), ViewMatrix, Center);

        for (auto ParticlesCenterIterator = CellEngineDataFileObjectPointer->GetParticlesCenters().begin(); ParticlesCenterIterator != CellEngineDataFileObjectPointer->GetParticlesCenters().end(); ++ParticlesCenterIterator)
        {
            auto ParticlesCenterObject = *ParticlesCenterIterator;

            vmath::vec3 FinalModelPosition = RenderObject(ParticlesCenterObject, ViewMatrix, Center, true, ParticlesCenterIterator == CellEngineDataFileObjectPointer->GetParticlesCenters().end() - 1, true, NumberOfAllRenderedAtoms);

            if (CellEngineDataFileObjectPointer->ShowDetailsInAtomScale == true)
                if (CheckDistanceToDrawDetailsInAtomScale(FinalModelPosition.X(), FinalModelPosition.Y(), FinalModelPosition.Z()) == true)
                {
                    NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                    for (uint64_t AtomObjectIndex = 0; AtomObjectIndex < CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex].size(); AtomObjectIndex += CellEngineDataFileObjectPointer->LoadOfAtomsStep)
                        RenderObject(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex][AtomObjectIndex], ViewMatrix, Center, false, false, false, NumberOfAllRenderedAtoms);
                }
        }
        CellEngineDataFileObjectPointer->ShowNextStructureFromActiveFilm();

        LoggersManagerObject.Log(STREAM("NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = " << to_string(NumberOfFoundParticlesCenterToBeRenderedInAtomDetails) << " NumberOfAllRenderedAtoms = " << to_string(NumberOfAllRenderedAtoms) << " ViewZ = " << to_string(ViewZ) << " AtomSize = " << to_string(CellEngineDataFileObjectPointer->SizeX) << endl));
    }
    CATCH("rendering cell visualization")
}

void CellEngineOpenGLVisualiser::onKey(int key, int action)
{
    try
    {
        if (action)
            switch (key)
            {
                case 'L': load_shaders(); break;
                case 'K': per_vertex = !per_vertex; break;

                case '1': CameraZPosition += CellEngineDataFileObjectPointer->CameraXMoveStep; break;
                case '2': CameraZPosition -= CellEngineDataFileObjectPointer->CameraXMoveStep; break;
                case '3': CameraXPosition += CellEngineDataFileObjectPointer->CameraYMoveStep; break;
                case '4': CameraXPosition -= CellEngineDataFileObjectPointer->CameraYMoveStep; break;
                case '5': CameraYPosition += CellEngineDataFileObjectPointer->CameraZMoveStep; break;
                case '6': CameraYPosition -= CellEngineDataFileObjectPointer->CameraZMoveStep; break;

                case 'Q': ViewX += CellEngineDataFileObjectPointer->ViewXMoveStep; break;
                case 'W': ViewX -= CellEngineDataFileObjectPointer->ViewXMoveStep; break;
                case 'E': ViewY += CellEngineDataFileObjectPointer->ViewYMoveStep; break;
                case 'R': ViewY -= CellEngineDataFileObjectPointer->ViewYMoveStep; break;
                case 'T': ViewZ += CellEngineDataFileObjectPointer->ViewZMoveStep; break;
                case 'Y': ViewZ -= CellEngineDataFileObjectPointer->ViewZMoveStep; break;

                case 'A': RotationAngle1 += 1; break;
                case 'S': RotationAngle1 -= 1; break;
                case 'D': RotationAngle2 += 1; break;
                case 'F': RotationAngle2 -= 1; break;
                case 'G': RotationAngle3 += 1; break;
                case 'H': RotationAngle3 -= 1; break;

                case 'B': CellEngineDataFileObjectPointer->DrawBonds = !CellEngineDataFileObjectPointer->DrawBonds; break;
                case 'J': CellEngineDataFileObjectPointer->ShowDetailsInAtomScale = !CellEngineDataFileObjectPointer->ShowDetailsInAtomScale; break;

                case 'O': CellEngineDataFileObjectPointer->StartFilmOfStructures(); break;
                case 'P': CellEngineDataFileObjectPointer->StopFilmOfStructures(); break;
                case 'N': CellEngineDataFileObjectPointer->ShowNextStructure(); break;
                case 'M': CellEngineDataFileObjectPointer->ShowPrevStructure(); break;
                case 'Z': CellEngineDataFileObjectPointer->SizeX -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeY -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeZ -= CellEngineDataFileObjectPointer->SizeStep; break;
                case 'X': CellEngineDataFileObjectPointer->SizeX += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeY += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeZ += CellEngineDataFileObjectPointer->SizeStep; break;

                case 'U': CellEngineDataFileObjectPointer->LoadOfAtomsStep == 100 ? CellEngineDataFileObjectPointer->LoadOfAtomsStep = 10 : CellEngineDataFileObjectPointer->LoadOfAtomsStep = 100;  break;
                case 'I': CellEngineDataFileObjectPointer->LoadOfAtomsStep = 1;  break;

                case '9': AtomGraphicsObject.load("..//objects//cube.sbm");  break;
                case '0': AtomGraphicsObject.load("..//objects//sphere.sbm");  break;

                default: break;
            }
    }
    CATCH("executing on key event for cell visualisation")
}

void CellEngineOpenGLVisualiser::onMouseWheel(int pos)
{
    try
    {
        ViewZ += static_cast<float>(pos) * CellEngineDataFileObjectPointer->ViewStep;
    }
    CATCH("executing on mouse wheel event for cell visualisation")
}

void CellEngineOpenGLVisualiser::onMouseButton(int button, int action)
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

void CellEngineOpenGLVisualiser::onMouseMove(int x, int y)
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

void CellEngineOpenGLVisualiser::InitArcBall()
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

void CellEngineOpenGLVisualiser::onResize(int w, int h)
{
    try
    {
        info.windowWidth = w;
        info.windowHeight = h;
        ArcBall->setBounds(static_cast<float>(info.windowWidth), static_cast<float>(info.windowHeight));
    }
    CATCH("executing window resize event - setting bounds of arc ball counting data for cell visualisation")
}

DECLARE_MAIN(CellEngineOpenGLVisualiser)