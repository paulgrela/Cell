
#include <omp.h>
#include <sb7.h>
#include <sb7color.h>
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

class CellEngineOpenGLVisualiser : public sb7::OpenGLApplication
{
private:
    GLuint LineVAO;
    GLuint LineDataBuffer[2];
private:
    GLuint ShaderProgramPhong = 0;
    GLuint ShaderProgramSimple = 0;
private:
    struct UniformsBlock
    {
        vmath::mat4 MoveMatrix;
        vmath::mat4 ProjMatrix;
        vmath::vec3 Color;
    };
private:
    GLuint UniformsBuffer{};
private:
    struct
    {
        GLint DiffuseAlbedo;
        GLint SpecularAlbedo;
        GLint SpecularPower;
    }
    Uniforms{};
private:
    sb7::GraphicObject AtomGraphicsObject;
private:
    Matrix3fT ArcBallPrevRotationMatrix{};
    Matrix3fT ArcBallActualRotationMatrix{};
    std::unique_ptr<ArcBallT> ArcBall;
    Point2fT MousePosition{};
private:
    float LengthUnit = 1;
private:
    float CameraXPosition = 0.0;
    float CameraYPosition = 0.0;
    float CameraZPosition = 0.0;
private:
    float ViewX = 0.0;
    float ViewY = 0.0;
    float ViewZ = 50.0;
private:
    vmath::mat4 RotationMatrix;
private:
    float RotationAngle1 = 0.0;
    float RotationAngle2 = 0.0;
    float RotationAngle3 = 0.0;
private:
    UnsignedIntType PressedRightMouseButton = 0;
private:
    bool RenderObjects = true;
private:
    vector<pair<UnsignedIntType, UnsignedIntType>> BondsBetweenParticlesCentersToDraw;
    vector<vector<pair<UnsignedIntType, UnsignedIntType>>> BondsBetweenAtomsToDraw;
private:
    std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;
public:
    CellEngineOpenGLVisualiser() = default;
protected:
    void Init() override
    {
        static const char title[] = "Cell Engine Visualizer";

        sb7::OpenGLApplication::Init();

        memcpy(Info.Title, title, sizeof(title));
    }
protected:
    void InitArcBall();
protected:
    void InitExternalData() override;
protected:
    void InitLineVertexes();
    void DeleteLineVertexes();
    static void FindBondsToDraw(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw);
    void DrawBonds(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix, const vmath::vec3& Center);
    void DrawBond(float x1, float y1, float z1, float x2, float y2, float z2);
protected:
    bool CheckVisibilityOfParticles(UnsignedIntType EntityId);
    void SetVisibilityOfAllParticles(bool VisibleParam);
    void SetVisibilityOfParticlesExcept(UnsignedIntType EntityId, bool VisibleParam);
protected:
    void LoadShadersPhong();
    void LoadShadersSimple();
    static void LoadShaders(const char* VertexShaderFileName, const char* FragmentShaderFileName, GLuint& ShaderProgram);
protected:
    void StartUp() override;
    void Render(double CurrentTime) override;
    void OnKey(int Key, int Action) override;
    void OnMouseWheel(int Pos) override;
    void OnMouseButton(int Button, int Action) override;
    void OnMouseMove(int X, int Y) override;
    void OnResize(int Width, int Height) override;
protected:
    inline vmath::vec3 GetColor(const CellEngineAtom& AtomObject);
    static inline void DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix);
    inline vmath::vec3 GetFinalModelPosition(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawOutsideBorder) const;
    inline vmath::vec3 RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedIntType& NumberOfAllRenderedAtoms);
    inline vmath::vec3 CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional);
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

void CellEngineOpenGLVisualiser::StartUp()
{
    try
    {
        LoadShadersSimple();
        LoadShadersPhong();

        glGenBuffers(1, &UniformsBuffer);
        glBindBuffer(GL_UNIFORM_BUFFER, UniformsBuffer);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformsBlock), nullptr, GL_DYNAMIC_DRAW);

        AtomGraphicsObject.Load("..//objects//sphere.sbm");
        InitLineVertexes();

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        InitArcBall();
    }
    CATCH("initiation of data for cell visualization")
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

        BondsBetweenAtomsToDraw.resize(CellEngineDataFileObjectPointer->GetParticlesCenters().size());
    }
    CATCH("reading of data file")
}

void CellEngineOpenGLVisualiser::LoadShaders(const char* VertexShaderFileName, const char* FragmentShaderFileName, GLuint& ShaderProgram)
{
    try
    {
        GLuint VertexShader = sb7::shader::Load(VertexShaderFileName, GL_VERTEX_SHADER);
        GLuint FragmentShader = sb7::shader::Load(FragmentShaderFileName, GL_FRAGMENT_SHADER);

        if (ShaderProgram)
            glDeleteProgram(ShaderProgram);

        ShaderProgram = glCreateProgram();
        glAttachShader(ShaderProgram, VertexShader);
        glAttachShader(ShaderProgram, FragmentShader);
        glLinkProgram(ShaderProgram);

        glDeleteShader(VertexShader);
        glDeleteShader(FragmentShader);
    }
    CATCH_AND_THROW("loading phong shaders for cell visualization")
}

void CellEngineOpenGLVisualiser::LoadShadersPhong()
{
    try
    {
        LoadShaders("..\\shaders\\per-fragment-phong.vs.glsl", "..\\shaders\\per-fragment-phong.fs.glsl", ShaderProgramPhong);

        Uniforms.DiffuseAlbedo = glGetUniformLocation(ShaderProgramPhong, "diffuse_albedo");
        Uniforms.SpecularAlbedo = glGetUniformLocation(ShaderProgramPhong, "specular_albedo");
        Uniforms.SpecularPower = glGetUniformLocation(ShaderProgramPhong, "specular_power");
    }
    CATCH("loading phong shaders for cell visualization")
}

void CellEngineOpenGLVisualiser::LoadShadersSimple()
{
    try
    {
        LoadShaders("..\\shaders\\per-fragment-simple.vs.glsl", "..\\shaders\\per-fragment-simple.fs.glsl", ShaderProgramSimple);
    }
    CATCH("loading simple shaders for cell visualization")
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

void CellEngineOpenGLVisualiser::FindBondsToDraw(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw)
{
    try
    {
        //..#pragma omp parallel for shared(Atoms) private(BondsToDraw)
        for (UnsignedIntType AtomObjectIndex1 = 0; AtomObjectIndex1 < Atoms.size(); AtomObjectIndex1++)
            for (UnsignedIntType AtomObjectIndex2 = 0; AtomObjectIndex2 < Atoms.size(); AtomObjectIndex2++)
                if (AtomObjectIndex1 != AtomObjectIndex2)
                {
                    const auto& ParticlesCenterObject1 = Atoms[AtomObjectIndex1];
                    const auto& ParticlesCenterObject2 = Atoms[AtomObjectIndex2];

                    float DiffX = ParticlesCenterObject2.X - ParticlesCenterObject1.X;
                    float DiffY = ParticlesCenterObject2.Y - ParticlesCenterObject1.Y;
                    float DiffZ = ParticlesCenterObject2.Z - ParticlesCenterObject1.Z;
                    float VectorLength = sqrt(DiffX * DiffX + DiffY * DiffY + DiffZ * DiffZ);
                    if (VectorLength < 1.5)
                        BondsToDraw.emplace_back(make_pair(AtomObjectIndex1, AtomObjectIndex2));
                }
    }
    CATCH("finding bonds")
}

void CellEngineOpenGLVisualiser::DrawBonds(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix, const vmath::vec3& Center)
{
    try
    {
        if (DrawBonds == true)
        {
            if (BondsToDraw.empty() == true)
                FindBondsToDraw(Atoms, BondsToDraw);

            for (const auto& BondToDrawObject : BondsToDraw)
            {
                const auto& AtomObject1 = Atoms[BondToDrawObject.first];
                const auto& AtomObject2 = Atoms[BondToDrawObject.second];

                CreateUniformBlockForVertexShader(vmath::vec3(0.0, 0.0, 0.0), FromVec4ToVec3(sb7::color::DarkCyan), ViewMatrix, vmath::translate(0.0f, 0.0f, 0.0f), Center, false, false, false, false);

                DrawBond(AtomObject1.X - CameraXPosition - Center.X(), AtomObject1.Y - CameraYPosition - Center.Y(), AtomObject1.Z - CameraZPosition - Center.Z(), AtomObject2.X - CameraXPosition - Center.X(), AtomObject2.Y - CameraYPosition - Center.Y(), AtomObject2.Z - CameraZPosition - Center.Z());
            }
        }
    }
    CATCH("drawing bonds")
}

inline bool CellEngineOpenGLVisualiser::CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew) const
{
    if (CellEngineDataFileObjectPointer->CheckAtomVisibility == true)
    {
        if (ViewZ > CellEngineDataFileObjectPointer->Distance)
            return sqrt((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew)) > CellEngineDataFileObjectPointer->Distance && ZNew > CellEngineDataFileObjectPointer->CutZ;
            //return ((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew) > CellEngineDataFileObjectPointer->Distance * CellEngineDataFileObjectPointer->Distance) && ZNew > CellEngineDataFileObjectPointer->CutZ;
        else
            return (ZNew > ViewZ - 300 && ZNew < ViewZ + 100 && XNew >- 200 && XNew < 200 && YNew > -200 && YNew < 200);
            //return (XNew >- 10 && XNew < 10 && YNew > -10 && YNew < 10 && ZNew > -10);
    }
    else
        return false;
}

inline void CellEngineOpenGLVisualiser::DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix)
{
    try
    {
        ModelMatrix = vmath::translate(0.0f, 0.0f, 0.0f) * vmath::scale(vmath::vec3(0.5, 0.5, 0.5));
        MatrixUniformBlockForVertexShaderPointer->Color = FromVec4ToVec3(sb7::color::Purple);
    }
    CATCH("drawing center point for data for cell visualization")
}

inline vmath::vec3 CellEngineOpenGLVisualiser::GetFinalModelPosition(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawOutsideBorder) const
{
    try
    {
        if (CountNewPosition == true)
        {
            float XNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][0] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][0] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][0] * (AtomPosition.Z() + CameraZPosition - Center.Z());
            float YNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][1] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][1] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][1] * (AtomPosition.Z() + CameraZPosition - Center.Z());
            float ZNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][2] * (AtomPosition.X() + CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][2] * (AtomPosition.Y() + CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][2] * (AtomPosition.Z() + CameraZPosition - Center.Z());

            if (DrawOutsideBorder == true)
                if (CheckDistanceToDrawDetailsInAtomScale(XNew, YNew, ZNew) == true)
                    MatrixUniformBlockForVertexShaderPointer->Color = FromVec4ToVec3(sb7::color::Purple);

            return vmath::vec3(XNew, YNew, ZNew);
        }
    }
    CATCH("getting final model position for data for cell visualization")

    return vmath::vec3(0, 0, 0);
}

inline vmath::vec3 CellEngineOpenGLVisualiser::CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional)
{
    vmath::vec3 FinalModelPosition;

    try
    {
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, UniformsBuffer);
        auto MatrixUniformBlockForVertexShaderPointer = (UniformsBlock*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(UniformsBlock), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

        MatrixUniformBlockForVertexShaderPointer->ProjMatrix = vmath::perspective(50.0f, (float)Info.WindowWidth / (float)Info.WindowHeight, 0.1f, 95000.0f);
        MatrixUniformBlockForVertexShaderPointer->Color = Color;

        if (DrawAdditional == true)
        {
            FinalModelPosition = GetFinalModelPosition(Position, MatrixUniformBlockForVertexShaderPointer, Center, CountNewPosition, DrawOutsideBorder);
            if (DrawCenter == true)
                DrawCenterPoint(MatrixUniformBlockForVertexShaderPointer, ModelMatrix);
        }

        MatrixUniformBlockForVertexShaderPointer->MoveMatrix = ViewMatrix * ModelMatrix;

        glUnmapBuffer(GL_UNIFORM_BUFFER);
    }
    CATCH("rendering object for data for cell visualization")

    return FinalModelPosition;
}

inline vmath::vec3 CellEngineOpenGLVisualiser::GetColor(const CellEngineAtom& AtomObject)
{
    vmath::vec3 FinalColor;

    try
    {
        if (CellEngineDataFileObjectPointer->DrawRandomColorForEveryParticle == true)
            FinalColor = AtomObject.RandomParticleColor;
        else
            FinalColor = (CellEngineDataFileObjectPointer->DrawColorForEveryAtom == false ? AtomObject.ParticleColor : AtomObject.AtomColor);
    }
    CATCH("getting color")

    return FinalColor;
}

inline vmath::vec3 CellEngineOpenGLVisualiser::RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const vmath::vec3& Center, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedIntType& NumberOfAllRenderedAtoms)
{
    vmath::vec3 FinalModelPosition;

    try
    {
        NumberOfAllRenderedAtoms++;

        vmath::vec3 AtomPosition = LengthUnit * AtomObject.Position();
        vmath::mat4 ModelMatrix = vmath::translate(AtomPosition.X() - CameraXPosition - Center.X(), AtomPosition.Y() + CameraYPosition - Center.Y(), AtomPosition.Z() + CameraZPosition - Center.Z()) * vmath::scale(vmath::vec3(CellEngineDataFileObjectPointer->SizeX, CellEngineDataFileObjectPointer->SizeY, CellEngineDataFileObjectPointer->SizeZ));

        FinalModelPosition = CreateUniformBlockForVertexShader(AtomPosition, GetColor(AtomObject), ViewMatrix, ModelMatrix, Center, CountNewPosition, DrawCenter, DrawOutsideBorder, true);

        if (RenderObjects == true)
            AtomGraphicsObject.Render();
    }
    CATCH("rendering object for data for cell visualization")

    return FinalModelPosition;
}

void CellEngineOpenGLVisualiser::Render(double CurrentTime)
{
    try
    {
        static const GLfloat gray[] = { 0.1f, 0.1f, 0.1f, 0.0f };
        static const GLfloat ones[] = { 1.0f };

        glDisable(GL_SCISSOR_TEST);

        glViewport(0, 0, Info.WindowWidth, Info.WindowHeight);

        glClearBufferfv(GL_COLOR, 0, gray);
        glClearBufferfv(GL_DEPTH, 0, ones);

        glClearColor(0, 0, 0, 0);
        glClearStencil(0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glEnable(GL_STENCIL_TEST);
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

        vmath::vec3 ViewPositionVector = vmath::vec3(ViewX, ViewY, ViewZ);
        vmath::mat4 ViewMatrix = vmath::lookat(ViewPositionVector, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix;

        vmath::vec3 Center = CellEngineDataFileObjectPointer->GetCenter(CellEngineDataFileObjectPointer->GetParticlesCenters());

        glUniform1f(Uniforms.SpecularPower, powf(2.0f, 3.0f));
        glUniform3fv(Uniforms.SpecularAlbedo, 1, vmath::vec3(1.0f / 9.0f + 1.0f / 9.0f));

        UnsignedIntType NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
        UnsignedIntType NumberOfAllRenderedAtoms = 0;

        glUseProgram(ShaderProgramSimple);
        DrawBonds(CellEngineDataFileObjectPointer->GetParticlesCenters(), BondsBetweenParticlesCentersToDraw, CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters, ViewMatrix, Center);

        glUseProgram(ShaderProgramPhong);
        GLuint PartOfStencilBufferIndex[3];
        for (uint64_t LoopCounter = 0; LoopCounter < CellEngineDataFileObjectPointer->NumberOfStencilBufferLoop; LoopCounter++)
        {
            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
            NumberOfAllRenderedAtoms = 0;

            if (LoopCounter > 0)
            {
                glEnable(GL_SCISSOR_TEST);
                glScissor(MousePosition.s.X, (float)Info.WindowHeight - MousePosition.s.Y - 1, 1, 1);
            }
            else
                glDisable(GL_SCISSOR_TEST);

            for (auto ParticlesCenterIterator = CellEngineDataFileObjectPointer->GetParticlesCenters().begin(); ParticlesCenterIterator != CellEngineDataFileObjectPointer->GetParticlesCenters().end(); ++ParticlesCenterIterator)
            {
                uint8_t ToInsert = (NumberOfAllRenderedAtoms) >> (8 * LoopCounter);
                glStencilFunc(GL_ALWAYS, ToInsert, -1);

                auto ParticlesCenterObject = *ParticlesCenterIterator;

                vmath::vec3 FinalModelPosition = RenderObject(ParticlesCenterObject, ViewMatrix, Center, true, ParticlesCenterIterator == CellEngineDataFileObjectPointer->GetParticlesCenters().end() - 1, true, NumberOfAllRenderedAtoms);

                if (CellEngineDataFileObjectPointer->ShowDetailsInAtomScale == true)
                    if (CheckDistanceToDrawDetailsInAtomScale(FinalModelPosition.X(), FinalModelPosition.Y(), FinalModelPosition.Z()) == true)
                        if (CheckVisibilityOfParticles(ParticlesCenterObject.EntityId) == true)
                        {
                            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                            glUseProgram(ShaderProgramSimple);
                            DrawBonds(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex], BondsBetweenAtomsToDraw[ParticlesCenterObject.AtomIndex], CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms, ViewMatrix, Center);

                            glUseProgram(ShaderProgramPhong);
                            for (UnsignedIntType AtomObjectIndex = 0; AtomObjectIndex < CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex].size(); AtomObjectIndex += CellEngineDataFileObjectPointer->LoadOfAtomsStep)
                                RenderObject(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex][AtomObjectIndex], ViewMatrix, Center, false, false, false, NumberOfAllRenderedAtoms);
                        }
            }
            CellEngineDataFileObjectPointer->ShowNextStructureFromActiveFilm();

            GLuint StencilIndex;
            glReadPixels(MousePosition.s.X, (float)Info.WindowHeight - MousePosition.s.Y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &StencilIndex);
            PartOfStencilBufferIndex[LoopCounter] = StencilIndex;
        }

        uint64_t ParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

        if (ParticleCenterIndex > CellEngineDataFileObjectPointer->GetParticlesCenters().size())
            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + to_string(ParticleCenterIndex) + " MAXIMAL NUMBER OF OBJECTS = " + to_string(CellEngineDataFileObjectPointer->GetParticlesCenters().size()));

        auto FoundParticleObjectMarked = CellEngineDataFileObjectPointer->GetParticlesCenters()[ParticleCenterIndex];

        glDisable(GL_SCISSOR_TEST);




        //RenderObject(*FoundParticleObjectMarked, ViewMatrix, Center, false, true, false, NumberOfAllRenderedAtoms);

        vmath::vec3 AtomPosition = LengthUnit * FoundParticleObjectMarked.Position();

        vmath::mat4 ModelMatrix = vmath::translate(AtomPosition.X() - CameraXPosition - Center.X(), AtomPosition.Y() + CameraYPosition - Center.Y(), AtomPosition.Z() + CameraZPosition - Center.Z()) * vmath::scale(vmath::vec3(CellEngineDataFileObjectPointer->SizeX, CellEngineDataFileObjectPointer->SizeY, CellEngineDataFileObjectPointer->SizeZ));

        glBindBufferBase(GL_UNIFORM_BUFFER, 0, UniformsBuffer);
        auto MatrixUniformBlockForVertexShaderPointer = (UniformsBlock*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(UniformsBlock), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

        MatrixUniformBlockForVertexShaderPointer->ProjMatrix = vmath::perspective(50.0f, (float)Info.WindowWidth / (float)Info.WindowHeight, 0.1f, 95000.0f);
        MatrixUniformBlockForVertexShaderPointer->Color = FromVec4ToVec3(sb7::color::Yellow);

        MatrixUniformBlockForVertexShaderPointer->MoveMatrix = ViewMatrix * ModelMatrix;

        glUnmapBuffer(GL_UNIFORM_BUFFER);

        AtomGraphicsObject.Render();
        LoggersManagerObject.Log(STREAM("FOUND CHOSEN OBJECT INDEX = " << ParticleCenterIndex << " POSITION[ X = " << AtomPosition.X() << " Y = " << AtomPosition.Y() << " Z = " << AtomPosition.Z() << " ] " << endl));


        LoggersManagerObject.Log(STREAM("NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = " << to_string(NumberOfFoundParticlesCenterToBeRenderedInAtomDetails) << " NumberOfAllRenderedAtoms = " << to_string(NumberOfAllRenderedAtoms) << " ViewZ = " << to_string(ViewZ) << " AtomSize = " << to_string(CellEngineDataFileObjectPointer->SizeX) << endl));
    }
    CATCH("rendering cell visualization")
}

void CellEngineOpenGLVisualiser::SetVisibilityOfAllParticles(bool VisibleParam)
{
    try
    {
        for (auto& ParticleKindObject : CellEngineDataFileObjectPointer->ParticlesKinds)
            ParticleKindObject.second.Visible = VisibleParam;
    }
    CATCH("setting visibility of all particles")
}

void CellEngineOpenGLVisualiser::SetVisibilityOfParticlesExcept(UnsignedIntType EntityId, bool VisibleParam)
{
    try
    {
        SetVisibilityOfAllParticles(VisibleParam);
        CellEngineDataFileObjectPointer->ParticlesKinds.find(EntityId)->second.Visible = !VisibleParam;
    }
    CATCH("setting visibility of particles except")
}

bool CellEngineOpenGLVisualiser::CheckVisibilityOfParticles(UnsignedIntType EntityId)
{
    bool Visible = false;

    try
    {
        Visible = CellEngineDataFileObjectPointer->ParticlesKinds.find(EntityId)->second.Visible;
    }
    CATCH("checking visibility of particles")

    return Visible;
}

void CellEngineOpenGLVisualiser::OnKey(int Key, int Action)
{
    try
    {
        if (Action)
            switch (Key)
            {
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

                case 'B': CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters = !CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters; break;
                case 'C': CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms = !CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms; break;

                case 'J': CellEngineDataFileObjectPointer->ShowDetailsInAtomScale = !CellEngineDataFileObjectPointer->ShowDetailsInAtomScale; break;

                case 'O': CellEngineDataFileObjectPointer->StartFilmOfStructures(); break;
                case 'P': CellEngineDataFileObjectPointer->StopFilmOfStructures(); break;
                case 'N': CellEngineDataFileObjectPointer->ShowNextStructure(); break;
                case 'M': CellEngineDataFileObjectPointer->ShowPrevStructure(); break;
                case 'Z': CellEngineDataFileObjectPointer->SizeX -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeY -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeZ -= CellEngineDataFileObjectPointer->SizeStep; break;
                case 'X': CellEngineDataFileObjectPointer->SizeX += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeY += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeZ += CellEngineDataFileObjectPointer->SizeStep; break;

                case 'U': CellEngineDataFileObjectPointer->LoadOfAtomsStep == 100 ? CellEngineDataFileObjectPointer->LoadOfAtomsStep = 10 : CellEngineDataFileObjectPointer->LoadOfAtomsStep = 100;  break;
                case 'I': CellEngineDataFileObjectPointer->LoadOfAtomsStep = 1;  break;

                case '9': AtomGraphicsObject.Load("..//objects//cube.sbm");  break;
                case '0': AtomGraphicsObject.Load("..//objects//sphere.sbm");  break;

                case GLFW_KEY_F7: RenderObjects = !RenderObjects; break;
                case GLFW_KEY_F8: CellEngineDataFileObjectPointer->NumberOfStencilBufferLoop == 1 ? CellEngineDataFileObjectPointer->NumberOfStencilBufferLoop = 3 : CellEngineDataFileObjectPointer->NumberOfStencilBufferLoop = 1; break;
                case GLFW_KEY_F9: CellEngineDataFileObjectPointer->DrawColorForEveryAtom = !CellEngineDataFileObjectPointer->DrawColorForEveryAtom; break;
                case GLFW_KEY_F11: CellEngineDataFileObjectPointer->DrawRandomColorForEveryParticle = !CellEngineDataFileObjectPointer->DrawRandomColorForEveryParticle; break;
                case GLFW_KEY_F10: SetVisibilityOfParticlesExcept(694, false); break;
                case GLFW_KEY_F12: SetVisibilityOfAllParticles(true); break;

                default: break;
            }
    }
    CATCH("executing on Key event for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnMouseWheel(int Pos)
{
    try
    {
        ViewZ += static_cast<float>(Pos) * CellEngineDataFileObjectPointer->ViewStep;
    }
    CATCH("executing on mouse wheel event for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnMouseButton(int Button, int Action)
{
    try
    {
        if (Button == 0)
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
    CATCH("executing on mouse Button event for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnMouseMove(int X, int Y)
{
    try
    {
        MousePosition.s.X = X;
        MousePosition.s.Y = Y;

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
        ArcBall->setBounds(static_cast<float>(Info.WindowWidth), static_cast<float>(Info.WindowHeight));

        RotationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);
    }
    CATCH("initiation of arc ball counting data for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnResize(int Width, int Height)
{
    try
    {
        Info.WindowWidth = Width;
        Info.WindowHeight = Height;
        ArcBall->setBounds(static_cast<float>(Info.WindowWidth), static_cast<float>(Info.WindowHeight));
    }
    CATCH("executing window resize event - setting bounds of arc ball counting data for cell visualisation")
}

DECLARE_MAIN(CellEngineOpenGLVisualiser)