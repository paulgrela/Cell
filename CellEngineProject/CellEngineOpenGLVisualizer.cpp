

#include <omp.h>

#include <sb7.h>
#include <sb7color.h>
#include <vmath.h>
#include <object.h>
#include <shader.h>
#include <sb7textoverlay.h>

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

#include "CellEngineConfigurationFileReaderWriter.h"

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
        vmath::mat4 ProjectionMatrix;
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
    sb7::TextOverlay TextOverlayObject;
private:
    Matrix3fT ArcBallPrevRotationMatrix{};
    Matrix3fT ArcBallActualRotationMatrix{};
    std::unique_ptr<ArcBallT> ArcBall;
    Point2fT MousePosition{};
private:
    float LengthUnit = 1;
private:
    vmath::vec3 Center;
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
    bool RenderObjectsBool = true;
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
        sb7::OpenGLApplication::Init();

        static const char title[] = "Cell Engine Visualizer";
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
    string GetEntityName(const uint64_t EntityId);
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
    void ShutDown() override;
    void Render(double CurrentTime) override;
    void OnKey(int Key, int Action) override;
    void OnMouseWheel(int Pos) override;
    void OnMouseButton(int Button, int Action) override;
    void OnMouseMove(int X, int Y) override;
    void OnResize(int Width, int Height) override;
protected:
    inline vmath::vec3 GetSize(const CellEngineAtom& AtomObject);
    inline vmath::vec3 GetColor(const CellEngineAtom& AtomObject, bool Chosen);
    static inline void DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix);
    inline bool GetFinalVisibilityInModelWorld(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, const bool CountNewPosition, const bool DrawOutsideBorder) const;
    inline bool CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional);
    inline bool RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedIntType& NumberOfAllRenderedAtoms, const bool Chosen, const bool RenderObjectParameter);
    inline void SetAutomaticParametersForRendering();
    inline void PrepareOpenGLToRenderObjectsOnScene();
    inline void PrintAtomDescriptionOnScreen(CellEngineAtom& ChosenParticleObject);
    inline void ChooseAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, const vector<pair<uint64_t, uint64_t>>& TemporaryRenderedAtomsList, UnsignedIntType& NumberOfAllRenderedAtoms);
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

        uint64_t ExecuteCellStateId;
        if (__argc > 1)
            ExecuteCellStateId = stoi(__argv[1]);
        else
            LoggersManagerObject.Log(STREAM("Lack of cell id to execute in program parameters"));

        CellEngineConfigurationFileReaderWriterObject.ReadChessConfigurationFile("CellEngineProjectConfig.xml", CellEngineDataFileObjectPointer, ExecuteCellStateId);
        CellEngineDataFileObjectPointer->ReadDataFromFile();

        BondsBetweenAtomsToDraw.resize(CellEngineDataFileObjectPointer->GetParticlesCenters().size());
    }
    CATCH("reading of data file")
}

void CellEngineOpenGLVisualiser::StartUp()
{
    try
    {
        LoadShadersSimple();
        LoadShadersPhong();

        TextOverlayObject.Init(160, 80, "..//textures//cp437_9x16.ktx");

        glGenBuffers(1, &UniformsBuffer);
        glBindBuffer(GL_UNIFORM_BUFFER, UniformsBuffer);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformsBlock), nullptr, GL_DYNAMIC_DRAW);

        AtomGraphicsObject.Load("..//objects//sphere.sbm");
        InitLineVertexes();

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        InitArcBall();

        Center = CellEngineDataFileObjectPointer->GetCenter(CellEngineDataFileObjectPointer->GetParticlesCenters());

        CameraXPosition = CellEngineDataFileObjectPointer->CameraXPosition;
        CameraYPosition = CellEngineDataFileObjectPointer->CameraYPosition;
        CameraZPosition = CellEngineDataFileObjectPointer->CameraZPosition;
    }
    CATCH("initiation of data for cell visualization")
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

void CellEngineOpenGLVisualiser::ShutDown()
{
    try
    {
        TextOverlayObject.TearDown();
    }
    CATCH("deleting of data for cell visualization")
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
        vector<vector<pair<UnsignedIntType, UnsignedIntType>>> BondsToDrawLocal;
        BondsToDrawLocal.resize(std::thread::hardware_concurrency());
        UnsignedIntType AtomObjectIndex1 = 0;
        UnsignedIntType AtomObjectIndex2 = 0;
        #pragma omp parallel for default(none) shared(BondsToDrawLocal, Atoms, LoggersManagerObject) private(AtomObjectIndex1, AtomObjectIndex2)
        for (AtomObjectIndex1 = 0; AtomObjectIndex1 < Atoms.size(); AtomObjectIndex1++)
            for (AtomObjectIndex2 = 0; AtomObjectIndex2 < Atoms.size(); AtomObjectIndex2++)
                if (AtomObjectIndex1 != AtomObjectIndex2)
                {
                    const auto& ParticlesCenterObject1 = Atoms[AtomObjectIndex1];
                    const auto& ParticlesCenterObject2 = Atoms[AtomObjectIndex2];

                    float DiffX = ParticlesCenterObject2.X - ParticlesCenterObject1.X;
                    float DiffY = ParticlesCenterObject2.Y - ParticlesCenterObject1.Y;
                    float DiffZ = ParticlesCenterObject2.Z - ParticlesCenterObject1.Z;
                    float VectorLength = sqrt(DiffX * DiffX + DiffY * DiffY + DiffZ * DiffZ);
                    if (VectorLength < 1.5)
                        BondsToDrawLocal[omp_get_thread_num()].emplace_back(make_pair(AtomObjectIndex1, AtomObjectIndex2));
                }

        for (const auto& BondsToDrawLocalObject : BondsToDrawLocal)
            for (const auto& BondsToDrawObject : BondsToDrawLocalObject)
                BondsToDraw.emplace_back(BondsToDrawObject);
    }
    CATCH("finding bonds")
}

void CellEngineOpenGLVisualiser::DrawBonds(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix, const vmath::vec3& Center)
{
    try
    {
        if (DrawBonds == true)
        {
            glUseProgram(ShaderProgramSimple);

            if (BondsToDraw.empty() == true)
                FindBondsToDraw(Atoms, BondsToDraw);

            for (const auto& BondToDrawObject : BondsToDraw)
            {
                const auto& AtomObject1 = Atoms[BondToDrawObject.first];
                const auto& AtomObject2 = Atoms[BondToDrawObject.second];

                CreateUniformBlockForVertexShader(vmath::vec3(0.0, 0.0, 0.0), sb7::FromVec4ToVec3(sb7::color::DarkCyan), ViewMatrix, vmath::translate(0.0f, 0.0f, 0.0f), false, false, false, false);

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
            return ZNew > CellEngineDataFileObjectPointer->CutZ && sqrt((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew)) > CellEngineDataFileObjectPointer->Distance;
        else
            return (ZNew > ViewZ + CellEngineDataFileObjectPointer->ZLowToDrawInAtomScale && ZNew < ViewZ + CellEngineDataFileObjectPointer->ZHighToDrawInAtomScale && XNew > CellEngineDataFileObjectPointer->XLowToDrawInAtomScale && XNew < CellEngineDataFileObjectPointer->XHighToDrawInAtomScale && YNew > CellEngineDataFileObjectPointer->YLowToDrawInAtomScale && YNew < CellEngineDataFileObjectPointer->YHighToDrawInAtomScale);
    }
    else
        return false;
}

inline void CellEngineOpenGLVisualiser::DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix)
{
    try
    {
        ModelMatrix = vmath::translate(0.0f, 0.0f, 0.0f) * vmath::scale(vmath::vec3(0.5, 0.5, 0.5));
        MatrixUniformBlockForVertexShaderPointer->Color = sb7::FromVec4ToVec3(sb7::color::Purple);
    }
    CATCH("drawing center point for data for cell visualization")
}

inline bool CellEngineOpenGLVisualiser::GetFinalVisibilityInModelWorld(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, const bool CountNewPosition, const bool DrawOutsideBorder) const
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
                {
                    MatrixUniformBlockForVertexShaderPointer->Color = sb7::FromVec4ToVec3(sb7::color::Purple);
                    return true;
                }

            return false;
        }
    }
    CATCH("getting final model position for data for cell visualization")

    return false;
}

inline bool CellEngineOpenGLVisualiser::CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional)
{
    bool FinalVisibilityInModelWorld;

    try
    {
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, UniformsBuffer);
        auto MatrixUniformBlockForVertexShaderPointer = (UniformsBlock*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(UniformsBlock), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

        MatrixUniformBlockForVertexShaderPointer->ProjectionMatrix = vmath::perspective(50.0f, (float)Info.WindowWidth / (float)Info.WindowHeight, 0.1f, 10000.0f);
        MatrixUniformBlockForVertexShaderPointer->Color = Color;

        if (DrawAdditional == true)
        {
            FinalVisibilityInModelWorld = GetFinalVisibilityInModelWorld(Position, MatrixUniformBlockForVertexShaderPointer, CountNewPosition, DrawOutsideBorder);
            if (DrawCenter == true)
                DrawCenterPoint(MatrixUniformBlockForVertexShaderPointer, ModelMatrix);
        }

        MatrixUniformBlockForVertexShaderPointer->MoveMatrix = ViewMatrix * ModelMatrix;

        glUnmapBuffer(GL_UNIFORM_BUFFER);
    }
    CATCH("rendering object for data for cell visualization")

    return FinalVisibilityInModelWorld;
}

inline vmath::vec3 CellEngineOpenGLVisualiser::GetColor(const CellEngineAtom& AtomObject, bool Chosen)
{
    vmath::vec3 FinalColor;

    try
    {
        if (Chosen == true)
            FinalColor = sb7::FromVec4ToVec3(sb7::color::Yellow);
        else
            switch (CellEngineDataFileObjectPointer->MakeColorsTypeObject)
            {
                case CellEngineDataFile::MakeColorsType::DrawColorForEveryAtom : FinalColor = AtomObject.AtomColor; break;
                case CellEngineDataFile::MakeColorsType::DrawColorForEveryParticle : FinalColor = AtomObject.ParticleColor; break;
                case CellEngineDataFile::MakeColorsType::DrawRandomColorForEveryParticle : FinalColor = AtomObject.RandomParticleColor; break;
                default : break;
            }
    }
    CATCH("getting color")

    return FinalColor;
}

inline vmath::vec3 CellEngineOpenGLVisualiser::GetSize(const CellEngineAtom& AtomObject)
{
    vmath::vec3 Size;

    try
    {
        switch(CellEngineDataFileObjectPointer->SizeOfAtomsDrawingTypesObject)
        {
            case CellEngineDataFile::SizeOfAtomsDrawingTypes::AtomSize : Size = vmath::vec3(AtomObject.SizeXAtom, AtomObject.SizeYAtom, AtomObject.SizeZAtom); break;
            case CellEngineDataFile::SizeOfAtomsDrawingTypes::ParticleSize : Size = vmath::vec3(AtomObject.SizeXParticle, AtomObject.SizeYParticle, AtomObject.SizeZParticle); break;
            case CellEngineDataFile::SizeOfAtomsDrawingTypes::AutomaticChangeSize : Size = vmath::vec3(CellEngineDataFileObjectPointer->SizeOfAtomX, CellEngineDataFileObjectPointer->SizeOfAtomY, CellEngineDataFileObjectPointer->SizeOfAtomZ); break;
            default : break;
        }
    }
    CATCH("getting size")

    return Size;
}

inline bool CellEngineOpenGLVisualiser::RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedIntType& NumberOfAllRenderedAtoms, const bool Chosen, const bool RenderObjectParameter)
{
    bool FinalVisibilityInModelWorld;

    try
    {
        if (RenderObjectParameter == true)
            NumberOfAllRenderedAtoms++;

        vmath::vec3 AtomPosition = LengthUnit * AtomObject.Position();
        vmath::vec3 SizeLocal = GetSize(AtomObject);
        vmath::mat4 ModelMatrix = vmath::translate(AtomPosition.X() - CameraXPosition - Center.X(), AtomPosition.Y() + CameraYPosition - Center.Y(), AtomPosition.Z() + CameraZPosition - Center.Z()) * vmath::scale(vmath::vec3(SizeLocal.X(), SizeLocal.Y(), SizeLocal.Z()));

        FinalVisibilityInModelWorld = CreateUniformBlockForVertexShader(AtomPosition, GetColor(AtomObject, Chosen), ViewMatrix, ModelMatrix, CountNewPosition, DrawCenter, DrawOutsideBorder, true);

        if (RenderObjectParameter == true)
            AtomGraphicsObject.Render();
    }
    CATCH("rendering object for data for cell visualization")

    return FinalVisibilityInModelWorld;
}

inline void CellEngineOpenGLVisualiser::SetAutomaticParametersForRendering()
{
    try
    {
        if (CellEngineDataFileObjectPointer->ShowDetailsInAtomScale == true)
        {
            if (ViewZ > CellEngineDataFileObjectPointer->Distance)
            {
                if (CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep == true)
                    CellEngineDataFileObjectPointer->LoadOfAtomsStep = 100;
                if (CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom == true)
                    CellEngineDataFileObjectPointer->SizeOfAtomX = CellEngineDataFileObjectPointer->SizeOfAtomY = CellEngineDataFileObjectPointer->SizeOfAtomZ = 3;
            }
            else
            {
                if (CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep == true)
                    CellEngineDataFileObjectPointer->LoadOfAtomsStep = 1;
                if (CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom == true)
                    CellEngineDataFileObjectPointer->SizeOfAtomX = CellEngineDataFileObjectPointer->SizeOfAtomY = CellEngineDataFileObjectPointer->SizeOfAtomZ = 1;
            }
        }
    }
    CATCH("setting automatic parameters for rendering")
}

inline void CellEngineOpenGLVisualiser::PrepareOpenGLToRenderObjectsOnScene()
{
    try
    {
        static const GLfloat gray[] = {0.1f, 0.1f, 0.1f, 0.0f};
        static const GLfloat ones[] = {1.0f};

        glUseProgram(ShaderProgramPhong);
        glDisable(GL_SCISSOR_TEST);

        glViewport(0, 0, Info.WindowWidth, Info.WindowHeight);

        glClearBufferfv(GL_COLOR, 0, gray);
        glClearBufferfv(GL_DEPTH, 0, ones);

        vmath::vec3 BackgroundColor = CellEngineDataFileObjectPointer->BackgroundColors[CellEngineDataFileObjectPointer->ChosenBackgroundColor];
        glClearColor(BackgroundColor.data[0], BackgroundColor.data[1], BackgroundColor.data[2], 0.0f);

        glClearStencil(0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glEnable(GL_STENCIL_TEST);
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

        glUniform1f(Uniforms.SpecularPower, CellEngineDataFileObjectPointer->SpecularPower);
        glUniform3fv(Uniforms.SpecularAlbedo, 1, vmath::vec3(CellEngineDataFileObjectPointer->SpecularAlbedo));
    }
    CATCH("preparing opengl to render objects on scene")
}

void CellEngineOpenGLVisualiser::Render(double CurrentTime)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        SetAutomaticParametersForRendering();

        PrepareOpenGLToRenderObjectsOnScene();

        vmath::vec3 ViewPositionVector = vmath::vec3(ViewX, ViewY, ViewZ);
        vmath::mat4 ViewMatrix = vmath::lookat(ViewPositionVector, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix;

        DrawBonds(CellEngineDataFileObjectPointer->GetParticlesCenters(), BondsBetweenParticlesCentersToDraw, CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters, ViewMatrix, Center);

        UnsignedIntType NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
        UnsignedIntType NumberOfAllRenderedAtoms = 0;

        vector<pair<uint64_t, uint64_t>> TemporaryRenderedAtomsList;
        glUseProgram(ShaderProgramPhong);
        GLuint PartOfStencilBufferIndex[3];

        for (uint64_t StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
            NumberOfAllRenderedAtoms = 0;

            if (StencilBufferLoopCounter > 0)
            {
                glEnable(GL_SCISSOR_TEST);
                glScissor(MousePosition.s.X, (float)Info.WindowHeight - MousePosition.s.Y - 1, 1, 1);
            }
            else
                glDisable(GL_SCISSOR_TEST);

            TemporaryRenderedAtomsList.clear();

            for (auto ParticlesCenterIterator = CellEngineDataFileObjectPointer->GetParticlesCenters().begin(); ParticlesCenterIterator != CellEngineDataFileObjectPointer->GetParticlesCenters().end(); ++ParticlesCenterIterator)
            {
                if (CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject == CellEngineDataFile::StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters)
                    glStencilFunc(GL_ALWAYS, uint8_t((NumberOfAllRenderedAtoms) >> (8 * StencilBufferLoopCounter)), -1);

                auto ParticlesCenterObject = *ParticlesCenterIterator;

                bool FinalVisibilityInModelWorld = RenderObject(ParticlesCenterObject, ViewMatrix, true, ParticlesCenterIterator == CellEngineDataFileObjectPointer->GetParticlesCenters().end() - 1, true, NumberOfAllRenderedAtoms, false, !CellEngineDataFileObjectPointer->ShowDetailsInAtomScale);

                if (CellEngineDataFileObjectPointer->ShowDetailsInAtomScale == true)
                    if (FinalVisibilityInModelWorld == true)
                        if (CheckVisibilityOfParticles(ParticlesCenterObject.EntityId) == true)
                        {
                            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                            DrawBonds(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex], BondsBetweenAtomsToDraw[ParticlesCenterObject.AtomIndex], CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms, ViewMatrix, Center);

                            glUseProgram(ShaderProgramPhong);
                            UnsignedIntType AtomObjectIndex;
                            for (AtomObjectIndex = 0; AtomObjectIndex < CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex].size(); AtomObjectIndex += CellEngineDataFileObjectPointer->LoadOfAtomsStep)
                            {
                                if (CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject == CellEngineDataFile::StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale)
                                {
                                    uint8_t ToInsert = (TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter);
                                    glStencilFunc(GL_ALWAYS, ToInsert, -1);
                                    TemporaryRenderedAtomsList.emplace_back(make_pair(ParticlesCenterObject.AtomIndex, AtomObjectIndex));
                                }
                                RenderObject(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex][AtomObjectIndex], ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                            }
                        }
            }
            CellEngineDataFileObjectPointer->ShowNextStructureFromActiveFilm();

            if (CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops > 1)
            {
                GLuint StencilIndex;
                glReadPixels(MousePosition.s.X, (float)Info.WindowHeight - MousePosition.s.Y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &StencilIndex);
                PartOfStencilBufferIndex[StencilBufferLoopCounter] = StencilIndex;
            }
        }

        ChooseAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, TemporaryRenderedAtomsList, NumberOfAllRenderedAtoms);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Time of one frame = ", "Exception in measuring time")));
        LoggersManagerObject.Log(STREAM("NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = " << to_string(NumberOfFoundParticlesCenterToBeRenderedInAtomDetails) << " NumberOfAllRenderedAtoms = " << to_string(NumberOfAllRenderedAtoms) << " ViewZ = " << to_string(ViewZ) << " CameraZPosition = " << to_string(CameraZPosition) << " AtomSize = " << to_string(CellEngineDataFileObjectPointer->SizeOfAtomX) << endl));
    }
    CATCH("rendering cell visualization")
}

inline void ClearRectangleOnScreen(const GLint XStart, const GLint YStart, const GLsizei  XWidth, const GLsizei YHeight)
{
    try
    {
        glEnable(GL_SCISSOR_TEST);
        glScissor(XStart, YStart, XWidth, YHeight);
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
    }
    CATCH("filling rectangle on screen")
}

inline void CellEngineOpenGLVisualiser::PrintAtomDescriptionOnScreen(CellEngineAtom& ChosenParticleObject)
{
    try
    {
        glDisable(GL_CULL_FACE);
        TextOverlayObject.Clear();
        string AtomDescription = "ATOM DATA: Serial = " + to_string(ChosenParticleObject.Serial) + " Name = " + ChosenParticleObject.Name + " ResName = " + ChosenParticleObject.ResName;
        if (CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject == CellEngineDataFile::StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale)
        AtomDescription += " Chain [" + string(ChosenParticleObject.Chain) + "] EntityId = " + to_string(ChosenParticleObject.EntityId) + " Entity Name = [" + GetEntityName(ChosenParticleObject.EntityId) + "]";
        TextOverlayObject.DrawText(AtomDescription.c_str(), 0, 0);
        TextOverlayObject.Draw();
        glEnable(GL_CULL_FACE);
    }
    CATCH("printing atom description on screen")
}

inline void CellEngineOpenGLVisualiser::ChooseAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, const vector<pair<uint64_t, uint64_t>>& TemporaryRenderedAtomsList, UnsignedIntType& NumberOfAllRenderedAtoms)
{
    try
    {
        if (CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops > 1)
        {
            uint64_t ChosenParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

            if (ChosenParticleCenterIndex > 0)
            {
                CellEngineAtom ChosenParticleObject;
                if (CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject == CellEngineDataFile::StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters)
                {
                    if (ChosenParticleCenterIndex > CellEngineDataFileObjectPointer->GetParticlesCenters().size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + to_string(ChosenParticleCenterIndex) + " MAXIMAL NUMBER OF OBJECTS = " + to_string(CellEngineDataFileObjectPointer->GetParticlesCenters().size()));
                    else
                        ChosenParticleObject = CellEngineDataFileObjectPointer->GetParticlesCenters()[ChosenParticleCenterIndex];
                }
                else
                {
                    if (TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first > CellEngineDataFileObjectPointer->GetAllAtoms().size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + to_string(TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first) + " MAXIMAL NUMBER OF OBJECTS = " + to_string(CellEngineDataFileObjectPointer->GetAllAtoms().size()));
                    else
                    {
                        if (TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second > CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first].size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + to_string(TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second) + " MAXIMAL NUMBER OF OBJECTS = " + to_string(CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first].size()));
                        else
                            ChosenParticleObject = CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first][TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second];
                    }
                }

                glDisable(GL_SCISSOR_TEST);

                RenderObject(ChosenParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

                ClearRectangleOnScreen(0, Info.WindowHeight - 16, Info.WindowWidth, 16);

                PrintAtomDescriptionOnScreen(ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

string CellEngineOpenGLVisualiser::GetEntityName(const uint64_t EntityId)
{
    string EntityName;

    try
    {
        auto EntityIterator = CellEngineDataFileObjectPointer->ParticlesKinds.find(EntityId);
        if (EntityIterator != CellEngineDataFileObjectPointer->ParticlesKinds.end())
            EntityName = EntityIterator->second.Name;
        else
            EntityName = "";
    }
    CATCH("getting entity name")

    return EntityName;
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
                case '1': CameraZPosition += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraZMoveShortStep : CellEngineDataFileObjectPointer->CameraZMoveLongStep); break;
                case '2': CameraZPosition -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraZMoveShortStep : CellEngineDataFileObjectPointer->CameraZMoveLongStep); break;
                case '3': CameraXPosition += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraXMoveShortStep : CellEngineDataFileObjectPointer->CameraXMoveLongStep); break;
                case '4': CameraXPosition -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraZMoveShortStep : CellEngineDataFileObjectPointer->CameraXMoveLongStep); break;
                case '5': CameraYPosition += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraYMoveShortStep : CellEngineDataFileObjectPointer->CameraYMoveLongStep); break;
                case '6': CameraYPosition -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->CameraYMoveShortStep : CellEngineDataFileObjectPointer->CameraYMoveLongStep); break;

                case 'Q': ViewX += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewXMoveShortStep : CellEngineDataFileObjectPointer->ViewXMoveLongStep); break;
                case 'W': ViewX -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewXMoveShortStep : CellEngineDataFileObjectPointer->ViewXMoveLongStep); break;
                case 'E': ViewY += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewYMoveShortStep : CellEngineDataFileObjectPointer->ViewYMoveLongStep); break;
                case 'R': ViewY -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewYMoveShortStep : CellEngineDataFileObjectPointer->ViewYMoveLongStep); break;
                case 'T': ViewZ += (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewZMoveShortStep : CellEngineDataFileObjectPointer->ViewZMoveLongStep); break;
                case 'Y': ViewZ -= (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewZMoveShortStep : CellEngineDataFileObjectPointer->ViewZMoveLongStep); break;

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
                case 'Z': CellEngineDataFileObjectPointer->SizeOfAtomX -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeOfAtomY -= CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeOfAtomZ -= CellEngineDataFileObjectPointer->SizeStep; break;
                case 'X': CellEngineDataFileObjectPointer->SizeOfAtomX += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeOfAtomY += CellEngineDataFileObjectPointer->SizeStep; CellEngineDataFileObjectPointer->SizeOfAtomZ += CellEngineDataFileObjectPointer->SizeStep; break;

                case 'U': CellEngineDataFileObjectPointer->LoadOfAtomsStep == 100 ? CellEngineDataFileObjectPointer->LoadOfAtomsStep = 10 : CellEngineDataFileObjectPointer->LoadOfAtomsStep = 100;  break;
                case 'I': CellEngineDataFileObjectPointer->LoadOfAtomsStep = 1;  break;

                case '9': AtomGraphicsObject.Load("..//objects//cube.sbm");  break;
                case '0': AtomGraphicsObject.Load("..//objects//sphere.sbm");  break;

                case GLFW_KEY_F1: CellEngineDataFileObjectPointer->ChosenBackgroundColor = 1; break;
                case GLFW_KEY_F2: CellEngineDataFileObjectPointer->ChosenBackgroundColor = 3; break;

                case GLFW_KEY_F3: RenderObjectsBool = !RenderObjectsBool; break;
                case GLFW_KEY_F4: CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom = !CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom; break;
                case GLFW_KEY_F5: CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep = !CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep; break;
                case GLFW_KEY_F6: CellEngineDataFileObjectPointer->ViewChangeUsingLongStep = !CellEngineDataFileObjectPointer->ViewChangeUsingLongStep; break;
                case GLFW_KEY_F7: CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops == 1 ? CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops = 3 : CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops = 1; break;

                case GLFW_KEY_F8: CellEngineDataFileObjectPointer->MakeColorsTypeObject = CellEngineDataFile::MakeColorsType::DrawColorForEveryAtom; break;
                case GLFW_KEY_F9: CellEngineDataFileObjectPointer->MakeColorsTypeObject = CellEngineDataFile::MakeColorsType::DrawColorForEveryParticle; break;
                case GLFW_KEY_F10: CellEngineDataFileObjectPointer->MakeColorsTypeObject = CellEngineDataFile::MakeColorsType::DrawRandomColorForEveryParticle; break;

                case GLFW_KEY_F11: SetVisibilityOfParticlesExcept(694, false); break;
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
        ViewZ += static_cast<float>(Pos) * (CellEngineDataFileObjectPointer->ViewChangeUsingLongStep == false ? CellEngineDataFileObjectPointer->ViewZMoveShortStep : CellEngineDataFileObjectPointer->ViewZMoveLongStep);
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