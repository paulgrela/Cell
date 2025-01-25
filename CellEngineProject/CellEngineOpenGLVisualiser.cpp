
#include <omp.h>

#include <sb7.h>
#include <vmath.h>
#include <shader.h>

#include <string>
#include <memory>
#include <iostream>

#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineUseful.h"
#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsEngine.h"

using namespace std;

void CellEngineOpenGLVisualiser::InitExternalData()
{
    try
    {
        CellEngineDataFileObjectPointer->ReadDataFromFile(true, true, CellEngineConfigDataObject.TypeOfFileToRead);

        GetMemoryForBondsBetweenAtomsToDraw();
    }
    CATCH("reading of data file")
}

void CellEngineOpenGLVisualiser::StartUp()
{
    try
    {
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

        GetStartCenterPoint();

        glUseProgram(ShaderProgramPhong);
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
        LoadShaders("..//shaders//per-fragment-phong.vs.glsl", "..//shaders//per-fragment-phong.fs.glsl", ShaderProgramPhong);

        Uniforms.DiffuseAlbedo = glGetUniformLocation(ShaderProgramPhong, "diffuse_albedo");
        Uniforms.SpecularAlbedo = glGetUniformLocation(ShaderProgramPhong, "specular_albedo");
        Uniforms.SpecularPower = glGetUniformLocation(ShaderProgramPhong, "specular_power");
    }
    CATCH("loading phong shaders for cell visualization")
}

void CellEngineOpenGLVisualiser::ShutDown()
{
    try
    {
        TextOverlayObject.TearDown();
    }
    CATCH("deleting of data for cell visualization")
}

void CellEngineOpenGLVisualiser::DeleteLineVertexes() const
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

        constexpr float LineVertexes[] = { 0.0, 0.0, 0.0,   1.0, 0.0, 0.0 };
        constexpr float LineNormals[] = { 1.0, 1.0 };

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

void CellEngineOpenGLVisualiser::DrawBond(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2) const
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

void CellEngineOpenGLVisualiser::FindBondsToDraw(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedInt, UnsignedInt>>& BondsToDraw)
{
    try
    {
        vector<vector<pair<UnsignedInt, UnsignedInt>>> BondsToDrawLocal;
        BondsToDrawLocal.resize(std::thread::hardware_concurrency());
        UnsignedInt AtomObjectIndex1 = 0;
        UnsignedInt AtomObjectIndex2 = 0;

        #pragma omp parallel for default(none) shared(BondsToDrawLocal, Atoms, LoggersManagerObject) private(AtomObjectIndex1, AtomObjectIndex2)
        for (AtomObjectIndex1 = 0; AtomObjectIndex1 < Atoms.size(); AtomObjectIndex1++)
            for (AtomObjectIndex2 = 0; AtomObjectIndex2 < Atoms.size(); AtomObjectIndex2++)
                if (AtomObjectIndex1 != AtomObjectIndex2)
                {
                    const auto& ParticlesCenterObject1 = Atoms[AtomObjectIndex1];
                    const auto& ParticlesCenterObject2 = Atoms[AtomObjectIndex2];
                    const float DiffX = ParticlesCenterObject2.X - ParticlesCenterObject1.X;
                    const float DiffY = ParticlesCenterObject2.Y - ParticlesCenterObject1.Y;
                    const float DiffZ = ParticlesCenterObject2.Z - ParticlesCenterObject1.Z;
                    const float VectorLength = sqrt(DiffX * DiffX + DiffY * DiffY + DiffZ * DiffZ);
                    if (VectorLength < 1.5)
                        BondsToDrawLocal[omp_get_thread_num()].emplace_back(make_pair(AtomObjectIndex1, AtomObjectIndex2));
                }

        for (const auto& BondsToDrawLocalObject : BondsToDrawLocal)
            for (const auto& BondsToDrawObject : BondsToDrawLocalObject)
                BondsToDraw.emplace_back(BondsToDrawObject);
    }
    CATCH("finding bonds")
}

void CellEngineOpenGLVisualiser::DrawBonds(const vector<CellEngineAtom>& Atoms, vector<pair<UnsignedInt, UnsignedInt>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix) const
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

                CreateUniformBlockForVertexShader(vmath::vec3(0.0, 0.0, 0.0), vmath::vec3(-1.0, -1.0, -1.0), ViewMatrix, vmath::translate(0.0f, 0.0f, 0.0f), false, false, false, false);

                DrawBond(AtomObject1.X - CellEngineConfigDataObject.CameraXPosition - Center.X(), AtomObject1.Y - CellEngineConfigDataObject.CameraYPosition - Center.Y(), AtomObject1.Z - CellEngineConfigDataObject.CameraZPosition - Center.Z(), AtomObject2.X - CellEngineConfigDataObject.CameraXPosition - Center.X(), AtomObject2.Y - CellEngineConfigDataObject.CameraYPosition - Center.Y(), AtomObject2.Z - CellEngineConfigDataObject.CameraZPosition - Center.Z());
            }
        }
    }
    CATCH("drawing bonds")
}

inline bool CellEngineOpenGLVisualiser::CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew)
{
    if (CellEngineConfigDataObject.CheckAtomVisibility == true)
    {
        if (CellEngineConfigDataObject.ViewPositionZ > CellEngineConfigDataObject.Distance)
        {
            if (CellEngineConfigDataObject.ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside == false)
                return ZNew > CellEngineConfigDataObject.CutZ && sqrt((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew)) > CellEngineConfigDataObject.Distance;
            else
                return true;
        }
        else
            return (ZNew > CellEngineConfigDataObject.ViewPositionZ + CellEngineConfigDataObject.ZLowToDrawInAtomScale && ZNew < CellEngineConfigDataObject.ViewPositionZ + CellEngineConfigDataObject.ZHighToDrawInAtomScale && XNew > CellEngineConfigDataObject.XLowToDrawInAtomScale && XNew < CellEngineConfigDataObject.XHighToDrawInAtomScale && YNew > CellEngineConfigDataObject.YLowToDrawInAtomScale && YNew < CellEngineConfigDataObject.YHighToDrawInAtomScale);
    }
    else
        return false;
}

inline void CellEngineOpenGLVisualiser::DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix)
{
    try
    {
        ModelMatrix = vmath::translate(0.0f, 0.0f, 0.0f) * vmath::scale(vmath::vec3(0.5, 0.5, 0.5));
        MatrixUniformBlockForVertexShaderPointer->Color = vmath::FromVec4ToVec3(sb7::color::Purple);
    }
    CATCH("drawing center point for data for cell visualization")
}

inline bool CellEngineOpenGLVisualiser::GetFinalVisibilityInModelWorld(const vmath::vec3& AtomPosition, UniformsBlock* MatrixUniformBlockForVertexShaderPointer, const bool CountNewPosition, const bool DrawOutsideBorder) const
{
    try
    {
        if (CountNewPosition == true)
        {
            const float XNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][0] * (AtomPosition.X() + CellEngineConfigDataObject.CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][0] * (AtomPosition.Y() + CellEngineConfigDataObject.CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][0] * (AtomPosition.Z() + CellEngineConfigDataObject.CameraZPosition - Center.Z());
            const float YNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][1] * (AtomPosition.X() + CellEngineConfigDataObject.CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][1] * (AtomPosition.Y() + CellEngineConfigDataObject.CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][1] * (AtomPosition.Z() + CellEngineConfigDataObject.CameraZPosition - Center.Z());
            const float ZNew = MatrixUniformBlockForVertexShaderPointer->MoveMatrix[0][2] * (AtomPosition.X() + CellEngineConfigDataObject.CameraXPosition - Center.X()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[1][2] * (AtomPosition.Y() + CellEngineConfigDataObject.CameraYPosition - Center.Y()) + MatrixUniformBlockForVertexShaderPointer->MoveMatrix[2][2] * (AtomPosition.Z() + CellEngineConfigDataObject.CameraZPosition - Center.Z());

            if (DrawOutsideBorder == true)
                if (CheckDistanceToDrawDetailsInAtomScale(XNew, YNew, ZNew) == true)
                {
                    MatrixUniformBlockForVertexShaderPointer->Color = vmath::FromVec4ToVec3(sb7::color::Purple);
                    return true;
                }

            return false;
        }
    }
    CATCH("getting final model position for data for cell visualization")

    return false;
}

inline bool CellEngineOpenGLVisualiser::CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional) const
{
    bool FinalVisibilityInModelWorld = false;

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

template <class T>
vector3_16 CellEngineOpenGLVisualiser::GetColor(const T& Object, const Particle& ParticleObject, const bool Chosen)
{
    vector3_16 FinalColor{};

    try
    {
        if (Chosen == true)
            FinalColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Yellow));
        else
            switch (CellEngineConfigDataObject.MakeColorsTypeObject)
            {
                case CellEngineConfigData::MakeColorsType::DrawColorForEveryAtom : FinalColor = Object.AtomColor; break;
                case CellEngineConfigData::MakeColorsType::DrawColorForEveryParticle : FinalColor = ParticleObject.ParticleColor; break;
                case CellEngineConfigData::MakeColorsType::DrawRandomColorForEveryParticleKind : FinalColor = ParticleObject.RandomParticleKindColor; break;
                case CellEngineConfigData::MakeColorsType::DrawRandomColorForEveryUniqueParticle : FinalColor = ParticleObject.UniqueParticleColor; break;
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
        switch(CellEngineConfigDataObject.SizeOfAtomsDrawingTypesObject)
        {
            #ifdef EXTENDED_RAM_MEMORY
            case CellEngineConfigData::SizeOfAtomsDrawingTypes::AtomSize : Size = vmath::vec3(AtomObject.SizeXAtom, AtomObject.SizeYAtom, AtomObject.SizeZAtom); break;
            case CellEngineConfigData::SizeOfAtomsDrawingTypes::ParticleSize : Size = vmath::vec3(AtomObject.SizeXParticle, AtomObject.SizeYParticle, AtomObject.SizeZParticle); break;
            #endif
            case CellEngineConfigData::SizeOfAtomsDrawingTypes::AutomaticChangeSize : Size = vmath::vec3(CellEngineConfigDataObject.SizeOfAtomX, CellEngineConfigDataObject.SizeOfAtomY, CellEngineConfigDataObject.SizeOfAtomZ); break;
            default : break;
        }
    }
    CATCH("getting size")

    return Size;
}

bool CellEngineOpenGLVisualiser::RenderObject(const CellEngineAtom& AtomObject, const Particle& ParticleObject, const vmath::mat4& ViewMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedInt& NumberOfAllRenderedAtoms, const bool Chosen, const bool RenderObjectParameter)
{
    bool FinalVisibilityInModelWorld{};

    try
    {
        if (RenderObjectParameter == true)
            NumberOfAllRenderedAtoms++;

        const vmath::vec3 AtomPosition = LengthUnit * AtomObject.Position();
        const vmath::vec3 SizeLocal = GetSize(AtomObject);
        const vmath::mat4 ModelMatrix = vmath::translate(AtomPosition.X() - CellEngineConfigDataObject.CameraXPosition - Center.X(), AtomPosition.Y() + CellEngineConfigDataObject.CameraYPosition - Center.Y(), AtomPosition.Z() + CellEngineConfigDataObject.CameraZPosition - Center.Z()) * vmath::scale(vmath::vec3(SizeLocal.X(), SizeLocal.Y(), SizeLocal.Z()));

        FinalVisibilityInModelWorld = CreateUniformBlockForVertexShader(AtomPosition, CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject, Chosen)), ViewMatrix, ModelMatrix, CountNewPosition, DrawCenter, DrawOutsideBorder, true);

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
        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
        {
            if (CellEngineConfigDataObject.ViewPositionZ > CellEngineConfigDataObject.Distance)
            {
                if (CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep == true)
                    CellEngineConfigDataObject.LoadOfAtomsStep = 100;
                if (CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom == true)
                    CellEngineConfigDataObject.SizeOfAtomX = CellEngineConfigDataObject.SizeOfAtomY = CellEngineConfigDataObject.SizeOfAtomZ = 3;
            }
            else
            {
                if (CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep == true)
                    CellEngineConfigDataObject.LoadOfAtomsStep = 1;
                if (CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom == true)
                    CellEngineConfigDataObject.SizeOfAtomX = CellEngineConfigDataObject.SizeOfAtomY = CellEngineConfigDataObject.SizeOfAtomZ = 1;
            }
        }
    }
    CATCH("setting automatic parameters for rendering")
}

inline void CellEngineOpenGLVisualiser::PrepareOpenGLToRenderObjectsOnScene() const
{
    try
    {
        static const GLfloat gray[] = {0.1f, 0.1f, 0.1f, 0.0f};
        static const GLfloat ones[] = {1.0f};

        glViewport(0, 0, Info.WindowWidth, Info.WindowHeight);

        glClearBufferfv(GL_COLOR, 0, gray);
        glClearBufferfv(GL_DEPTH, 0, ones);

        vmath::vec3 BackgroundColor = CellEngineConfigDataObject.BackgroundColors[CellEngineConfigDataObject.ChosenBackgroundColor];
        glClearColor(BackgroundColor.data[0], BackgroundColor.data[1], BackgroundColor.data[2], 0.0f);

        glClearStencil(0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glEnable(GL_STENCIL_TEST);
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

        glUniform1f(Uniforms.SpecularPower, CellEngineConfigDataObject.SpecularPower);
        glUniform3fv(Uniforms.SpecularAlbedo, 1, vmath::vec3(CellEngineConfigDataObject.SpecularAlbedo));
    }
    CATCH("preparing opengl to render objects on scene")
}

inline void CellEngineOpenGLVisualiser::LoadShapeOfAtomsWhenChanged()
{
    try
    {
        static int PrevShapesOfAtoms = 1;

        if (CellEngineConfigDataObject.ChosenShapeOfAtoms != PrevShapesOfAtoms)
        {
            switch (CellEngineConfigDataObject.ChosenShapeOfAtoms)
            {
                case 1 : AtomGraphicsObject.Load("..//objects//sphere.sbm"); break;
                case 2 : AtomGraphicsObject.Load("..//objects//cube.sbm"); break;
                case 3 : AtomGraphicsObject.Load("..//objects//torus.sbm"); break;
                default : break;
            }
            PrevShapesOfAtoms = CellEngineConfigDataObject.ChosenShapeOfAtoms;
        }
    }
    CATCH("preparing opengl to render objects on scene")
}

void CellEngineOpenGLVisualiser::CopyMousePositionWhenButtonPressed()
{
    if (PressedRightMouseButtonBool == true)
    {
        MousePositionLocal = MousePosition;
        PressedRightMouseButtonBool = false;
    }
};

void CellEngineOpenGLVisualiser::Render(double CurrentTime)
{
    try
    {
        glUseProgram(ShaderProgramPhong);

        CopyMousePositionWhenButtonPressed();

        CellEngineConfigDataObject.UseStencilBuffer == true ? CellEngineConfigDataObject.NumberOfStencilBufferLoops = 3 : CellEngineConfigDataObject.NumberOfStencilBufferLoops = 1;

        const auto start_time = chrono::high_resolution_clock::now();

        LoadShapeOfAtomsWhenChanged();

        SetAutomaticParametersForRendering();

        PrepareOpenGLToRenderObjectsOnScene();

        vmath::vec3 ViewPositionVector = vmath::vec3(CellEngineConfigDataObject.ViewPositionX, CellEngineConfigDataObject.ViewPositionY, CellEngineConfigDataObject.ViewPositionZ);
        vmath::mat4 ViewMatrix = vmath::lookat(ViewPositionVector, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(CellEngineConfigDataObject.RotationAngle1, CellEngineConfigDataObject.RotationAngle2, CellEngineConfigDataObject.RotationAngle3) * RotationMatrix;

        DrawBondsForParticlesCenters(BondsBetweenParticlesCentersToDraw, CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters, ViewMatrix);

        UnsignedInt NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
        UnsignedInt NumberOfAllRenderedAtoms = 0;

        RenderSpace(NumberOfAllRenderedAtoms, NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, ViewMatrix);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineDataFileObjectPointer->ShowNextStructureFromActiveFilm();

        CellEngineConfigDataObject.TimeParametersOfRenderingStr = GetDurationTimeInOneLineStr(start_time, stop_time, "Time of one frame = ", "Exception in measuring time");
        CellEngineConfigDataObject.NumberOfRenderedAtomsParametersOfRenderingStr = "NumberOfFoundParticlesCenter = " + to_string(NumberOfFoundParticlesCenterToBeRenderedInAtomDetails) + " NumberOfAllRenderedAtoms = " + to_string(NumberOfAllRenderedAtoms);
        if (CellEngineConfigDataObject.LogParametersOfRenderingToFile == true)
        {
            LoggersManagerObject.Log(STREAM(CellEngineConfigDataObject.TimeParametersOfRenderingStr));
            LoggersManagerObject.Log(STREAM(CellEngineConfigDataObject.NumberOfRenderedAtomsParametersOfRenderingStr << " ViewZ = " << to_string(CellEngineConfigDataObject.ViewPositionZ) << " CameraZPosition = " << to_string(CellEngineConfigDataObject.CameraZPosition) << " AtomSize = " << to_string(CellEngineConfigDataObject.SizeOfAtomX) << endl));
        }
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

void CellEngineOpenGLVisualiser::PrintAtomDescriptionOnScreen(CellEngineAtom& ChosenAtomObject, const Particle& ChosenParticleObject)
{
    try
    {
        glDisable(GL_CULL_FACE);

        TextOverlayObject.Clear();

        string LocalTextStr = CellEngineUseful::AtomDescriptionTextsObject.Texts[0] = string("Name = ") + ChosenAtomObject.Name + " ResName = " + ChosenAtomObject.ResName;
        if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale)
        {
            CellEngineUseful::AtomDescriptionTextsObject.Texts[1] = "Chain [" + string(ChosenAtomObject.Chain) + "]" + "(" + ChosenParticleObject.Nucleotide + ")";
            CellEngineUseful::AtomDescriptionTextsObject.Texts[1] += " SEQ = [" + ChosenParticleObject.SequenceStr + "]";
            CellEngineUseful::AtomDescriptionTextsObject.Texts[2] = "EntityId = " + to_string(ChosenParticleObject.EntityId);
            CellEngineUseful::AtomDescriptionTextsObject.Texts[3] = "Entity Name = [" + GetEntityName(ChosenParticleObject.EntityId) + "]";
            CellEngineUseful::AtomDescriptionTextsObject.Texts[4] = "Gen Index = [" + to_string(ChosenParticleObject.GenomeIndex) + "]";
            CellEngineUseful::AtomDescriptionTextsObject.Texts[5] = "Gen Index Prev = [" + to_string(ChosenParticleObject.GenomeIndexPrev) + "] Gen Index Next = [" + to_string(ChosenParticleObject.GenomeIndexNext) + "]";
            LocalTextStr += " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[1] + " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[2] + " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[3];
        }

        if (CellEngineConfigDataObject.PrintAtomDescriptionOnScreen == true)
        {
            ClearRectangleOnScreen(0, Info.WindowHeight - 16, Info.WindowWidth, 16);

            TextOverlayObject.DrawText(string("ATOM DATA: " + LocalTextStr).c_str(), 0, 0);
            TextOverlayObject.Draw();
        }

        glEnable(GL_CULL_FACE);
    }
    CATCH("printing atom description on screen")
}

string CellEngineOpenGLVisualiser::GetEntityName(const UnsignedInt EntityId)
{
    string EntityName;

    try
    {
        if (auto EntityIterator = ParticlesKindsManagerObject.ParticlesKinds.find(EntityId); EntityIterator != ParticlesKindsManagerObject.ParticlesKinds.end())
            EntityName = (EntityIterator->second.IdStr.empty() == false ? EntityIterator->second.IdStr + " " + EntityIterator->second.GraphicData.NameFromDataFile : EntityIterator->second.GraphicData.NameFromDataFile);
        else
            EntityName = "";
    }
    CATCH("getting entity name")

    return EntityName;
}

void CellEngineOpenGLVisualiser::SetVisibilityOfAllParticles(const bool VisibleParam)
{
    try
    {
        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            ParticleKindObject.second.GraphicData.Visible = VisibleParam;
    }
    CATCH("setting visibility of all particles")
}

void CellEngineOpenGLVisualiser::SetVisibilityOfParticlesExcept(const UnsignedInt EntityId, const bool VisibleParam)
{
    try
    {
        SetVisibilityOfAllParticles(VisibleParam);
        ParticlesKindsManagerObject.GetGraphicParticleKind(EntityId).Visible = !VisibleParam;
    }
    CATCH("setting visibility of particles except")
}

void CellEngineOpenGLVisualiser::OnKey(int Key, int Action)
{
    try
    {
        if (Action)
            switch (Key)
            {
                case GLFW_KEY_F1: CellEngineConfigDataObject.ImGuiDemoWindowMenu = !CellEngineConfigDataObject.ImGuiDemoWindowMenu; break;
                case GLFW_KEY_F12: CellEngineConfigDataObject.ViewChangeUsingLongStep = !CellEngineConfigDataObject.ViewChangeUsingLongStep; break;
                default: break;
            }
    }
    CATCH("executing on Key event for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnMouseWheel(int Pos)
{
    try
    {
        CellEngineConfigDataObject.ViewPositionZ += static_cast<float>(Pos) * (CellEngineConfigDataObject.ViewChangeUsingLongStep == false ? CellEngineConfigDataObject.ViewZMoveShortStep : CellEngineConfigDataObject.ViewZMoveLongStep);
    }
    CATCH("executing on mouse wheel event for cell visualisation")
}

void CellEngineOpenGLVisualiser::OnMouseButton(int Button, int Action)
{
    try
    {
        if (Button == 0)
        {
            PressedRightMouseButtonBool = true;
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
        MousePosition.s.X = static_cast<GLfloat>(X);
        MousePosition.s.Y = static_cast<GLfloat>(Y);

        if (PressedRightMouseButton == true)
        {
            Quat4fT ThisQuaternion;
            ArcBall->drag(&MousePosition, &ThisQuaternion);
            Matrix3fSetRotationFromQuat4f(&ArcBallActualRotationMatrix, &ThisQuaternion);
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
