
#ifndef CELL_ENGINE_OPENGL_VISUALISER_H
#define CELL_ENGINE_OPENGL_VISUALISER_H

#include <object.h>
#include <sb7color.h>
#include <sb7textoverlay.h>

#include "ArcBall.h"

class CellEngineOpenGLVisualiser : public sb7::OpenGLApplication
{
public:
    sb7::GraphicObject AtomGraphicsObject;
    sb7::TextOverlay TextOverlayObject;
private:
    GLuint LineVAO;
    GLuint LineDataBuffer[2];
private:
    GLuint ShaderProgramPhong = 0;
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
    Matrix3fT ArcBallPrevRotationMatrix{};
    Matrix3fT ArcBallActualRotationMatrix{};
    std::unique_ptr<ArcBallT> ArcBall;
    Point2fT MousePosition{};
private:
    float LengthUnit = 1;
private:
    vmath::vec3 Center;
private:
    vmath::mat4 RotationMatrix;
private:
    UnsignedIntType PressedRightMouseButton = 0;
private:
    std::vector<std::pair<UnsignedIntType, UnsignedIntType>> BondsBetweenParticlesCentersToDraw;
    std::vector<std::vector<std::pair<UnsignedIntType, UnsignedIntType>>> BondsBetweenAtomsToDraw;
private:
    std::vector<std::pair<UnsignedIntType, UnsignedIntType>> TemporaryRenderedAtomsList;
public:
    bool RenderObjectsBool = true;
public:
    CellEngineOpenGLVisualiser() = default;
protected:
    void Init(int WindowWidth, int WindowHeight) override
    {
        sb7::OpenGLApplication::Init(WindowWidth, WindowHeight);

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
    static void FindBondsToDraw(const std::vector<CellEngineAtom>& Atoms, std::vector<std::pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw);
    void DrawBonds(const std::vector<CellEngineAtom>& Atoms, std::vector<std::pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix);
    void DrawBond(float x1, float y1, float z1, float x2, float y2, float z2);
public:
    static std::string GetEntityName(const UnsignedIntType EntityId);
    static bool CheckVisibilityOfParticles(UnsignedIntType EntityId);
    static void SetVisibilityOfAllParticles(bool VisibleParam);
    static void SetVisibilityOfParticlesExcept(UnsignedIntType EntityId, bool VisibleParam);
protected:
    void LoadShadersPhong();
    static void LoadShaders(const char* VertexShaderFileName, const char* FragmentShaderFileName, GLuint& ShaderProgram);
protected:
    void StartUp() override;
    void ShutDown() override;
    void Render(double CurrentTime) override;
public:
    void OnKey(int Key, int Action) override;
    void OnMouseWheel(int Pos) override;
    void OnMouseButton(int Button, int Action) override;
    void OnMouseMove(int X, int Y) override;
    void OnResize(int Width, int Height) override;
protected:
    static inline vmath::vec3 GetSize(const CellEngineAtom& AtomObject);
    static inline vmath::vec3 GetColor(const CellEngineAtom& AtomObject, bool Chosen);
    static inline void DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix);
    inline bool GetFinalVisibilityInModelWorld(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, const bool CountNewPosition, const bool DrawOutsideBorder) const;
    inline bool CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4 ModelMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional);
    inline bool RenderObject(const CellEngineAtom& AtomObject, const vmath::mat4& ViewMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, UnsignedIntType& NumberOfAllRenderedAtoms, const bool Chosen, const bool RenderObjectParameter);
    inline void SetAutomaticParametersForRendering();
    inline void PrepareOpenGLToRenderObjectsOnScene();
    inline void LoadShapeOfAtomsWhenChanged();
    inline void PrintAtomDescriptionOnScreen(CellEngineAtom& ChosenParticleObject);
    inline void ChooseAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, const std::vector<std::pair<UnsignedIntType, UnsignedIntType>>& TemporaryRenderedAtomsList, UnsignedIntType& NumberOfAllRenderedAtoms);
protected:
    [[nodiscard]] static inline bool CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew);
};

#endif
