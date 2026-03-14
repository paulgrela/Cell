
#ifndef CELL_ENGINE_OPENGL_VISUALISER_H
#define CELL_ENGINE_OPENGL_VISUALISER_H

#include "../Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include <sb7.h>
#include <object.h>
#include <sb7textoverlay.h>

#include "ArcBall.h"

struct GPUAtom
{
public:
    float X{};
    float Y{};
    float Z{};
    float ColorR{};
    float ColorG{};
    float ColorB{};
    float _padding1[2]{};
};

struct GPUParticle
{
public:
    uint32_t AtomOffset{};
    uint32_t AtomCount{};
    uint32_t _padding[3]{};
};

struct GPUAtomLocal
{
public:
    UniqueIdInt Index{};
    UniqueIdInt AtomOffset{};
    uint16_t ParticleSectorXIndex{};
    uint16_t ParticleSectorYIndex{};
    uint16_t ParticleSectorZIndex{};
};

static constexpr size_t MAX_PARTICLES = 100'000'000;
static constexpr size_t MAX_ATOMS = 100'000'000;

class CellEngineOpenGLVisualiser : public sb7::OpenGLApplication
{
public:
    sb7::GraphicObject AtomGraphicsObject{};
    sb7::TextOverlay TextOverlayObject{};

public:
    GLuint FrameBufferObject{};
    GLuint ScreenColorTexture{};
    GLuint ScreenBufferInstanceTexture{};
    GLuint ScreenDepthRenderBuffer{};

    void CreateFramebuffer();
    void DiagnoseFrameBufferObject() const;
protected:
    GLuint AtomSSBO{};
    GLuint ParticleSSBO{};
    GLuint ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO{};
private:
    GLuint LineVAO{};
    GLuint LineDataBuffer[2]{};
protected:
    GLuint ComputeShaderProgramPhong = 0;
    GLuint ShaderProgramPhong = 0;
protected:
    struct UniformsBlock
    {
        vmath::mat4 MoveMatrix{};
        vmath::vec3 Color{};
        float padding{};
    };
protected:
    std::vector<UniformsBlock> UniformsBlocks{};
    static inline std::mutex UniformsBlocksFullAtomSimulationSpaceMutexObject;
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
    std::unique_ptr<ArcBallT> ArcBall{};
protected:
    Point2fT MousePosition{};
    Point2fT MousePositionLocal{};
private:
    float LengthUnit = 1;
protected:
    vmath::vec3 Center{};
private:
    vmath::mat4 RotationMatrix{};
protected:
    UnsignedInt PressedRightMouseButton = 0;
    bool PressedRightMouseButtonBool = false;
public:
    bool RenderObjectsBool = true;
public:
    vmath::mat4 ProjectionMatrixGlobal{};
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
    void CopyMousePositionWhenButtonPressed();
protected:
    void InitArcBall();
protected:
    void InitExternalData() override;
protected:
    void InitLineVertexes();
    void DeleteLineVertexes() const;
    static void FindBondsToDraw(const std::vector<CellEngineAtom>& Atoms, std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw);
    void DrawBonds(const Particle& ParticleObject, std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw, bool DrawBonds, const vmath::mat4& ViewMatrix) const;
    void DrawBond(float x1, float y1, float z1, float x2, float y2, float z2) const;
public:
    static std::string GetEntityName(UnsignedInt EntityId);
    static void SetVisibilityOfAllParticles(bool VisibleParam);
    static void SetVisibilityOfParticlesExcept(UnsignedInt EntityId, bool VisibleParam);
protected:
    void LoadShadersPhong();
    static void LoadShaders(const char* ComputeShaderFileName, const char* VertexShaderFileName, const char* FragmentShaderFileName, GLuint& ComputeShaderProgram, GLuint& ShaderProgram);
protected:
    void StartUp() override;
    void ShutDown() override;
public:
    void Render(double CurrentTime) override;
public:
    void OnKey(int Key, int Action) override;
    void OnMouseWheel(int Pos) override;
    void OnMouseButton(int Button, int Action) override;
    void OnMouseMove(int X, int Y) override;
    void OnResize(int Width, int Height) override;
protected:
    static inline vmath::vec3 GetSize(const CellEngineAtom& AtomObject);
    template <class T> static vector3_16 GetColor(const T& Object, const Particle& ParticleObject, bool Chosen);
    static inline void DrawCenterPoint(UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, vmath::mat4& ModelMatrix);
    inline bool GetFinalVisibilityInModelWorld(const vmath::vec3& AtomPosition, UniformsBlock*  MatrixUniformBlockForVertexShaderPointer, bool CountNewPosition, bool DrawOutsideBorder) const;
    inline bool CreateUniformBlockForVertexShader(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4& ModelMatrix, bool CountNewPosition, bool DrawCenter, bool DrawOutsideBorder, bool DrawAdditional) const;
    inline bool CreateUniformBlockForVertexShader1(const vmath::vec3& Position, const vmath::vec3& Color, const vmath::mat4& ViewMatrix, vmath::mat4& ModelMatrix, bool CountNewPosition, bool DrawCenter, bool DrawOutsideBorder, bool DrawAdditional) const;
    bool RenderObject(const CellEngineAtom& AtomObject, const Particle& ParticleObject, const vmath::mat4& ViewMatrix, bool CountNewPosition, bool DrawCenter, bool DrawOutsideBorder, bool Chosen, bool RenderObjectParameter);
    static inline void SetAutomaticParametersForRendering();
    inline void PrepareOpenGLToRenderObjectsOnScene() const;
    inline void LoadShapeOfAtomsWhenChanged();
    void PrintAtomDescriptionOnScreen(CellEngineAtom& ChosenAtomObject, const Particle& ChosenParticleObject);
protected:
    virtual void RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix) = 0;
protected:
    virtual void GetStartCenterPoint() = 0;
protected:
    [[nodiscard]] static inline bool CheckDistanceToDrawDetailsInAtomScale(float XNew, float YNew, float ZNew);
    //[[nodiscard]] inline bool CheckDistanceToDrawDetailsInAtomScale(const vmath::vec3& AtomPosition, float XNew, float YNew, float ZNew) const;
};

#endif

#endif
