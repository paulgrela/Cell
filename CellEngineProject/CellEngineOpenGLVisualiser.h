
#ifndef CELL_ENGINE_OPENGL_VISUALISER_H
#define CELL_ENGINE_OPENGL_VISUALISER_H

#include "ArcBall.h"

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
    vmath::mat4 RotationMatrix;
private:
    UnsignedIntType PressedRightMouseButton = 0;
private:
    bool RenderObjectsBool = true;
private:
    std::vector<std::pair<UnsignedIntType, UnsignedIntType>> BondsBetweenParticlesCentersToDraw;
    std::vector<std::vector<std::pair<UnsignedIntType, UnsignedIntType>>> BondsBetweenAtomsToDraw;
//private:
//    std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;
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
    void DrawBonds(const std::vector<CellEngineAtom>& Atoms, std::vector<std::pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix, const vmath::vec3& Center);
    void DrawBond(float x1, float y1, float z1, float x2, float y2, float z2);
protected:
    std::string GetEntityName(const uint64_t EntityId);
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
public:
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
    inline void ChooseAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, const std::vector<std::pair<uint64_t, uint64_t>>& TemporaryRenderedAtomsList, UnsignedIntType& NumberOfAllRenderedAtoms);
protected:
    [[nodiscard]] inline bool CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew) const;
};

#endif
