
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
#include "AdditionalFunctions.h"

using namespace std;

class phonglighting_app : public sb7::application
{
public:
    phonglighting_app() : per_fragment_program(0), per_vertex_program(0), per_vertex(false)
    {
    }

protected:
    void init()
    {
        static const char title[] = "Cell Engine Visualizer - OpenGL";

        sb7::application::init();

        memcpy(info.title, title, sizeof(title));
    }

    void InitArcBall();

    void startup();
    void load_shaders();
    void render(double currentTime);
    void onKey(int key, int action);
    void onMouseWheel(int pos);    
    void onMouseButton(int button, int action);
    void onMouseMove(int x, int y);
    void onResize(int w, int h);


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

    //float CameraXPosition = 40.0;
    float CameraXPosition = 0.0;
    float CameraYPosition = 0.0;
    float CameraZPosition = 0.0;
    //float CameraZPosition = 50.0;

    float ViewZ = 50.0;

    vmath::mat4 RotationMatrix;

    float RotationAngle1 = 0.0;
    float RotationAngle2 = 0.0;
    float RotationAngle3 = 0.0;
};

const uint64_t MaxX = 170 / 1;
const uint64_t MaxY = 170 / 1;
const uint64_t MaxZ = 7;

vmath::vec3 Colors[MaxX][MaxY][MaxZ];

void GenereteRandomColors()
{
    for (float a = 0.1f, j = 0; j < MaxX; j++, a < 1 ? a++ : a = 0.1f)
        for (float b = 0.1f, i = 0; i < MaxY; i++, b < 1 ? b++ : b = 0.1f)
            for (float c = 0.1f, z = 0; z < MaxZ; z++, c < 1 ? c++ : c = 0.1f)
                Colors[uint64_t(j)][uint64_t(i)][uint64_t(z)] = vmath::vec3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);
}

void phonglighting_app::load_shaders()
{
    GLuint vs;
    GLuint fs;

    vs = sb7::shader::load("per-fragment-phong.vs.glsl", GL_VERTEX_SHADER);
    fs = sb7::shader::load("per-fragment-phong.fs.glsl", GL_FRAGMENT_SHADER);

    if (per_fragment_program)
        glDeleteProgram(per_fragment_program);

    per_fragment_program = glCreateProgram();
    glAttachShader(per_fragment_program, vs);
    glAttachShader(per_fragment_program, fs);
    glLinkProgram(per_fragment_program);

    uniforms[0].diffuse_albedo = glGetUniformLocation(per_fragment_program, "diffuse_albedo");
    uniforms[0].specular_albedo = glGetUniformLocation(per_fragment_program, "specular_albedo");
    uniforms[0].specular_power = glGetUniformLocation(per_fragment_program, "specular_power");

    vs = sb7::shader::load("per-vertex-phong.vs.glsl", GL_VERTEX_SHADER);
    fs = sb7::shader::load("per-vertex-phong.fs.glsl", GL_FRAGMENT_SHADER);

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

void phonglighting_app::startup()
{
    load_shaders();

    glGenBuffers(1, &uniforms_buffer);
    glBindBuffer(GL_UNIFORM_BUFFER, uniforms_buffer);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(uniforms_block), NULL, GL_DYNAMIC_DRAW);

    object.load("sphere.sbm");

    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    GenereteRandomColors();

    InitArcBall();
}

void phonglighting_app::render(double currentTime)
{
    static const GLfloat zeros[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    static const GLfloat gray[] = { 0.1f, 0.1f, 0.1f, 0.0f };
    static const GLfloat ones[] = { 1.0f };
    const float f = (float)currentTime;

    glUseProgram(per_vertex ? per_vertex_program : per_fragment_program);
    glViewport(0, 0, info.windowWidth, info.windowHeight);

    glClearBufferfv(GL_COLOR, 0, gray);
    glClearBufferfv(GL_DEPTH, 0, ones);

    //vmath::vec3 view_position = vmath::vec3(0.0f, 0.0f, 0.0f);
    vmath::vec3 view_position = vmath::vec3(0.0f, 0.0f, ViewZ);
    //vmath::vec3 view_position = vmath::vec3(CameraXPosition, CameraYPosition, CameraZPosition);
    //vmath::mat4 view_matrix = vmath::lookat(view_position, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix * vmath::translate(CameraXPosition, CameraYPosition, CameraZPosition);
    vmath::mat4 view_matrix = vmath::lookat(view_position, vmath::vec3(0.0f, 0.0f, 0.0f), vmath::vec3(0.0f, 1.0f, 0.0f)) * vmath::rotate(RotationAngle1, RotationAngle2, RotationAngle3) * RotationMatrix;

    vmath::vec3 light_position = vmath::vec3(0.0f, 0.0f, 100.0f);

    vmath::mat4 light_proj_matrix = vmath::frustum(-1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 200.0f);
    vmath::mat4 light_view_matrix = vmath::lookat(light_position, vmath::vec3(0.0f), vmath::vec3(0.0f, 1.0f, 0.0f));

    for (float a = 0.1f, j = 0; j < MaxX; j++, a < 1 ? a++ : a = 0.1f)
        for (float b = 0.1f, i = 0; i < MaxY; i++, b < 1 ? b++ : b = 0.1f)
            for (float c = 0.1f, z = 0; z < MaxZ; z++, c < 1 ? c++ : c = 0.1f)
            {
                glBindBufferBase(GL_UNIFORM_BUFFER, 0, uniforms_buffer);
                uniforms_block* block = (uniforms_block*)glMapBufferRange(GL_UNIFORM_BUFFER, 0, sizeof(uniforms_block), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);

                vmath::mat4 model_matrix = vmath::translate(i * 3.0f - CameraXPosition, j * 3.0f + CameraYPosition, z * 3.0f + CameraZPosition);

                block->mv_matrix = view_matrix * model_matrix;
                block->view_matrix = view_matrix;
                block->proj_matrix = vmath::perspective(50.0f, (float)info.windowWidth / (float)info.windowHeight, 0.1f, 2000.0f);

                block->color = Colors[uint64_t(j)][uint64_t(i)][uint64_t(z)];

                glUnmapBuffer(GL_UNIFORM_BUFFER);

                glUniform1f(uniforms[per_vertex ? 1 : 0].specular_power, powf(2.0f, (float)j + 2.0f));
                glUniform3fv(uniforms[per_vertex ? 1 : 0].specular_albedo, 1, vmath::vec3((float)i / 9.0f + 1.0f / 9.0f));

                object.render();
            }
}

void phonglighting_app::onKey(int key, int action)
{
    if (action)
        switch (key)
        {
            case 'M': load_shaders(); break;
            case 'V': per_vertex = !per_vertex; break;


            //case '1': CameraZPosition += 1; CameraXPosition -= 1; break;
            //case '2': CameraZPosition -= 1; CameraXPosition += 1; break;

            case '1': CameraZPosition += 1; break;
            case '2': CameraZPosition -= 1; break;
            case '3': CameraXPosition += 1; break;
            case '4': CameraXPosition -= 1; break;
            case '5': CameraYPosition += 1; break;
            case '6': CameraYPosition -= 1; break;

            case 'Q': ViewZ += 1; break;
            case 'W': ViewZ -= 1; break;
            
            case 'E': RotationAngle1 += 1; break;
            case 'R': RotationAngle1 -= 1; break;
            case 'T': RotationAngle2 += 1; break;
            case 'Y': RotationAngle2 -= 1; break;
            case 'U': RotationAngle3 += 1; break;
            case 'I': RotationAngle3 -= 1; break;
        }
}

void phonglighting_app::onMouseWheel(int pos)
{
    if (pos > 0)
        CameraZPosition += 1;
    else
        CameraZPosition -= 1;
}

bool Pressed = false;

void phonglighting_app::onMouseButton(int button, int action)
{
    if (button == 0)
    {
        ArcBallPrevRotationMatrix = ArcBallActualRotationMatrix;
        ArcBall->click(&MousePosition);
        Pressed = true;
    }
    else
        Pressed = false;
}

void phonglighting_app::onMouseMove(int x, int y)
{    
    MousePosition.s.X = x;
    MousePosition.s.Y = y;

    if (Pressed == true)
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

void phonglighting_app::InitArcBall()
{
    Matrix3fSetIdentity(&ArcBallPrevRotationMatrix);
    Matrix3fSetIdentity(&ArcBallActualRotationMatrix);

    ArcBall = make_unique<ArcBallT>(ArcBallT(640.0f, 480.0f));
    ArcBall->setBounds(static_cast<float>(info.windowWidth), static_cast<float>(info.windowHeight));

    RotationMatrix = vmath::rotate(0.0f, 0.0f, 0.0f);
}

void phonglighting_app::onResize(int w, int h)
{
    info.windowWidth = w;
    info.windowHeight = h;
    ArcBall->setBounds(static_cast<float>(info.windowWidth), static_cast<float>(info.windowHeight));
}

//DECLARE_MAIN(phonglighting_app)

sb7::application *app = 0;
int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    phonglighting_app *app = new phonglighting_app;
    app->run(app);
    delete app;
    return 0;
}