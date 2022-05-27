
#include <sb7textoverlay.h>

#include <sb7ktx.h>

namespace sb7
{
    void TextOverlay::Init(int width, int height, const char* font)
    {
        GLuint VertexShader;
        GLuint FragmentShader;

        BufferWidth = width;
        BufferHeight = height;

        VertexShader = glCreateShader(GL_VERTEX_SHADER);
        FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

        static const char * VertexShaderSource[] =
        {
            "#version 440 core\n"
            "void main(void)\n"
            "{\n"
            "    gl_Position = vec4(float((gl_VertexID >> 1) & 1) * 2.0 - 1.0,\n"
            "                       float((gl_VertexID & 1)) * 2.0 - 1.0,\n"
            "                       0.0, 1.0);\n"
            "}\n"
        };

        static const char * FragmentShaderSource[] =
        {
            "#version 440 core\n"
            "layout (origin_upper_left) in vec4 gl_FragCoord;\n"
            "layout (location = 0) out vec4 o_color;\n"
            "layout (binding = 0) uniform isampler2D TextBuffer;\n"
            "layout (binding = 1) uniform isampler2DArray FontTexture;\n"
            "void main(void)\n"
            "{\n"
            "    ivec2 frag_coord = ivec2(gl_FragCoord.xy);\n"
            "    ivec2 char_size = textureSize(FontTexture, 0).xy;\n"
            "    ivec2 char_location = frag_coord / char_size;\n"
            "    ivec2 texel_coord = frag_coord % char_size;\n"
            "    int character = texelFetch(TextBuffer, char_location, 0).x;\n"
            "    float val = texelFetch(FontTexture, ivec3(texel_coord, character), 0).x;\n"
            "    if (val == 0.0)\n"
            "        discard;\n"
            "    o_color = vec4(1.0);\n"
            "}\n"
        };

        glShaderSource(VertexShader, 1, VertexShaderSource, nullptr);
        glCompileShader(VertexShader);

        glShaderSource(FragmentShader, 1, FragmentShaderSource, nullptr);
        glCompileShader(FragmentShader);

        TextProgram = glCreateProgram();
        glAttachShader(TextProgram, VertexShader);
        glAttachShader(TextProgram, FragmentShader);
        glLinkProgram(TextProgram);

        glDeleteShader(FragmentShader);
        glDeleteShader(VertexShader);

        // glCreateVertexArrays(1, &VAO);
        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);

        // glCreateTextures(GL_TEXTURE_2D, 1, &TextBuffer);
        glGenTextures(1, &TextBuffer);
        glBindTexture(GL_TEXTURE_2D, TextBuffer);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_R8UI, width, height);

        FontTexture = sb7::ktx::File::Load(font ? font : "media/textures/cp437_9x16.ktx");

        ScreenBuffer = new char[width * height];
        memset(ScreenBuffer, 0, width * height);
    }

    void TextOverlay::TearDown()
    {
        delete[] ScreenBuffer;
        glDeleteTextures(1, &FontTexture);
        glDeleteTextures(1, &TextBuffer);
        glDeleteVertexArrays(1, &VAO);
        glDeleteProgram(TextProgram);
    }

    void TextOverlay::Draw()
    {
        glUseProgram(TextProgram);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, TextBuffer);
        if (Dirty)
        {
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, BufferWidth, BufferHeight, GL_RED_INTEGER, GL_UNSIGNED_BYTE, ScreenBuffer);
            Dirty = false;
        }
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D_ARRAY, FontTexture);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }

    void TextOverlay::DrawText(const char* StringParameter, int X, int Y)
    {
        char* Destination = ScreenBuffer + Y * BufferWidth + X;
        strcpy(Destination, StringParameter);
        Dirty = true;
    }

    void TextOverlay::Scroll(int Lines)
    {
        const char* Source = ScreenBuffer + Lines * BufferWidth;
        char * Destination = ScreenBuffer;

        memmove(Destination, Source, (BufferHeight - Lines) * BufferWidth);

        Dirty = true;
    }

    void TextOverlay::Print(const char* StringParameter)
    {
        const char* StringPtr = StringParameter;
        char Character;
        char* Destination = ScreenBuffer + CursorY * BufferWidth + CursorX;

        while (*StringPtr != 0)
        {
            Character = *StringPtr++;
            if (Character == '\n')
            {
                CursorY++;
                CursorX = 0;
                if (CursorY >= BufferHeight)
                {
                    CursorY--;
                    Scroll(1);
                }
                Destination = ScreenBuffer + CursorY * BufferWidth + CursorX;
            }
            else
            {
                *Destination++ = Character;
                CursorX++;
                if (CursorX >= BufferWidth)
                {
                    CursorY++;
                    CursorX = 0;
                    if (CursorY >= BufferHeight)
                    {
                        CursorY--;
                        Scroll(1);
                    }
                    Destination = ScreenBuffer + CursorY * BufferWidth + CursorX;
                }
            }
        }

        Dirty = true;
    }

    void TextOverlay::MoveCursor(int x, int y)
    {
        CursorX = x;
        CursorY = y;
    }

    void TextOverlay::Clear()
    {
        memset(ScreenBuffer, 0, BufferWidth * BufferHeight);
        Dirty = true;
        CursorX = 0;
        CursorY = 0;
    }
}
