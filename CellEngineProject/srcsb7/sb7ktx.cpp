
#include "sb7ktx.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#endif 

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "GL/gl3w.h"

namespace sb7
{
    namespace ktx
    {
        namespace File
        {
            static const unsigned char Identifier[] =
            {
                0xAB, 0x4B, 0x54, 0x58, 0x20, 0x31, 0x31, 0xBB, 0x0D, 0x0A, 0x1A, 0x0A
            };

            static const unsigned int Swap32(const unsigned int u32)
            {
                union
                {
                    unsigned int u32;
                    unsigned char u8[4];
                } 
                a, b;

                a.u32 = u32;
                b.u8[0] = a.u8[3];
                b.u8[1] = a.u8[2];
                b.u8[2] = a.u8[1];
                b.u8[3] = a.u8[0];

                return b.u32;
            }

            static const unsigned short Swap16(const unsigned short u16)
            {
                union
                {
                    unsigned short u16;
                    unsigned char u8[2];
                } 
                a, b;

                a.u16 = u16;
                b.u8[0] = a.u8[1];
                b.u8[1] = a.u8[0];

                return b.u16;
            }

            static unsigned int CalculateStride(const Header& HeaderParameter, unsigned int Width, unsigned int Pad = 4)
            {
                unsigned int Channels = 0;

                switch (HeaderParameter.glBaseInternalFormat)
                {
                    case GL_RED: Channels = 1; break;
                    case GL_RG: Channels = 2; break;
                    case GL_BGR:
                    case GL_RGB: Channels = 3; break;
                    case GL_BGRA:
                    case GL_RGBA: Channels = 4; break;
                }

                unsigned int Stride = HeaderParameter.glTypeSize * Channels * Width;

                Stride = (Stride + (Pad - 1)) & ~(Pad - 1);

                return Stride;
            }

            static unsigned int CalculateFaceSize(const Header& HeaderParameter)
            {
                unsigned int Stride = CalculateStride(HeaderParameter, HeaderParameter.PixelWidth);

                return Stride * HeaderParameter.PixelHeight;
            }

            extern unsigned int Load(const char * FileName, unsigned int Tex)
            {
                FILE* ktxFile;
                GLuint Temporary = 0;
                GLuint ReturnValue = 0;
                Header HeaderObject;
                size_t DataStart, DataEnd;
                unsigned char* Data;
                GLenum Target = GL_NONE;

                ktxFile = fopen(FileName, "rb");

                if (!ktxFile)
                    return 0;

                if (fread(&HeaderObject, sizeof(HeaderObject), 1, ktxFile) != 1)
                    goto FailReadLabel;

                if (memcmp(HeaderObject.Identifier, Identifier, sizeof(Identifier)) != 0)
                    goto FailHeaderLabel;

                if (HeaderObject.Endianness == 0x04030201)
                {
                }
                else 
                if (HeaderObject.Endianness == 0x01020304)
                {
                    HeaderObject.Endianness = Swap32(HeaderObject.Endianness);
                    HeaderObject.glType = Swap32(HeaderObject.glType);
                    HeaderObject.glTypeSize = Swap32(HeaderObject.glTypeSize);
                    HeaderObject.glFormat = Swap32(HeaderObject.glFormat);
                    HeaderObject.glInternalFormat = Swap32(HeaderObject.glInternalFormat);
                    HeaderObject.glBaseInternalFormat = Swap32(HeaderObject.glBaseInternalFormat);
                    HeaderObject.PixelWidth = Swap32(HeaderObject.PixelWidth);
                    HeaderObject.PixelHeight = Swap32(HeaderObject.PixelHeight);
                    HeaderObject.PixelDepth = Swap32(HeaderObject.PixelDepth);
                    HeaderObject.ArrayElements = Swap32(HeaderObject.ArrayElements);
                    HeaderObject.Faces = Swap32(HeaderObject.Faces);
                    HeaderObject.mipLevels = Swap32(HeaderObject.mipLevels);
                    HeaderObject.KeyPairBytes = Swap32(HeaderObject.KeyPairBytes);
                }
                else
                    goto FailHeaderLabel;

                if (HeaderObject.PixelHeight == 0)
                {
                    if (HeaderObject.ArrayElements == 0)
                        Target = GL_TEXTURE_1D;
                    else
                        Target = GL_TEXTURE_1D_ARRAY;
                }
                else 
                if (HeaderObject.PixelDepth == 0)
                {
                    if (HeaderObject.ArrayElements == 0)
                    {
                        if (HeaderObject.Faces == 0)
                            Target = GL_TEXTURE_2D;
                        else
                            Target = GL_TEXTURE_CUBE_MAP;
                    }
                    else
                    {
                        if (HeaderObject.Faces == 0)
                            Target = GL_TEXTURE_2D_ARRAY;
                        else
                            Target = GL_TEXTURE_CUBE_MAP_ARRAY;
                    }
                }
                else
                    Target = GL_TEXTURE_3D;

                if (Target == GL_NONE || (HeaderObject.PixelWidth == 0) || (HeaderObject.PixelHeight == 0 && HeaderObject.PixelDepth != 0))
                    goto FailHeaderLabel;

                Temporary = Tex;
                if (Tex == 0)
                    glGenTextures(1, &Tex);

                glBindTexture(Target, Tex);

                DataStart = ftell(ktxFile) + HeaderObject.KeyPairBytes;
                fseek(ktxFile, 0, SEEK_END);
                DataEnd = ftell(ktxFile);
                fseek(ktxFile, DataStart, SEEK_SET);

                Data = new unsigned char [DataEnd - DataStart];
                memset(Data, 0, DataEnd - DataStart);

                fread(Data, 1, DataEnd - DataStart, ktxFile);

                if (HeaderObject.mipLevels == 0)
                    HeaderObject.mipLevels = 1;

                switch (Target)
                {
                    case GL_TEXTURE_1D:
                        glTexStorage1D(GL_TEXTURE_1D, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth);
                        glTexSubImage1D(GL_TEXTURE_1D, 0, 0, HeaderObject.PixelWidth, HeaderObject.glFormat, HeaderObject.glInternalFormat, Data);
                        break;
                    case GL_TEXTURE_2D:
                        // glTexImage2D(GL_TEXTURE_2D, 0, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight, 0, HeaderObject.glFormat, HeaderObject.glType, Data);
                        if (HeaderObject.glType == GL_NONE)
                            glCompressedTexImage2D(GL_TEXTURE_2D, 0, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight, 0, 420 * 380 / 2, Data);
                        else
                        {
                            glTexStorage2D(GL_TEXTURE_2D, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight);
                            {
                                unsigned char* Ptr = Data;
                                unsigned int Height = HeaderObject.PixelHeight;
                                unsigned int Width = HeaderObject.PixelWidth;
                                glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
                                for (unsigned int mipLevel = 0; mipLevel < HeaderObject.mipLevels; mipLevel++)
                                {
                                    glTexSubImage2D(GL_TEXTURE_2D, mipLevel, 0, 0, Width, Height, HeaderObject.glFormat, HeaderObject.glType, Ptr);
                                    Ptr += Height * CalculateStride(HeaderObject, Width, 1);
                                    Height >>= 1;
                                    Width >>= 1;
                                    if (!Height)
                                        Height = 1;
                                    if (!Width)
                                        Width = 1;
                                }
                            }
                        }
                        break;
                    case GL_TEXTURE_3D:
                        glTexStorage3D(GL_TEXTURE_3D, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.PixelDepth);
                        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.PixelDepth, HeaderObject.glFormat, HeaderObject.glType, Data);
                        break;
                    case GL_TEXTURE_1D_ARRAY:
                        glTexStorage2D(GL_TEXTURE_1D_ARRAY, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.ArrayElements);
                        glTexSubImage2D(GL_TEXTURE_1D_ARRAY, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.ArrayElements, HeaderObject.glFormat, HeaderObject.glType, Data);
                        break;
                    case GL_TEXTURE_2D_ARRAY:
                        glTexStorage3D(GL_TEXTURE_2D_ARRAY, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.ArrayElements);
                        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.ArrayElements, HeaderObject.glFormat, HeaderObject.glType, Data);
                        break;
                    case GL_TEXTURE_CUBE_MAP:
                        glTexStorage2D(GL_TEXTURE_CUBE_MAP, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight);
                        // glTexSubImage3D(GL_TEXTURE_CUBE_MAP, 0, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.Faces, HeaderObject.glFormat, HeaderObject.glType, Data);
                        {
                            unsigned int face_size = CalculateFaceSize(HeaderObject);
                            for (unsigned int i = 0; i < HeaderObject.Faces; i++)
                                glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.glFormat, HeaderObject.glType, Data + face_size * i);
                        }
                        break;
                    case GL_TEXTURE_CUBE_MAP_ARRAY:
                        glTexStorage3D(GL_TEXTURE_CUBE_MAP_ARRAY, HeaderObject.mipLevels, HeaderObject.glInternalFormat, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.ArrayElements);
                        glTexSubImage3D(GL_TEXTURE_CUBE_MAP_ARRAY, 0, 0, 0, 0, HeaderObject.PixelWidth, HeaderObject.PixelHeight, HeaderObject.Faces * HeaderObject.ArrayElements, HeaderObject.glFormat, HeaderObject.glType, Data);
                        break;
                    default:
                        goto fail_target;
                }

                if (HeaderObject.mipLevels == 1)
                    glGenerateMipmap(Target);

                ReturnValue = Tex;

            fail_target:
                delete [] Data;

            FailHeaderLabel:;
            FailReadLabel:;
                fclose(ktxFile);

                return ReturnValue;
            }

            bool Save(const char * FileName, unsigned int Target, unsigned int Tex)
            {
                Header HeaderObject;

                memset(&HeaderObject, 0, sizeof(HeaderObject));
                memcpy(HeaderObject.Identifier, Identifier, sizeof(Identifier));
                HeaderObject.Endianness = 0x04030201;

                glBindTexture(Target, Tex);

                glGetTexLevelParameteriv(Target, 0, GL_TEXTURE_WIDTH, (GLint *)&HeaderObject.PixelWidth);
                glGetTexLevelParameteriv(Target, 0, GL_TEXTURE_HEIGHT, (GLint *)&HeaderObject.PixelHeight);
                glGetTexLevelParameteriv(Target, 0, GL_TEXTURE_DEPTH, (GLint *)&HeaderObject.PixelDepth);

                return true;
            }
        }
    }
}
