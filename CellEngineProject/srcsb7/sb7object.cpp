
#define _CRT_SECURE_NO_WARNINGS 1

#include "GL/gl3w.h"
#include <object.h>

#include <cstdio>

namespace sb7
{
    GraphicObject::GraphicObject() : DataBuffer(0), IndexType(0), VAO(0)
    {
    }

    GraphicObject::~GraphicObject() = default;

    void GraphicObject::Load(const char* FileName)
    {
        FILE* GraphicObjectFile = fopen(FileName, "rb");
        size_t FileSize;
        char* FileData;

        this->Free();

        fseek(GraphicObjectFile, 0, SEEK_END);
        FileSize = ftell(GraphicObjectFile);
        fseek(GraphicObjectFile, 0, SEEK_SET);

        FileData = new char[FileSize];

        fread(FileData, FileSize, 1, GraphicObjectFile);

        char* Ptr = FileData;
        auto Header = (SB6M_HEADER*)Ptr;
        Ptr += Header->Size;

        SB6M_VERTEX_ATTRIB_CHUNK* VertexAttribChunk = nullptr;
        SB6M_CHUNK_VERTEX_DATA* VertexDataChunk = nullptr;
        SB6M_CHUNK_INDEX_DATA* IndexDataChunk = nullptr;
        SB6M_CHUNK_SUB_OBJECT_LIST* SubGraphicObjectChunk = nullptr;
        SB6M_DATA_CHUNK* DataChunk = nullptr;

        for (unsigned int ChunkIndex = 0; ChunkIndex < Header->NumberOfChunks; ChunkIndex++)
        {
            auto Chunk = (SB6M_CHUNK_HEADER*)Ptr;
            Ptr += Chunk->Size;
            switch (Chunk->ChunkType)
            {
                case SB6M_CHUNK_TYPE_VERTEX_ATTRIBS: VertexAttribChunk = (SB6M_VERTEX_ATTRIB_CHUNK*)Chunk; break;
                case SB6M_CHUNK_TYPE_VERTEX_DATA: VertexDataChunk = (SB6M_CHUNK_VERTEX_DATA*)Chunk; break;
                case SB6M_CHUNK_TYPE_INDEX_DATA: IndexDataChunk = (SB6M_CHUNK_INDEX_DATA*)Chunk; break;
                case SB6M_CHUNK_TYPE_SUB_OBJECT_LIST: SubGraphicObjectChunk = (SB6M_CHUNK_SUB_OBJECT_LIST*)Chunk; break;
                case SB6M_CHUNK_TYPE_DATA: DataChunk = (SB6M_DATA_CHUNK*)Chunk; break;
                default: break;
            }
        }

        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);

        if (DataChunk != nullptr)
        {
            glGenBuffers(1, &DataBuffer);
            glBindBuffer(GL_ARRAY_BUFFER, DataBuffer);
            glBufferData(GL_ARRAY_BUFFER, DataChunk->DataLength, (unsigned char*)DataChunk + DataChunk->DataOffset, GL_STATIC_DRAW);
        }
        else
        {
            unsigned int DataSize = 0;
            unsigned int SizeUsed = 0;

            if (VertexDataChunk != nullptr)
                DataSize += VertexDataChunk->DataSize;

            if (IndexDataChunk != nullptr)
                DataSize += IndexDataChunk->IndexCount * (IndexDataChunk->IndexType == GL_UNSIGNED_SHORT ? sizeof(GLushort) : sizeof(GLubyte));

            glGenBuffers(1, &DataBuffer);
            glBindBuffer(GL_ARRAY_BUFFER, DataBuffer);
            glBufferData(GL_ARRAY_BUFFER, DataSize, nullptr, GL_STATIC_DRAW);

            if (VertexDataChunk != nullptr)
            {
                glBufferSubData(GL_ARRAY_BUFFER, 0, VertexDataChunk->DataSize, FileData + VertexDataChunk->DataOffset);
                SizeUsed += VertexDataChunk->DataOffset;
            }

            if (IndexDataChunk != nullptr)
                glBufferSubData(GL_ARRAY_BUFFER, SizeUsed, IndexDataChunk->IndexCount * (IndexDataChunk->IndexType == GL_UNSIGNED_SHORT ? sizeof(GLushort) : sizeof(GLubyte)), FileData + IndexDataChunk->IndexDataOffset);
        }

        for (unsigned int attrib_index = 0; attrib_index < VertexAttribChunk->AttribCount; attrib_index++)
        {
            SB6M_VERTEX_ATTRIB_DECL &attrib_decl = VertexAttribChunk->AttribData[attrib_index];
            glVertexAttribPointer(attrib_index, attrib_decl.Size, attrib_decl.Type, attrib_decl.Flags & SB6M_VERTEX_ATTRIB_FLAG_NORMALIZED ? GL_TRUE : GL_FALSE, attrib_decl.Stride, (GLvoid*)(uintptr_t)attrib_decl.DataOffset);
            glEnableVertexAttribArray(attrib_index);
        }

        if (IndexDataChunk != nullptr)
        {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, DataBuffer);
            IndexType = IndexDataChunk->IndexType;
        }
        else
            IndexType = GL_NONE;

        if (SubGraphicObjectChunk != nullptr)
        {
            if (SubGraphicObjectChunk->Count > MAX_SUB_GRAPHIC_OBJECTS)
                SubGraphicObjectChunk->Count = MAX_SUB_GRAPHIC_OBJECTS;

            for (unsigned int chunk_index = 0; chunk_index < SubGraphicObjectChunk->Count; chunk_index++)
                SubGraphicObjects[chunk_index] = SubGraphicObjectChunk->SubObject[chunk_index];

            NumberOfSubGraphicObjects = SubGraphicObjectChunk->Count;
        }
        else
        {
            SubGraphicObjects[0].First = 0;
            SubGraphicObjects[0].Count = IndexType != GL_NONE ? IndexDataChunk->IndexCount : VertexDataChunk->TotalVertices;
            NumberOfSubGraphicObjects = 1;
        }

        delete[] FileData;

        fclose(GraphicObjectFile);

        glBindVertexArray(0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    void GraphicObject::Free()
    {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &DataBuffer);

        VAO = 0;
        DataBuffer = 0;
    }

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wint-to-pointer-cast"
    void GraphicObject::RenderSubGraphicObject(unsigned int GraphicObjectIndex, unsigned int InstanceCount, unsigned int BaseInstance)
    {
        glBindVertexArray(VAO);

        if (IndexType != GL_NONE)
            glDrawElementsInstancedBaseInstance(GL_TRIANGLES, SubGraphicObjects[GraphicObjectIndex].Count, IndexType, (void*)SubGraphicObjects[GraphicObjectIndex].First, InstanceCount, BaseInstance);
        else
            glDrawArraysInstancedBaseInstance(GL_TRIANGLES, SubGraphicObjects[GraphicObjectIndex].First, SubGraphicObjects[GraphicObjectIndex].Count, InstanceCount, BaseInstance);
    }
    #pragma GCC diagnostic pop
}
