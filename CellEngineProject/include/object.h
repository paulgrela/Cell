
#ifndef OBJECT_H
#define OBJECT_H

#include "sb6mfile.h"

#ifndef SB6M_FILETYPES_ONLY

#include <GL/glcorearb.h>

namespace sb7
{
    class GraphicObject
    {
    public:
        GraphicObject();
        ~GraphicObject();

        inline void Render(unsigned int InstanceCount = 1, unsigned int BaseInstance = 0)
        {
            RenderSubGraphicObject(0, InstanceCount, BaseInstance);
        }

        void RenderSubGraphicObject(unsigned int GraphicObjectIndex, unsigned int InstanceCount = 1, unsigned int BaseInstance = 0);

        void GetSubGraphicObjectInfo(unsigned int Index, GLuint &First, GLuint &NumberOfGraphicObjectsParameter)
        {
            if (Index >= NumberOfSubGraphicObjects)
            {
                First = 0;
                NumberOfGraphicObjectsParameter = 0;
            }
            else
            {
                First = SubGraphicObjects[Index].First;
                NumberOfGraphicObjectsParameter = SubGraphicObjects[Index].Count;
            }
        }

        [[nodiscard]] unsigned int GetNumberOfSubGraphicObject() const
        { 
            return NumberOfSubGraphicObjects;
        }

        [[nodiscard]] GLuint GetVAO() const
        { 
            return VAO;
        }

        void Load(const char* FileName);
        void Free();

    private:
        GLuint DataBuffer;
        GLuint VAO;
        GLuint IndexType;

        enum { MAX_SUB_GRAPHIC_OBJECTS = 256 };

        unsigned int NumberOfSubGraphicObjects;
        SB6M_SUB_OBJECT_DECL SubGraphicObjects[MAX_SUB_GRAPHIC_OBJECTS];
    };
}

#endif

#endif

