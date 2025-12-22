#version 450

bool ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside = false;
bool ShowDetailsInAtomScale = true;
bool CheckAtomVisibility = true;
float CutZ = 200;
float Distance = 2300;
float ViewPositionZ = 50;
float CameraXPosition = 0;
float CameraYPosition = 0;
float CameraZPosition = 0;
float XLowToDrawInAtomScale = -200;
float XHighToDrawInAtomScale = 200;
float YLowToDrawInAtomScale = -200;
float YHighToDrawInAtomScale = 200;
float ZLowToDrawInAtomScale = -650;
float ZHighToDrawInAtomScale = -50;

uniform vec3 Center;

layout(local_size_x = 256) in;

struct GPUAtom
{
    float X;
    float Y;
    float Z;
    short ColorR;
    short ColorG;
    short ColorB;
    uint _padding1[7];
};

struct GPUParticle
{
    uint EntityId;
    uint ChainId;
    uint Index;
    uint AtomOffset;
    uint AtomCount;
    uint _padding[3];
};

layout(std430, binding = 0) readonly buffer ParticleBuffer
{
    GPUParticle GPUParticles[];
};

layout(std430, binding = 1) readonly buffer AtomBuffer
{
    GPUAtom GPUAtoms[];
};

layout(std430, binding = 1) buffer VisibleIndices
{
    uint VisibleList[];
};

layout(std430, binding = 2) buffer IndirectDraw
{
    uint instanceCount;
}
drawCmd;




//OUT

struct UniformsBlock
{
    mat4 MoveMatrix;
    vec3 Color;
    float padding;
};

layout(std430, binding = 0) buffer ParticleBuffer
{
    UniformsBlock ParticlesOut[];
};

bool CheckDistanceToDrawDetailsInAtomScale(const float XNew, const float YNew, const float ZNew)
{
    if (CheckAtomVisibility == true)
    {
        if (ViewPositionZ > Distance)
        {
            if (ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside == false)
                return ZNew > CutZ && sqrt((XNew * XNew) + (YNew * YNew) + (ZNew * ZNew)) > Distance;
            else
                return true;
        }
        else
            return (ZNew > ViewPositionZ + ZLowToDrawInAtomScale && ZNew < ViewPositionZ + ZHighToDrawInAtomScale && XNew > XLowToDrawInAtomScale && XNew < XHighToDrawInAtomScale && YNew > YLowToDrawInAtomScale && YNew < YHighToDrawInAtomScale);
    }
    else
        return false;
}

bool GetFinalVisibilityInModelWorld(const vec3& AtomPosition, UniformsBlock* MatrixUniformBlock, const bool CountNewPosition, const bool DrawOutsideBorder) const
{
    if (CountNewPosition == true)
    {
        const float XNew = MatrixUniformBlock->MoveMatrix[0][0] * (AtomPosition.x + CameraXPosition - Center.x) + MatrixUniformBlock->MoveMatrix[1][0] * (AtomPosition.y + CameraYPosition - Center.y) + MatrixUniformBlock->MoveMatrix[2][0] * (AtomPosition.z + CameraZPosition - Center.z);
        const float YNew = MatrixUniformBlock->MoveMatrix[0][1] * (AtomPosition.x + CameraXPosition - Center.x) + MatrixUniformBlock->MoveMatrix[1][1] * (AtomPosition.y + CameraYPosition - Center.y) + MatrixUniformBlock->MoveMatrix[2][1] * (AtomPosition.z + CameraZPosition - Center.z);
        const float ZNew = MatrixUniformBlock->MoveMatrix[0][2] * (AtomPosition.x + CameraXPosition - Center.x) + MatrixUniformBlock->MoveMatrix[1][2] * (AtomPosition.y + CameraYPosition - Center.y) + MatrixUniformBlock->MoveMatrix[2][2] * (AtomPosition.z + CameraZPosition - Center.z);

        if (DrawOutsideBorder == true)
    	    if (CheckDistanceToDrawDetailsInAtomScale(XNew, YNew, ZNew) == true)
                return true;

    	return false;
    }
}

bool CreateUniformBlockForVertexShader(const vec3& Position, const vec3& Color, const mat4& ViewMatrix, mat4 ModelMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, bool DrawAdditional, const UnsignedInt ParticleSectorXIndex)
{
    bool FinalVisibilityInModelWorld = false;

    UniformsBlock MatrixUniformBlock;
    MatrixUniformBlock.MoveMatrix = ViewMatrix * ModelMatrix;
    MatrixUniformBlock.Color = Color;

//CZY WPISUJE DO TABELI UniformsBlocks bo to problem bo musza byc rownolegle czy wpisuje VisibleIndices i tylko wtedy przekazuje dalej
//ParticlesOut musi byc tablica zalezna od zmiennej wykonanej przez gl_GlobalInvocationID.x
//ale jak to skompresowac

    //UniformsBlocks.emplace_back(MatrixUniformBlock)
    ParticlesOut.emplace_back(MatrixUniformBlock);

    if (DrawAdditional == true)
    	FinalVisibilityInModelWorld = GetFinalVisibilityInModelWorld(Position, &MatrixUniformBlock, CountNewPosition, DrawOutsideBorder);

    return FinalVisibilityInModelWorld;
}

bool RenderObject(const GPUAtom& AtomObject, const Particle& ParticleObject, const mat4& ViewMatrix, const bool CountNewPosition, const bool DrawCenter, const bool DrawOutsideBorder, const bool RenderObjectParameter)
{
    bool FinalVisibilityInModelWorld{};

    //if (RenderObjectParameter == true)
    //    NumberOfAllRenderedAtoms++;

    const vec3 AtomPosition = LengthUnit * vec3(AtomObject.X, AtomObject.Y, AtomObject.Z);
    const vec3 SizeLocal = GetSize(AtomObject);
    const mat4 ModelMatrix = translate(AtomPosition.x - CameraXPosition - Center.x, AtomPosition.y + CameraYPosition - Center.y, AtomPosition.z + CameraZPosition - Center.z));

    FinalVisibilityInModelWorld = CreateUniformBlockForVertexShader(AtomPosition, vec3(AtomObject.ColorR, AtomObject.ColorG, AtomObject.ColorB), ViewMatrix, ModelMatrix, CountNewPosition, DrawCenter, DrawOutsideBorder, true);

    return FinalVisibilityInModelWorld;
}





void main()
{
    uint particleId = gl_GlobalInvocationID.x;
    if (particleId >= particles.length())
       return;

    GPUParticle particle = particles[particleId];

    for (uint i = 0; i < particle.AtomCount; i++)
    {
        GPUAtom atom = atoms[particle.AtomOffset + i];

        RenderObject(atom, particle, ViewMatrix, false, false, false, true);
    }






//    uint id = gl_GlobalInvocationID.x;
//    if (id >= particles.length())
//        return;

//    Particle p = particles[id];
    vec3 pos = p.posRadius.xyz;
    float radius = p.posRadius.w;

    bool visible = true;
    for (int i = 0; i < 6; i++)
    {
        float dist = dot(vec4(pos, 1.0), frustumPlanes[i]);
        if (dist < -radius)
        {
            visible = false;
            break;
        }
    }

    float distToView = distance(pos, viewPos);
    if (distToView > maxRenderDistance)
    {
        visible = false;
    }

    // LOD selection based on distance
    uint lod = 0;
    if (distToView > 100.0)
        lod = 1;
    if (distToView > 500.0)
        lod = 2;
    if (distToView > 1000.0)
        lod = 3;

    particles[id].visible = visible ? (lod + 1) : 0;

    if (visible)
    {
        uint index = atomicAdd(drawCmd.instanceCount, 1);
        visibleList[index] = id;
    }
}
