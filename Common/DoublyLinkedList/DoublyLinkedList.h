
#ifndef DOUBLY_LINKED_LIST
#define DOUBLY_LINKED_LIST

using UniqueIdInt = std::uint32_t;

class Particle;
class CellEngineVoxelSimulationSpace;

const UniqueIdInt NoParticleIndex = 0;

typedef Particle& (CellEngineVoxelSimulationSpace::*GetFromPointerType)(UniqueIdInt ParticleIndex);

template <class T>
struct DoublyLinkedListNode
{
public:
    T Prev = NoParticleIndex;
    T Next = NoParticleIndex;
public:
    static GetFromPointerType GetFromPointer;
    static CellEngineVoxelSimulationSpace* ObjectForGettingPointer;
public:
    static void SetGetFromPointer(CellEngineVoxelSimulationSpace* ObjectForGettingPointerParam, GetFromPointerType GetFromPointerParam)
    {
        ObjectForGettingPointer = ObjectForGettingPointerParam;
        GetFromPointer = GetFromPointerParam;
    }
public:
    static void InsertAtFront(T Front, T InsertedParticle);
    static void InsertAtEnd(T End, T InsertedParticle);
    static void InsertBeforeGivenNode(T ParticleInList, T InsertedParticle);
    static void InsertAfterGivenNode(T ParticleInList, T InsertedParticle);
    static void DeleteNode(T Front, T DeletedParticle);
    static T Search(T Front, int IndexToFind);
};

template <class T>
void DoublyLinkedListNode<T>::InsertAtFront(T Front, T InsertedParticle)
{
    (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = NoParticleIndex;

    (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next = Front;

    if (Front != NoParticleIndex)
        (ObjectForGettingPointer->*GetFromPointer)(Front).Prev = InsertedParticle;

    Front = InsertedParticle;
}

template <class T>
void DoublyLinkedListNode<T>::InsertAtEnd(T End, T InsertedParticle)
{
    (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next = NoParticleIndex;

    if (End == NoParticleIndex)
        (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = NoParticleIndex;
    else
        (ObjectForGettingPointer->*GetFromPointer)(End).Next = InsertedParticle;

    (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = End;
}

template <class T>
void DoublyLinkedListNode<T>::InsertAfterGivenNode(T ParticleInList, T InsertedParticle)
{
    if (ParticleInList != NoParticleIndex)
    {
        (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next = (ObjectForGettingPointer->*GetFromPointer)(ParticleInList).Next;

        (ObjectForGettingPointer->*GetFromPointer)(ParticleInList).Next = InsertedParticle;

        (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = ParticleInList;

        if ((ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next != NoParticleIndex)
            (ObjectForGettingPointer->*GetFromPointer)((ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next).Prev = InsertedParticle;
    }
    else
        InsertAtEnd(NoParticleIndex, InsertedParticle);
}

template <class T>
void DoublyLinkedListNode<T>::InsertBeforeGivenNode(T ParticleInList, T InsertedParticle)
{
    if (ParticleInList != NoParticleIndex)
    {
        (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = (ObjectForGettingPointer->*GetFromPointer)(ParticleInList).Prev;

        (ObjectForGettingPointer->*GetFromPointer)(ParticleInList).Prev = InsertedParticle;

        (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Next = ParticleInList;

        if ((ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev != NoParticleIndex)
            (ObjectForGettingPointer->*GetFromPointer)((ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev).Next = InsertedParticle;
        else
            (ObjectForGettingPointer->*GetFromPointer)(InsertedParticle).Prev = NoParticleIndex;
    }
    else
        InsertAtFront(NoParticleIndex, InsertedParticle);
}

template <class T>
void DoublyLinkedListNode<T>::DeleteNode(T Front, T DeletedParticle)
{
    if (DeletedParticle == Front)
        Front = (ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Next;

    (ObjectForGettingPointer->*GetFromPointer)((ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Prev).Next = (ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Next;

    if ((ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Next == NoParticleIndex)
        return;

    (ObjectForGettingPointer->*GetFromPointer)((ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Next).Prev = (ObjectForGettingPointer->*GetFromPointer)(DeletedParticle).Prev;
}

template <class T>
T DoublyLinkedListNode<T>::Search(T Front, int IndexToFind)
{
    T Temp = Front;

    int Index = 0;

    while (Temp != NoParticleIndex)
    {
        if (Index == IndexToFind)
            return Temp;

        Temp = (ObjectForGettingPointer->*GetFromPointer)(Temp).Next;

        Index++;
    }

    return NoParticleIndex;
}

#endif