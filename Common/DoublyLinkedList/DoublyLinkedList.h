
#ifndef DOUBLY_LINKED_LIST
#define DOUBLY_LINKED_LIST

template <class T, class B>
struct DoublyLinkedListNode
{
public:
    T* Prev = nullptr;
    T* Next = nullptr;
public:
    B PrevTemporary = 0;
    B NextTemporary = 0;
public:
    static void InsertAtFront(T* Front, T* InsertedParticle);
    static void InsertAtEnd(T* End, T* InsertedParticle);
    static void InsertBeforeGivenNode(T* ParticleInList, T* InsertedParticle);
    static void InsertAfterGivenNode(T* ParticleInList, T* InsertedParticle);
    static void DeleteNode(T* Front, T* DeletedParticle);
    static T* Search(T* Front, int IndexToFind);
};

template <class T, class B>
void DoublyLinkedListNode<T, B>::InsertAtFront(T* Front, T* InsertedParticle)
{
    InsertedParticle->Prev = nullptr;

    InsertedParticle->next = Front;

    if (Front != nullptr)
        Front->Prev = InsertedParticle;

    Front = InsertedParticle;
}

template <class T, class B>
void DoublyLinkedListNode<T, B>::InsertAtEnd(T* End, T* InsertedParticle)
{
    InsertedParticle->Next = nullptr;

    if (End == nullptr)
        InsertedParticle->Prev = nullptr;
    else
        End->Next = InsertedParticle;

    InsertedParticle->Prev = End;
}

template <class T, class B>
void DoublyLinkedListNode<T, B>::InsertAfterGivenNode(T* ParticleInList, T* InsertedParticle)
{
    if (ParticleInList != nullptr)
    {
        InsertedParticle->Next = ParticleInList->Next;

        ParticleInList->Next = InsertedParticle;

        InsertedParticle->Prev = ParticleInList;

        if (InsertedParticle->Next != nullptr)
            InsertedParticle->Next->Prev = InsertedParticle;
    }
    else
        InsertAtEnd(nullptr, InsertedParticle);
}

template <class T, class B>
void DoublyLinkedListNode<T, B>::InsertBeforeGivenNode(T* ParticleInList, T* InsertedParticle)
{
    if (ParticleInList != nullptr)
    {
        InsertedParticle->Prev = ParticleInList->Prev;

        ParticleInList->Prev = InsertedParticle;

        InsertedParticle->Next = ParticleInList;

        if(InsertedParticle->Prev != nullptr)
            InsertedParticle->Prev->Next = InsertedParticle;
        else
            InsertedParticle->Prev = nullptr;
    }
    else
        InsertAtFront(nullptr, InsertedParticle);
}

template <class T, class B>
void DoublyLinkedListNode<T, B>::DeleteNode(T* Front, T* DeletedParticle)
{
    if (DeletedParticle == Front)
        Front = DeletedParticle->Next;

    DeletedParticle->Prev->Next = DeletedParticle->Next;

    if (DeletedParticle->Next == nullptr)
        return;

    DeletedParticle->Next->Prev = DeletedParticle->Prev;
}

template <class T, class B>
T* DoublyLinkedListNode<T, B>::Search(T* Front, int IndexToFind)
{
    T* Temp = Front;

    int Index = 0;

    while (Temp != nullptr)
    {
        if (Index == IndexToFind)
            return Temp;

        Temp = Temp->Next;

        Index++;
    }

    return nullptr;
}

#endif