
#include <memory>
#include <exception>

#include "CellEnginePDBWindowGL.h"

using namespace std;

#pragma region PDB
PDBWindowGL::PDBWindowGL() : WindowGL(), PDBDataFileObjectPointer(NULL), ShowBonds(true), ShowPDBSize(1.00f), RefreshListOfDrawing(false), ProjectionType(false), ChosenAtomIndex(-1)
{
	char* FileName = NULL;
	
	if (__argc > 1)
		FileName = __argv[1];
	else
	{
		MessageBox(NULL, "Brak nazwy pliku w parametrze linii komend", "Przegladarka PDB", MB_OK | MB_ICONWARNING);
		PostQuitMessage(0);
		return;
	}

	if (!OpenPDBFile(FileName))
		PostQuitMessage(0);
}

PDBWindowGL::~PDBWindowGL()
{
	delete PDBDataFileObjectPointer;
	PDBDataFileObjectPointer = NULL;
}

bool PDBWindowGL::OpenPDBFile(const char* FileName)
{
	try
	{
		PDBDataFileObjectPointer = new PDBDataFile(FileName);
	}
	catch (const std::exception& exc)
	{
		char komunikat[256] = "Error with loading PDB object:\n";
		strcat(komunikat, exc.what());
		MessageBox(HandleWindow, komunikat, "PDB VIEWER", MB_OK | MB_ICONERROR);
		return false;
	}

	return true;
}


void PDBWindowGL::DrawAtoms(PDBDataFile* PDBDataFileObjectPtr, const double LengthUnit, const double AtomSizeLengthUnit, bool MakeColors) const
{
	if (PDBDataFileObjectPtr == nullptr)
		return;

	glPushMatrix();

	//zero ukladu wspolrzednych w srodku masy	
	DoubleVectorType MassCenter = LengthUnit * PDBDataFileObjectPtr->MassCenter();
	glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

	GLUquadricObj* Quadriga = gluNewQuadric();
	gluQuadricDrawStyle(Quadriga, GLU_FILL);
	glShadeModel(GL_SMOOTH);

	const IntType HowManyPointsInEachDimension = 10; //jakosc sfery (moze byc zalezna np. od tego, czy trwa przeciaganie)

	glInitNames(); //tworzenie stosu nazw
	glPushName(-1); //kladzenie "pustego" elementu (nazwy) na stosie

	for (const unique_ptr<Atom>& AtomObjectPtr : PDBDataFileObjectPtr->GetAtoms())
	{
		DoubleVectorType AtomPosition = LengthUnit * AtomObjectPtr->Position();

		//tu jest dobre miejsce na ewentalne filtrowanie rysowanych atomow

		if (MakeColors == true)
		{
			char AtomSymbol[3];
			Atom::GetAtomSymbol(AtomObjectPtr->Name, AtomSymbol);
			ChooseAtomColor(AtomSymbol, 1.0f);
		}

		glPushMatrix();
		glTranslated(AtomPosition.X, AtomPosition.Y, AtomPosition.Z);
		glLoadName(AtomObjectPtr->AtomIndex); //zastepowanie nazwy z wierzchu stosu
		gluSphere(Quadriga, LengthUnit * AtomSizeLengthUnit, HowManyPointsInEachDimension, HowManyPointsInEachDimension);
		glPopMatrix();
	}

	glPopName(); //usuwanie elementu ze stosu

	gluDeleteQuadric(Quadriga);
	glPopMatrix();
}


void PDBWindowGL::DrawChosenAtom(PDBDataFile* PDBDataFileObject, const double LengthUnit, const double AtomSizeLengthUnit, IntType ChosenAtomIndex, bool UseGrid) const
{
	if (ChosenAtomIndex < 0)
		return;

	glPushMatrix();

	GLUquadricObj* Quadric = gluNewQuadric();
	if (UseGrid)
		gluQuadricDrawStyle(Quadric, GLU_LINE);
	glLineWidth(1.0f);
	glShadeModel(GL_SMOOTH);
	DoubleVectorType polozenieWskazanegoAtomu = LengthUnit * (PDBDataFileObject->GetAtom(ChosenAtomIndex)->Position() - PDBDataFileObject->MassCenter());
	glTranslated(polozenieWskazanegoAtomu.X, polozenieWskazanegoAtomu.Y, polozenieWskazanegoAtomu.Z);
	gluSphere(Quadric, LengthUnit * AtomSizeLengthUnit, 15, 15);
	gluDeleteQuadric(Quadric);

	glPopMatrix();
}

							
							void Print(char* Text, IntType CharactersNumber, UnsignedIntType Font, IntType FirstCharCode)
							{
								if (Text == NULL || Text == "")
									return;

								//MessageBox(HandleWindow, Text, "PDB VIEWER", MB_OK);

								glPushAttrib(GL_LIST_BIT);   //Odklada na stos atrybuty wyswietlania
								//glListBase(Font - FirstCharCode);    //Ustawia podstawe znakow
								glListBase(0);    //Ustawia podstawe znakow
								glCallLists(CharactersNumber, GL_UNSIGNED_BYTE, Text);
								//Wyswietla kolejno listy liter (Text)
								glPopAttrib();               //Przywraca ze stosu atrybuty wyswietlania
							}

void PDBWindowGL::DrawChosenAtomDescription(IntType ChosenAtomIndex, IntType BitmapFont)
{
	char ChosenAtomDescription[256] = "";
	_itoa(PDBDataFileObjectPointer->GetAtom(ChosenAtomIndex)->Serial, ChosenAtomDescription, 10);
	strcat(ChosenAtomDescription, ". ");
	strcat(ChosenAtomDescription, PDBDataFileObjectPointer->GetAtom(ChosenAtomIndex)->Name);
	strcat(ChosenAtomDescription, " (");
	strcat(ChosenAtomDescription, PDBDataFileObjectPointer->GetAtom(ChosenAtomIndex)->ResName);
	strcat(ChosenAtomDescription, ")");

	//MessageBox(HandleWindow, ChosenAtomDescription, "PDB VIEWER", MB_OK);
	
							glPushMatrix();
							glLoadIdentity();
							//napis umieszczony w lewym dolnym rogu			
							float Coordinates = UserAreaHeight / (float)UserAreaWidth;
							float TextPosition = (!ProjectionType) ? -0.095f : -2.85f; //porownaj z parametrami frustum (argumenty wywolania glFrustum/glOrtho)

							glRasterPos3f(TextPosition, TextPosition * Coordinates, -0.3f);
							Print(ChosenAtomDescription, 256, BitmapFont, 32);

							glPopMatrix();
	
}


void PDBWindowGL::ChooseAtomColor(const char* AtomSymbol, const float Alpha = 1.0f) const
{
	switch (AtomSymbol[0]) //na razie tylko jednoliterowe
	{
		case 'C': glColor4f(0.25f, 0.75f, 0.75f, Alpha); break;
		case 'O': glColor4f(1.00f, 0.00f, 0.00f, Alpha); break;
		case 'H': glColor4f(1.00f, 1.00f, 1.00f, Alpha); break;
		case 'N': glColor4f(0.00f, 0.00f, 1.00f, Alpha); break;
		case 'P': glColor4f(0.50f, 0.50f, 0.20f, Alpha); break;
		default: glColor4f(0.50f, 0.50f, 0.50f, Alpha); break;
	}
}

void PDBWindowGL::DrawBonds(PDBDataFile* PDBDataFileObject, const double LengthUnit, bool MakeColors) const
{
	if (PDBDataFileObject == NULL)
		return;

	glPushMatrix();

	//zero ukladu wspolrzednych w srodku masy
	DoubleVectorType MassCenter = LengthUnit * PDBDataFileObject->MassCenter();;
	glTranslated(-MassCenter.X, -MassCenter.Y, -MassCenter.Z);

	glColor3f(0, 0, 0);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

	glShadeModel(GL_SMOOTH);
	glLineWidth(2);
	glBegin(GL_LINES);

	for (IntType AtomIndex1 = 0; AtomIndex1 < PDBDataFileObject->GetNumberOfAtoms(); AtomIndex1++)
		for (IntType AtomIndex2 = AtomIndex1 + 1; AtomIndex2 < PDBDataFileObject->GetNumberOfAtoms(); AtomIndex2++)
		{
			DoubleVectorType Atom1Position = LengthUnit * PDBDataFileObject->GetAtom(AtomIndex1)->Position();
			DoubleVectorType Atom2Position = LengthUnit * PDBDataFileObject->GetAtom(AtomIndex2)->Position();
			DoubleVectorType Position12 = Atom2Position - Atom1Position;

			if (Position12.Length() < 0.17) //w AA
			{
				if (MakeColors)
				{
					char AtomSymbol[3];
					Atom::GetAtomSymbol(PDBDataFileObject->GetAtom(AtomIndex1)->Name, AtomSymbol);
					ChooseAtomColor(AtomSymbol, 1.0f);
				}
				glVertex3d(Atom1Position.X, Atom1Position.Y, Atom1Position.Z);
				if (MakeColors)
				{
					char AtomSymbol[3];
					Atom::GetAtomSymbol(PDBDataFileObject->GetAtom(AtomIndex2)->Name, AtomSymbol);
					ChooseAtomColor(AtomSymbol, 1.0f);
				}
				glVertex3d(Atom2Position.X, Atom2Position.Y, Atom2Position.Z);
			}
		}
	glEnd();

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

	glPopMatrix();
}

//lista wyswietlania
UnsignedIntType PDBWindowGL::CreateListOfDrawing(const double LengthUnit) const
{
	GLuint DrawList = glGenLists(1);

	glNewList(DrawList, GL_COMPILE);
	glColor3f(1, 1, 1);

	if (ShowBonds)
		DrawBonds(PDBDataFileObjectPointer, LengthUnit, true);

	if (ShowPDBSize != 0)
		DrawAtoms(PDBDataFileObjectPointer, LengthUnit, ShowPDBSize, true);

	glEndList();

	return DrawList;
}

LRESULT PDBWindowGL::WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	long Result = WindowGL::WndProc(hWnd, message, wParam, lParam);

	switch (message)
	{
		case WM_KEYDOWN:
		{
			bool CallDrawStage = false;
			switch (wParam)
			{
				case VK_F1: MessageBox(HandleWindow, "Shortcut Keys: \ nZ - show / hide bonds \ nX - size of atoms \ nC - projection mode", "PDB VIEWER", MB_OK); break;
				case 'Z': ShowBonds = !ShowBonds; CallDrawStage = true; break;
				case 'X':
				{
					float NewSize = -1.0f;
				
					if (ShowPDBSize == 0.0f)
						NewSize = 0.25f;
					if (ShowPDBSize == 0.25f)
						NewSize = 0.5f;
					if (ShowPDBSize == 0.5f)
						NewSize = 1.0f;
					if (ShowPDBSize == 1.0f)
						NewSize = 0.0f;
					if (NewSize == -1.0f)
						NewSize = 1.0f;

					ShowPDBSize = NewSize;
				}
				CallDrawStage = true;
				break;
			
				case 'C': ProjectionType = !ProjectionType; SetStage(ProjectionType); CallDrawStage = true; break;
			}

			if (CallDrawStage)
			{
				RefreshListOfDrawing = true;
				DrawStage();
			}
			break;
		}

		case WM_LBUTTONDOWN:

			if (ShowPDBSize == 0)
				break;

			POINT CursorMouseIndex;
			CursorMouseIndex.x = LOWORD(lParam);
			CursorMouseIndex.y = HIWORD(lParam);

			IntType PrevIndex = ChosenAtomIndex;
			IntType NextIndex = ChooseAtom(CursorMouseIndex);

			if (NextIndex != -1)
				ChosenAtomIndex = NextIndex;

			if (PrevIndex != -1 || ChosenAtomIndex != -1)
				DrawStage();

			break;
	}

	return Result;
}

//selekcja
IntType PDBWindowGL::ChooseAtom(POINT MouseCursorPosition)
{
	//przygotowanie bufora zaznaczania
	const IntType MarkingBufferSize = 1024;
	unsigned int MarkingBuffer[MarkingBufferSize];
	ZeroMemory(MarkingBuffer, MarkingBufferSize);
	glSelectBuffer(MarkingBufferSize, MarkingBuffer);

	//przygotowanie promienia pod myszka = b. waskie frustum
	//niestety powtorzenie kodu z CWindowGL::SetStage
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	//przez ta funkcje nie mozna po prostu wywolac metody SetStage (mozna ja tam wstrzyknac, ale wowczas zamieszanie w tekscie)
	int Viewport[4];
	glGetIntegerv(GL_VIEWPORT, Viewport);
	gluPickMatrix(MouseCursorPosition.x, UserAreaHeight - MouseCursorPosition.y, 1, 1, Viewport);
	float wsp = UserAreaHeight / (float)UserAreaWidth;
	if (!ProjectionType)
		//left,right,bottom,top,znear,zfar (clipping) 	
		//mnozenie macierzy rzutowania przez macierz perspektywy - ustalanie frustum 	
		glFrustum(-0.1, 0.1, wsp * -0.1, wsp * 0.1, 0.3, 100.0);
	else
		glOrtho(-3, 3, wsp * -3, wsp * 3, 0.3, 100.0); //ZLE DZIALA Z glOrtho

	glMatrixMode(GL_MODELVIEW);

	//przelaczenie w tryb selekcji HitIndex renderowanie sceny
	glRenderMode(GL_SELECT); //zakomentuj to polecenie, zeby zobaczyc co "widzi" myszka
	DrawStage();

	//powrot do normalnego trybu renderowania
	IntType HitsNumber = glRenderMode(GL_RENDER);
	//IntType nrBledu=glGetError();

	//przywracanie oryginalnej macierzy rzutowania
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);

	//char bufor[256];
	//_itoa(HitsNumber,bufor,10);
	//SetWindowText(HandleWindow,bufor);

	//interpretacja zawartosci bufora zaznaczania
	if (HitsNumber > 0)
	{
		//zwracam obiekt najblizszy kamery
		UnsignedIntType ClosestAtomIndex = MarkingBuffer[3];
		UnsignedIntType ClosestAtomDistance = MarkingBuffer[1];
		IntType CurrentIndex = 0;
		for (IntType HitIndex = 0; HitIndex < HitsNumber; HitIndex++)
		{
			if (MarkingBuffer[CurrentIndex + 1] < ClosestAtomDistance)
			{
				ClosestAtomDistance = MarkingBuffer[CurrentIndex + 1];
				if (MarkingBuffer[CurrentIndex] > 0)
					ClosestAtomIndex = MarkingBuffer[CurrentIndex + 3]; //zdaza sie, ze ilosc nazw=0; w naszym przypadku moze byc tylko jedna nazwa w trafieniu
			}
			CurrentIndex += 3 + MarkingBuffer[CurrentIndex]; //ilosc nazw,z_min,z_max HitIndex nazwy (w ilosci zadanej przez ilosc nazw)
		}
		//SetWindowText(HandleWindow,PDBDataFileObjectPtr->GetAtom(ClosestAtomIndex)->PDBRecord);
		return ClosestAtomIndex;
	}
	else return -1;
}
#pragma endregion

							UnsignedIntType CreateFont1(HWND HandleWindow, const char* FontName, IntType HeightInPixels, bool Bold, bool Italics, IntType FirstCharCode, IntType LastCharCode)
							{
								UnsignedIntType FirstListIndex = glGenLists(LastCharCode + 1 - FirstCharCode);   //Tworzy liste na pelen zestaw czcionek

								HFONT FontHandle = CreateFont(HeightInPixels, 0, 0, 0, Bold ? FW_BOLD : FALSE, Italics ? TRUE : FALSE, FALSE, FALSE, ANSI_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY, FF_DONTCARE | DEFAULT_PITCH, FontName);

								HDC HandleDC = GetDC(HandleWindow);
								HFONT DCFontHandle = (HFONT)SelectObject(HandleDC, FontHandle);   //zwiazanie z naszym kontekstem okna
								//Funkcja WGL (specyficzna dla Windows); 32-spacja; ilosc: 96 - alfabet lacinski, 256 - gdy tez polskie znaki

								wglUseFontBitmaps(HandleDC, FirstCharCode, LastCharCode + 1 - FirstCharCode, FirstListIndex);

								SelectObject(HandleDC, DCFontHandle);   //Wybor czcionki, ktora stworzylismy
								DeleteObject(FontHandle);        //Usuwanie pomocniczego

								return FirstListIndex; //=indeks tablicy list
							}

void PDBWindowGL::DrawActors()
{
	ShowRenderingFrequency();

	const double LengthUnit = 0.1;

	//tworzenie list wyswietlania
	static GLuint DrawingList = NULL;
	if (!glIsList(DrawingList) || RefreshListOfDrawing)
	{
		if (RefreshListOfDrawing)
			glDeleteLists(DrawingList, 1);
		DrawingList = CreateListOfDrawing(LengthUnit);
		RefreshListOfDrawing = false;
	}

								//tworzenie czcionki bitmapowej	
								static UnsignedIntType BitmapFont = NULL;
								//XXX from Varia
								if (BitmapFont == NULL)
									//BitmapFont = CreateFont1(false, HandleWindow, "Calibri", 20, true, false, 32, 255);
									BitmapFont = CreateFont1(HandleWindow, "Calibri", 20, true, false, 0, 255);

	//DrawBonds(PDBDataFileObjectPtr,LengthUnit,true);		
	//DrawAtoms(PDBDataFileObjectPtr,LengthUnit,1,true);

	glCallList(DrawingList);
	int RenderingMode = 0;
	glGetIntegerv(GL_RENDER_MODE, &RenderingMode);
	if (RenderingMode == GL_RENDER)
	{
		if (ChosenAtomIndex >= 0 && ShowPDBSize > 0)
		{
			glColor3f(1, 1, 0);
			DrawChosenAtom(PDBDataFileObjectPointer, LengthUnit, 1.05 * ShowPDBSize, ChosenAtomIndex, false);
			glColor3f(5, 0, 5); //jasniejszy niz bialy!
			DrawChosenAtomDescription(ChosenAtomIndex, BitmapFont);
		}
	}
}

#pragma region Source of lights
//zrodla swiatla
void PDBWindowGL::LightSources()
{
	BackgroundLightIntensity = 0.5f;
	MilkyBulb(0.5f);

	//mgla
	//glEnable(GL_FOG);
	//const float biel[4]={1.0,1.0,1.0,1.0};
	//glFogfv(GL_FOG_COLOR,biel);
	//glFogf(GL_FOG_START,0.0);
	//glFogf(GL_FOG_END,100.0);
	//glFogf(GL_FOG_MODE,GL_LINEAR);
}

void PDBWindowGL::MilkyBulb(float jasnosc)
{
	const float Color[4] = { jasnosc, jasnosc, jasnosc, 1.0f };
	const float Position[4] = { 5.0f, 0.0f, 5.0f, 1.0f };
	glLightfv(GL_LIGHT1, GL_POSITION, Position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, Color);
	glEnable(GL_LIGHT1);
}

void PDBWindowGL::YellowGreenMilkyBulbs()
{
	//zolta mleczna zarowka
	const float YellowColor[4] = { 1.0f, 1.0f, 0.0f, 1.0f };
	const float YellowPosition[4] = { -2.0f, 0.0f, 1.0f, 1.0f };
	glLightfv(GL_LIGHT2, GL_POSITION, YellowPosition);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, YellowColor);
	//glEnable(GL_LIGHT2);

	//zielona mleczna zarowka
	const float GreenColor[4] = { 0.0f, 1.0f, 0.0f, 1.0f };
	const float GreenPosition[4] = { 2.0f, 0.0f, 1.0f, 1.0f };
	glLightfv(GL_LIGHT3, GL_POSITION, GreenPosition);
	glLightfv(GL_LIGHT3, GL_DIFFUSE, GreenColor);
	//glEnable(GL_LIGHT3);
}

void PDBWindowGL::Reflector(float FlareBirghtness, float BrigthnessDisperesed)
{
	const float ColorDispersed[4] = { BrigthnessDisperesed, BrigthnessDisperesed, BrigthnessDisperesed, 1.0f };
	const float ColorFlare[4] = { FlareBirghtness, FlareBirghtness, FlareBirghtness, 1.0 };
	const float Position[4] = { -10.0f, -10.0f, 10.0f, 1.0f };
	const float Direction[4] = { 1.0, 1.0, -1.0, 1.0 };
	const float BundleWidth = 30.0f; //w stopniach
	const float Extinction = 1.0f;

	glLightfv(GL_LIGHT4, GL_POSITION, Position);
	glLightfv(GL_LIGHT4, GL_DIFFUSE, ColorDispersed);

	glLightfv(GL_LIGHT4, GL_SPECULAR, ColorFlare);
	glLightfv(GL_LIGHT4, GL_SPOT_DIRECTION, Direction);
	glLightf(GL_LIGHT4, GL_SPOT_CUTOFF, BundleWidth);
	glLightf(GL_LIGHT4, GL_SPOT_EXPONENT, Extinction);
	//glEnable(GL_LIGHT4);
}
#pragma endregion
