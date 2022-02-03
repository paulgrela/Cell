
#include <windows.h>											
#include <gl\gl.h>												

#include "math.h"                                            

#include "ArcBall.h"

void ArcBall_t::_mapToSphere(const Point2fT* NewPt, Vector3fT* NewVec) const
{
	Point2fT TempPt;
	GLfloat length;

	TempPt = *NewPt;

	//Adjust point coords and scale down to range of [-1 ... 1]
	TempPt.s.X = (TempPt.s.X * this->AdjustWidth) - 1.0f;
	TempPt.s.Y = 1.0f - (TempPt.s.Y * this->AdjustHeight);

	//Compute the square of the length of the vector to the point from the center
	length = (TempPt.s.X * TempPt.s.X) + (TempPt.s.Y * TempPt.s.Y);

	//If the point is mapped outside of the sphere... (length > radius squared)
	if (length > 1.0f)
	{
		GLfloat norm;

		//Compute a normalizing factor (radius / sqrt(length))
		norm = 1.0f / FuncSqrt(length);

		//Return the "normalized" vector, a point on the sphere
		NewVec->s.X = TempPt.s.X * norm;
		NewVec->s.Y = TempPt.s.Y * norm;
		NewVec->s.Z = 0.0f;
	}
	else
	{
		//Return a vector to a point mapped inside the sphere sqrt(radius squared - length)
		NewVec->s.X = TempPt.s.X;
		NewVec->s.Y = TempPt.s.Y;
		NewVec->s.Z = FuncSqrt(1.0f - length);
	}
}

ArcBall_t::ArcBall_t(GLfloat NewWidth, GLfloat NewHeight)
{
	this->StVec.s.X = 0.0f;
	this->StVec.s.Y = 0.0f;
	this->StVec.s.Z = 0.0f;

	this->EnVec.s.X = 0.0f;
	this->EnVec.s.Y = 0.0f;
	this->EnVec.s.Z = 0.0f;

	this->setBounds(NewWidth, NewHeight);
}

void ArcBall_t::click(const Point2fT* NewPt)
{
	this->_mapToSphere(NewPt, &this->StVec);
}

//Mouse drag, calculate rotation
void ArcBall_t::drag(const Point2fT* NewPt, Quat4fT* NewRot)
{
	//Map the point to the sphere
	this->_mapToSphere(NewPt, &this->EnVec);

	//Return the quaternion equivalent to the rotation
	if (NewRot)
	{
		Vector3fT  Perp;

		//Compute the vector perpendicular to the begin and end vectors
		Vector3fCross(&Perp, &this->StVec, &this->EnVec);

		//Compute the length of the perpendicular vector
		if (Vector3fLength(&Perp) > Epsilon)
		{
			//We're ok, so return the perpendicular vector as the transform after all
			NewRot->s.X = Perp.s.X;
			NewRot->s.Y = Perp.s.Y;
			NewRot->s.Z = Perp.s.Z;
			//In the quaternion values, w is cosine (theta / 2), where theta is rotation angle
			NewRot->s.W = Vector3fDot(&this->StVec, &this->EnVec);
		}
		else
		{
			//The begin and end vectors coincide, so return an identity transform
			NewRot->s.X = NewRot->s.Y = NewRot->s.Z = NewRot->s.W = 0.0f;
		}
	}
}

