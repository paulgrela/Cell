
#ifndef VECTOR_TYPE_H
#define VECTOR_TYPE_H

#include <stdexcept>
#include <exception>
#include <math.h>

template<typename T> 
struct VectorType
{
public:
	T X, Y, Z, W;

	VectorType() : X(0), Y(0), Z(0), W(1)
	{
	};

	VectorType(T X, T Y, T Z, T W = 1) : X(X), Y(Y), Z(Z), W(W)
	{
	}

	T& operator[](const int Index)
	{
		switch (Index)
		{
			case 0: return X; break;
			case 1: return Y; break;
			case 2: return Z; break;
			default: throw std::invalid_argument("Index out of range"); break;
		}
	};

	bool operator==(const VectorType& v) const
	{
		return X == v.X && Y == v.Y && Z == v.Z && W == v.W;
	}

	bool operator!=(const VectorType& v) const
	{
		return !(*this == v);
	}

	bool Equivalent(const VectorType& v) const
	{
		if (W == 0 || v.W == 0)
			throw std::runtime_error("One or two points are in infinitt (W=0)");

		return X / W == v.X / v.W && Y / W == v.Y / v.W && Z / W == v.Z / v.W;
	}

	VectorType& operator+=(const VectorType& v)
	{
		if (W == 1)
		{
			X += v.X;
			Y += v.Y;
			Z += v.Z;
		}
		else
		{
			X = X * v.W + v.X * W;
			Y = Y * v.W + v.Y * W;
			Z = Z * v.W + v.Z * W;
			W *= v.W;
		}
		return *this;
	};

	VectorType operator+(const VectorType& v) const
	{
		return VectorType(*this) += v;
	};

	VectorType& operator-=(const VectorType& v)
	{
		if (W == 1)
		{
			X -= v.X;
			Y -= v.Y;
			Z -= v.Z;
		}
		else
		{
			X = X * v.W - v.X * W;
			Y = Y * v.W - v.Y * W;
			Z = Z * v.W - v.Z * W;
			W *= v.W;
		}
		return *this;
	};

	VectorType operator-(const VectorType& v) const
	{
		return VectorType(*this) -= v;
	};

	VectorType& operator*=(const T a)
	{
		X *= a;
		Y *= a;
		Z *= a;

		return *this;
	};

	VectorType operator*(const T a) const
	{
		return VectorType(*this) *= a;
	};

	VectorType& operator/=(const T a)
	{
		if (a == 0.0)
			throw std::invalid_argument("Dzielenie przez zero");

		X /= a;
		Y /= a;
		Z /= a;

		return *this;
	};

	VectorType operator/(const T a) const
	{
		return VectorType(*this) /= a;
	};

	VectorType operator-() const
	{
		return VectorType(-X, -Y, -Z, W);
	}
	VectorType operator+() const
	{
		return *this;
	};

	T operator*(const VectorType& v) const
	{
		return (X * v.X + Y * v.Y + Z * v.Z) / (W * v.W);
	};

	VectorType<T> operator^(const VectorType& v) const
	{
		if (W != 1 || v.W != 1)
			throw std::runtime_error("The vector product is not specified in homogeneous coordinates (the W coordinate must be 1)");

		VectorType Result;
		Result.X = Y * v.Z - Z * v.Y;
		Result.Y = -(X * v.Z - Z * v.X);
		Result.Z = X * v.Y - Y * v.X;
		return Result;
	}

	VectorType<T> operator%(const VectorType& v) const
	{
		return *this ^ v;
	}

	T LengthSquare() const
	{
		if (W == 0)
			throw std::runtime_error("Point is in infinity (W=0)");

		return (X * X + Y * Y + Z * Z) / (W * W);
	}

	T Length() const
	{
		return sqrt(LengthSquare());
	}

	void Normalize()
	{
		T Length = Length();
		if (Length != 0)
			*this /= Length;

		else throw std::invalid_argument("Normalisation of zero vector is impossible");
	}

	VectorType<T> Unormowany() const
	{
		VectorType<T> v(*this);
		v.Normalize();
		return v;
	}

	T* WriteInArray(T Array[3]) const
	{
		if (W == 0)
			throw std::runtime_error("Point is in infinity (W = 0)");

		Array[0] = X / W;
		Array[1] = Y / W;
		Array[2] = Z / W;
		return Array;
	}

	T* ZapiszWTablicy4(T tablica[4]) const
	{
		tablica[0] = X;
		tablica[1] = Y;
		tablica[2] = Z;
		tablica[3] = W;
		return tablica;
	}

	static VectorType<T> Zero()
	{
		return VectorType(0, 0, 0);
	}

	static VectorType<T> WersorX()
	{
		return VectorType(1, 0, 0);
	}

	static VectorType<T> WersorY()
	{
		return VectorType(0, 1, 0);
	}

	static VectorType<T> WersorZ()
	{
		return VectorType(0, 0, 1);
	}

	static VectorType<T> Jeden()
	{
		return VectorType(1, 1, 1);
	}
};

template<typename T> 
VectorType<T> inline operator*(const T a, const VectorType<T>& v)
{
	return v * a;
}

template<typename T> 
VectorType<T> IloczynWektorowy(VectorType<T> a, VectorType<T> b)
{
	return a ^ b;
}

typedef VectorType<float> FloatVectorType;

#endif
