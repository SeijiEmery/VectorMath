//
//  Vec3.h
//
//  Created by Seiji Emery on 9/12/14.
//  Copyright (c) 2014 Seiji Emery. All rights reserved.
//
//  Defines a generic Vec3 class with full operator overloading (with both mutable and
//  immutable versions of vector operations), and utility methods like distance() and
//  project().
//
//  Implementation is header-only, so no additional source files are required.
//

#ifndef Math_Vec3_h
#define Math_Vec3_h

// Includes - I/O, sqrt, acos (can't just predeclare, since this is a full header impl).
#include <ostream>
#define _USE_MATH_DEFINES
#include <cmath>


// Yes, macros are evil, but when used to reduce boilerplate and code repetition they can
// be a good thing imo – DRY is important and by defining the semantics for most vector
// operations in one place (or just once) we can reduce the potential for human error, etc.
// (More important for bigger projects / large classes w/ tons of overloads and lots of
//  boilerplate).

// Performs a binary operation on the components of two vectors in sequence
#define DO_VECTOR_OP(lhs, op, rhs) \
    (((lhs).x op (rhs).x), \
     ((lhs).y op (rhs).y), \
     ((lhs).z op (rhs).z))

// Performs a binary operation between the components of a vector and a scalar
#define DO_SCALAR_OP(lhs, op, rhs) \
    (((lhs).x op (rhs)), \
     ((lhs).y op (rhs)), \
     ((lhs).z op (rhs)))

// Note: normally, I'd only ever use macros in the implementation file and keep the header
// clean and readable, but in this case we're using a header implementation (for simplicity),
// so...


/// Defines a generic 3-dimensional vector class with operator overloads and methods for
/// magnitude, dot/cross product, etc.
template <typename T>
class Vec3 {
public:
    // Constructors (virtual destructor not needed since we won't be subclassing Vec3)
    Vec3 () {}
    Vec3 (T x, T y, T z) : x(x), y(y), z(z) {}
    Vec3 (const Vec3 & other) : x(other.x), y(other.y), z(other.z) {}
    
    // Assignment operator
    Vec3 & operator = (const Vec3 & other) {
        return DO_VECTOR_OP(*this, =, other), *this;
    }
    
    // Friend functions used to impement immutable binary operations.
    
    /// Provides ostream << writing/printing semantics for this vector.
    friend std::ostream & operator << (std::ostream & os, const Vec3 &v) {
        return os << "Vec3 { " << v.x << ", " << v.y << ", " << v.z << " }", os;
    }
    
    /// Performs immutable vector addition between two vectors.
    friend Vec3 operator + (const Vec3 & a, const Vec3 & b) {
        return { a.x + b.x, a.y + b.y, a.z + b.z };
    }
    
    /// Performs immutable vector subtraction between two vectors.
    friend Vec3 operator - (const Vec3 & a, const Vec3 & b) {
        return { a.x - b.x, a.y - b.y, a.z - b.z };
    }
    
    /// Performs immutable scalar multiplication between a vector and a constant.
    friend Vec3 operator * (const Vec3 & a, const T c) {
        return { a.x * c, a.y * c, a.z * c };
    }
    
    // Overload to provide commutative multiplication semantics (otherwise myVec * 3.0
    // is defined, but 3.0 * myVec is not). Not an issue with the other vector operations,
    // where the lhs and rhs types are the same.
    friend Vec3 operator * (const T c, const Vec3 & v) { return v * c; }
    

    /// Performs immutable scalar division between a vector and a constant.
    ///
    /// Note: will *not* prevent divide by zero errors. Not an issue for floating point
    /// (since divide by zero is well-defined by IEEE), but this is a near irrecoverable
    /// error for integers – so be aware of that if using Vec3<int>.
    friend Vec3 operator / (const Vec3 & a, const T c) {
        return a * (1.0f / c);
    }
    
    // Use methods for mutable (in-place) operations
    
    /// Perform in-place vector addition on this vector.
    Vec3 & operator += (const Vec3 & other) {
        return DO_VECTOR_OP(*this, +=, other), *this;
    }
    
    /// Perform in-place vector subtraction on this vector.
    Vec3 & operator -= (const Vec3 & other) {
        return DO_VECTOR_OP(*this, -=, other), *this;
    }
    
    /// Perform in-place scalar multiplication on this vector.
    Vec3 & operator *= (const T s) {
        return DO_SCALAR_OP(*this, *=, s), *this;
    }
    
    /// Perform in-place scalar division on this vector.
    ///
    /// Note: will *not* prevent divide by zero errors. Not an issue for floating point
    /// (since divide by zero is well-defined by IEEE), but this is a near irrecoverable
    /// error for integers – so be aware of that if using Vec3<int>.
    Vec3 & operator /= (const T s) {
        return *this *= (1.0f / s);
    }
    
    /// Returns the cross product of this and another vector (does not modify self).
    Vec3 cross (const Vec3 & other) const {
        return Vec3 {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }
    
    /// Returns a normalized vector (does not modify self).
    Vec3 normalize () const {
        // return *this / m;
        
        auto m = magnitude();
        
        // Have to handle two cases:
        //   m = 0 => return zero vector (self)
        //   m = 1 => can just return self (since self / 1 = self)

        if (m == 0 || m == 1)
            return { *this };
        return *this / m;
    }
    
    /// Normalizes this vector and returns a reference for operator chaining.
    Vec3 & normalized () {
        // return *this /= m;
    
        auto m = magnitude();
        
        // Same case as above
        if (m == 0 || m == 1)
            return *this;
        return *this /= m;
    }
    
    /// Returns the dot product of this and another vector (does not modify self).
    T dot (const Vec3 & other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    /// Returns this vector projected onto another vector. Does not modify self.
    Vec3 project (const Vec3 & other) const {
        return (dot(other) / magitudeSquared(other)) * other;
    }
    
    /// Returns the rejection of this vector onto another vector (does not modify self).
    Vec3 reject (const Vec3 & other) const {
        return *this - project(other);
    }
    
    // Returns the cosTheta value between this and another vector.
    T cosTheta (const Vec3 & other) const {
        // Note: had to increase complexity to handle divide by zero cases
        // return dot(other) / (magnitude() * other.magnitude());
        
        auto m = magnitude() * other.magnitude();
        auto dp = dot(other);
        
        // If product of magnitudes = 1, just return dot product (since would be dividing by 1)
        // If either vector is zero (magnitude(s) = 0), return 0.
        //      Since zero vector <=> m = 0 <=> dot product = 0, can just return the
        //      dot product using same case as m = 1.
        if (m == 0 || m == 1)
            return dp;
        return dp / m;
    }
    
    // Returns the angle (in radians) between this and another vector.
    auto angle (const Vec3 & other) const -> decltype(acos(cosTheta(other))) {
        return acos(cosTheta(other));
    }
    
    /// Returns the magnitude squared of this vector.
    /// Slightly cheaper than magnitude() (no sqrt()).
    T magnitudeSquared () const {
        return x * x + y * y + z * z;
    }
    
    /// Returns the magnitude of this vector.
    T magnitude () const {
        return sqrt(magnitudeSquared());
    }
    
    /// Alias for vector magnitude().
    T length () const { return magnitude(); }
    
    /// Returns the distance squared between two vector positions.
    /// Slightly cheaper than distance() (no sqrt).
    T distanceSquared (const Vec3 & other) const {
        auto dx = x - other.x;
        auto dy = y - other.y;
        auto dz = z - other.z;
        return dx * dx + dy * dy + dz * dz;
    }
    
    /// Returns the distance between two vector positions.
    T distance (const Vec3 & other) const {
        return sqrt(distanceSquared(other));
    }
    
protected:
    T x = 0;
    T y = 0;
    T z = 0;
};

// Don't forget to remove our helper macros (important since we're in a header file)
#undef DO_VECTOR_OP
#undef DO_SCALAR_OP

// Common typedefs
typedef Vec3<float>  Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int>    Vec3i; // Probably shouldn't really use (no real point to use integer vectors, and undefined stuff (eg divide by zero, etc) is headache inducing), but whatever...

void foo () {
    auto bar = Vec3f { 1, 2, 3 } .normalize() * 3.0 + Vec3f { 1, 2, 3 };
    auto baz = bar * 3.0;
//    auto borg = bar * baz;
    
}



#endif
