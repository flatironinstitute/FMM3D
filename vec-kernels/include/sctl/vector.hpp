#ifndef _SCTL_VECTOR_HPP_
#define _SCTL_VECTOR_HPP_

#include SCTL_INCLUDE(mem_mgr.hpp)
#include SCTL_INCLUDE(common.hpp)

#include <vector>
#include <cstdlib>
#include <cstdint>

namespace SCTL_NAMESPACE {

template <class ValueType> class Vector {
 public:
  typedef ValueType value_type;
  typedef ValueType& reference;
  typedef const ValueType& const_reference;
  typedef Iterator<ValueType> iterator;
  typedef ConstIterator<ValueType> const_iterator;
  typedef Long difference_type;
  typedef Long size_type;

  Vector();

  Vector(Long dim_, Iterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  Vector(const Vector& V);

  Vector(const std::vector<ValueType>& V);

  ~Vector();

  void Swap(Vector<ValueType>& v1);

  void ReInit(Long dim_, Iterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  void Write(const char* fname) const;

  void Read(const char* fname);

  Long Dim() const;

  Long Capacity() const;

  void SetZero();

  Iterator<ValueType> begin();

  ConstIterator<ValueType> begin() const;

  Iterator<ValueType> end();

  ConstIterator<ValueType> end() const;

  void PushBack(const ValueType& x);

  // Element access

  ValueType& operator[](Long j);

  const ValueType& operator[](Long j) const;

  // Vector-Vector operations

  Vector& operator=(const std::vector<ValueType>& V);

  Vector& operator=(const Vector& V);

  Vector& operator+=(const Vector& V);

  Vector& operator-=(const Vector& V);

  Vector& operator*=(const Vector& V);

  Vector& operator/=(const Vector& V);

  Vector operator+(const Vector& V) const;

  Vector operator-(const Vector& V) const;

  Vector operator*(const Vector& V) const;

  Vector operator/(const Vector& V) const;

  // Vector-Scalar operations

  Vector& operator=(ValueType s);

  Vector& operator+=(ValueType s);

  Vector& operator-=(ValueType s);

  Vector& operator*=(ValueType s);

  Vector& operator/=(ValueType s);

  Vector operator+(ValueType s) const;

  Vector operator-(ValueType s) const;

  Vector operator*(ValueType s) const;

  Vector operator/(ValueType s) const;

 private:
  void Init(Long dim_, Iterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  Long dim;
  Long capacity;
  Iterator<ValueType> data_ptr;
  bool own_data;
};

template <class ValueType> Vector<ValueType> operator+(ValueType s, const Vector<ValueType>& V);

template <class ValueType> Vector<ValueType> operator-(ValueType s, const Vector<ValueType>& V);

template <class ValueType> Vector<ValueType> operator*(ValueType s, const Vector<ValueType>& V);

template <class ValueType> Vector<ValueType> operator/(ValueType s, const Vector<ValueType>& V);

template <class ValueType> std::ostream& operator<<(std::ostream& output, const Vector<ValueType>& V);

}  // end namespace

#include SCTL_INCLUDE(vector.txx)

#endif  //_SCTL_VECTOR_HPP_
