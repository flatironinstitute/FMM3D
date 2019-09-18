#ifndef _SCTL_VECTOR_HPP_
#define _SCTL_VECTOR_HPP_

#include SCTL_INCLUDE(mem_mgr.hpp)
#include SCTL_INCLUDE(common.hpp)

#include <vector>
#include <cstdlib>
#include <cstdint>

// TODO: Implement dynamic-vector which can be created with variable size from a stack memory pool

namespace SCTL_NAMESPACE {

template <class ValueType> class ConstVector {
 public:
  typedef ValueType value_type;
  typedef ValueType& reference;
  typedef const ValueType& const_reference;
  typedef Iterator<ValueType> iterator;
  typedef ConstIterator<ValueType> const_iterator;
  typedef Long difference_type;
  typedef Long size_type;

  ConstVector();

  ConstVector(Long dim_, ConstIterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  ConstVector(const ConstVector& V, bool own_data_ = true);

  ConstVector(const std::vector<ValueType>& V, bool own_data_ = true);

  ~ConstVector();

  void Swap(ConstVector<ValueType>& v1);

  void ReInit(Long dim_, ConstIterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  void Write(const char* fname) const;

  Long Dim() const;

  ConstIterator<ValueType> begin() const;

  ConstIterator<ValueType> end() const;

  // Element access

  const ValueType& operator[](Long j) const;

  // ConstVector-ConstVector operations

  Vector<ValueType> operator+(const ConstVector& V) const;

  Vector<ValueType> operator-(const ConstVector& V) const;

  Vector<ValueType> operator*(const ConstVector& V) const;

  Vector<ValueType> operator/(const ConstVector& V) const;

  // ConstVector-Scalar operations

  Vector<ValueType> operator+(ValueType s) const;

  Vector<ValueType> operator-(ValueType s) const;

  Vector<ValueType> operator*(ValueType s) const;

  Vector<ValueType> operator/(ValueType s) const;

 protected:
  void Init(Long dim_, Iterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  Long dim;
  Long capacity;
  Iterator<ValueType> data_ptr;
  bool own_data;
};

template <class ValueType> class Vector : public ConstVector<ValueType> {
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

  //~Vector();

  void Swap(Vector<ValueType>& v1);

  void ReInit(Long dim_, Iterator<ValueType> data_ = NullIterator<ValueType>(), bool own_data_ = true);

  void Read(const char* fname);

  void SetZero();

  ConstIterator<ValueType> begin() const { return ConstVector<ValueType>::begin(); }

  Iterator<ValueType> begin();

  ConstIterator<ValueType> end() const { return ConstVector<ValueType>::end(); }

  Iterator<ValueType> end();

  void PushBack(const ValueType& x);

  // Element access

  const ValueType& operator[](Long j) const { return this->ConstVector<ValueType>::operator[](j); }

  ValueType& operator[](Long j);

  // Vector-Vector operations

  Vector& operator=(const std::vector<ValueType>& V);

  Vector& operator=(const Vector& V);

  Vector& operator+=(const Vector& V);

  Vector& operator-=(const Vector& V);

  Vector& operator*=(const Vector& V);

  Vector& operator/=(const Vector& V);

  // Vector-Scalar operations

  Vector& operator=(ValueType s);

  Vector& operator+=(ValueType s);

  Vector& operator-=(ValueType s);

  Vector& operator*=(ValueType s);

  Vector& operator/=(ValueType s);

};

template <class ValueType> ConstVector<ValueType> operator+(ValueType s, const ConstVector<ValueType>& V);

template <class ValueType> ConstVector<ValueType> operator-(ValueType s, const ConstVector<ValueType>& V);

template <class ValueType> ConstVector<ValueType> operator*(ValueType s, const ConstVector<ValueType>& V);

template <class ValueType> ConstVector<ValueType> operator/(ValueType s, const ConstVector<ValueType>& V);

template <class ValueType> std::ostream& operator<<(std::ostream& output, const ConstVector<ValueType>& V);

}  // end namespace

#include SCTL_INCLUDE(vector.txx)

#endif  //_SCTL_VECTOR_HPP_
