#ifndef _SCTL_COMM_HPP_
#define _SCTL_COMM_HPP_

#include SCTL_INCLUDE(common.hpp)

#include <map>
#include <stack>
#ifdef SCTL_HAVE_MPI
#include <mpi.h>
#endif

namespace SCTL_NAMESPACE {

template <class ValueType> class Vector;

class Comm {

 public:
  enum class CommOp {
    SUM,
    MIN,
    MAX
  };

  Comm();

#ifdef SCTL_HAVE_MPI
  Comm(const MPI_Comm mpi_comm) { Init(mpi_comm); }
#endif

  Comm(const Comm& c);

  static Comm Self();

  static Comm World();

  Comm& operator=(const Comm& c);

  ~Comm();

#ifdef SCTL_HAVE_MPI
  MPI_Comm GetMPI_Comm() { return mpi_comm_; }
#endif

  Comm Split(Integer clr) const;

  Integer Rank() const;

  Integer Size() const;

  void Barrier() const;

  template <class SType> void* Isend(ConstIterator<SType> sbuf, Long scount, Integer dest, Integer tag = 0) const;

  template <class RType> void* Irecv(Iterator<RType> rbuf, Long rcount, Integer source, Integer tag = 0) const;

  void Wait(void* req_ptr) const;

  template <class SType, class RType> void Allgather(ConstIterator<SType> sbuf, Long scount, Iterator<RType> rbuf, Long rcount) const;

  template <class SType, class RType> void Allgatherv(ConstIterator<SType> sbuf, Long scount, Iterator<RType> rbuf, ConstIterator<Long> rcounts, ConstIterator<Long> rdispls) const;

  template <class SType, class RType> void Alltoall(ConstIterator<SType> sbuf, Long scount, Iterator<RType> rbuf, Long rcount) const;

  template <class SType, class RType> void* Ialltoallv_sparse(ConstIterator<SType> sbuf, ConstIterator<Long> scounts, ConstIterator<Long> sdispls, Iterator<RType> rbuf, ConstIterator<Long> rcounts, ConstIterator<Long> rdispls, Integer tag = 0) const;

  template <class Type> void Alltoallv(ConstIterator<Type> sbuf, ConstIterator<Long> scounts, ConstIterator<Long> sdispls, Iterator<Type> rbuf, ConstIterator<Long> rcounts, ConstIterator<Long> rdispls) const;

  template <class Type> void Allreduce(ConstIterator<Type> sbuf, Iterator<Type> rbuf, Long count, CommOp op) const;

  template <class Type> void Scan(ConstIterator<Type> sbuf, Iterator<Type> rbuf, int count, CommOp op) const;

  template <class Type> void PartitionW(Vector<Type>& nodeList, const Vector<Long>* wts_ = nullptr) const;

  template <class Type> void PartitionN(Vector<Type>& v, Long N) const;

  template <class Type> void PartitionS(Vector<Type>& nodeList, const Type& splitter) const;

  template <class Type> void HyperQuickSort(const Vector<Type>& arr_, Vector<Type>& SortedElem) const;

  template <class Type> void SortScatterIndex(const Vector<Type>& key, Vector<Long>& scatter_index, const Type* split_key_ = nullptr) const;

  template <class Type> void ScatterForward(Vector<Type>& data_, const Vector<Long>& scatter_index) const;

  template <class Type> void ScatterReverse(Vector<Type>& data_, const Vector<Long>& scatter_index_, Long loc_size_ = 0) const;

 private:
  template <typename A, typename B> struct SortPair {
    int operator<(const SortPair<A, B>& p1) const { return key < p1.key; }
    A key;
    B data;
  };

#ifdef SCTL_HAVE_MPI
  void Init(const MPI_Comm mpi_comm);

  Vector<MPI_Request>* NewReq() const;

  void DelReq(Vector<MPI_Request>* req_ptr) const;

  mutable std::stack<void*> req;

  int mpi_rank_;
  int mpi_size_;
  MPI_Comm mpi_comm_;

  /**
   * \class CommDatatype
   * \brief An abstract class used for communicating messages using user-defined
   * datatypes. The user must implement the static member function "value()" that
   * returns the MPI_Datatype corresponding to this user-defined datatype.
   * \author Hari Sundar, hsundar@gmail.com
   */
  template <class Type> class CommDatatype {
   public:
    static MPI_Datatype value() {
      static bool first = true;
      static MPI_Datatype datatype;
      if (first) {
        first = false;
        MPI_Type_contiguous(sizeof(Type), MPI_BYTE, &datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }

    static MPI_Op sum() {
      static bool first = true;
      static MPI_Op myop;

      if (first) {
        first = false;
        int commune = 1;
        MPI_Op_create(sum_fn, commune, &myop);
      }

      return myop;
    }

    static MPI_Op min() {
      static bool first = true;
      static MPI_Op myop;

      if (first) {
        first = false;
        int commune = 1;
        MPI_Op_create(min_fn, commune, &myop);
      }

      return myop;
    }

    static MPI_Op max() {
      static bool first = true;
      static MPI_Op myop;

      if (first) {
        first = false;
        int commune = 1;
        MPI_Op_create(max_fn, commune, &myop);
      }

      return myop;
    }

   private:
    static void sum_fn(void* a_, void* b_, int* len_, MPI_Datatype* datatype) {
      Type* a = (Type*)a_;
      Type* b = (Type*)b_;
      int len = *len_;
      for (int i = 0; i < len; i++) {
        b[i] = a[i] + b[i];
      }
    }

    static void min_fn(void* a_, void* b_, int* len_, MPI_Datatype* datatype) {
      Type* a = (Type*)a_;
      Type* b = (Type*)b_;
      int len = *len_;
      for (int i = 0; i < len; i++) {
        if (a[i] < b[i]) b[i] = a[i];
      }
    }

    static void max_fn(void* a_, void* b_, int* len_, MPI_Datatype* datatype) {
      Type* a = (Type*)a_;
      Type* b = (Type*)b_;
      int len = *len_;
      for (int i = 0; i < len; i++) {
        if (a[i] > b[i]) b[i] = a[i];
      }
    }
  };
#else
  mutable std::multimap<Integer, ConstIterator<char>> send_req;
  mutable std::multimap<Integer, Iterator<char>> recv_req;
#endif
};

}  // end namespace

#include SCTL_INCLUDE(comm.txx)

#endif  //_SCTL_COMM_HPP_
