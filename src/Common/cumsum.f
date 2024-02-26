c  This file contains the following user callable
c  routines:
c
c  * cumsum:
c      computes the cumulative sum (aka prefix
c      sum or scan) of an array of integers. This code
c      does some memory management and basic decision
c      making around parallelism, which is hidden from
c      the caller
c      
c     Advanced user interfaces:
c******************************
c     Note for developers: the prefix sum for integer
c     arrays on a shared memory machine can be memory
c     bound and performance changes significantly depending
c     on the speed and size of the L1,L2, and L3 cache
c     on chip. Thus, the performance of the parallel code
c     typically exhibits non-monotonic behavior. The efficiency
c     increases with the length of the vectors as the
c     overhead is ammortized but then the memory bound
c     effects kicks in for larger vectors. The best-case
c     improvement is nthreads/2, though it is observed
c     only for the appropriate length vectors.
c******************************      
c     
c     cumsum1 - same as cumsum but forces the serial
c       code
c
c     cumsum_para - openmp implementation of parallel
c       cumulative sum
c
c


      
      subroutine cumsum(n,a,b)
c
c--------------------------
c  This subroutine computes the cumulative sum
c  (aka prefix sum aka scan) of the vector a and
c  returns in the result in b. i.e.
c
c  b(i) = sum_{j\leq i} a(j) for i = 1,...,n
c
c  This is the memory management routine with
c  the work done by either cumsum1 in serial or
c  cumsum_para in parallel
c
c  Input arguments:
c   
c    - n: integer
c        length of a and b
c    - a: integer(n)
c        the array to perform cumulative sum on
c
c  Output arguments:
c 
c    - b: integer(n)
c        the resulting cumulative sum array
c--------------------------
c
c
      implicit none
      integer, intent(in) :: n,a(n)
      integer, intent(out) :: b(n)
      
      integer lend, nsmall, maxth, offset
      parameter (lend = 200)
      integer d(lend)
      integer, allocatable :: d2(:)
c$    integer omp_get_max_threads
      data nsmall / 10000 /

c     if small problem, don't do parallel
      if (n .lt. nsmall) goto 300

c     get upper bound of number of threads on hand

      maxth = 1      
c$    maxth = omp_get_max_threads()

c     no benefit to parallelize if only 2 threads
      if (maxth .le. 2) goto 300

c     if not a ton of processors, use d on stack
      if (maxth .le. lend) goto 200

c     if tons of processors, allocate d2
      allocate(d2(maxth))
      call cumsum_para(n,a,b,maxth,d2)
      return

 200  continue

      call cumsum_para(n,a,b,lend,d)
      return
      
 300  continue
            
      call cumsum1(n,a,b)
      return
      
      end
c
c
      subroutine cumsum1(n,a,b)
c
c--------------------------
c  This subroutine computes the cumulative sum
c  (aka prefix sum aka scan) of the vector a and
c  returns in the result in b. i.e.
c
c  b(i) = sum_{j\leq i} a(j) for i = 1,...,n
c  
c  This is the serial version of the code.
c 
c  Input arguments:
c   
c    - n: integer
c        length of a and b
c    - a: integer(n)
c        the array to perform cumulative sum on
c
c  Output arguments:
c 
c    - b: integer(n)
c        the resulting cumulative sum array
c--------------------------
c      
      
      implicit none
      integer, intent(in) :: n,a(n)
      integer, intent(out) :: b(n)

      integer i,isum

      isum = 0

      do i=1,n
        isum = isum + a(i)
        b(i) = isum
      enddo
      
      return
      end
c
      
      subroutine cumsum_para(n,a,b,nd,d)
c
c
c--------------------------
c  This subroutine computes the cumulative sum
c  (aka prefix sum aka scan) of the vector a and
c  returns in the result in b. i.e.
c
c  b(i) = sum_{j\leq i} a(j) for i = 1,...,n
c  
c  This is the openmp parallel version of the code
c
c  Input arguments:
c   
c    - n: integer
c        length of a and b
c    - a: integer(n)
c        the array to perform cumulative sum on
c    - nd: integer
c        length of the work array d, nd should be larger
c        than the number of threads to be used
c
c  Output arguments:
c
c    - b: integer(n)
c        the resulting cumulative sum array
c    - d: integer(nd)
c        work array of length nd, nd should be larger than
c        the number of threads to be used
c--------------------------
c     
c  WARNING: this subroutine depends on the properties
c  of the "static" schedule in the OPENMP
c  specification. Specifically, the static schedule
c  with no chunk_size specified assigns at most one
c  block of iterates to each thread in thread order.
c  The assignment of indices is the same for any do
c  loops of the same cardinality within a parallel
c  region. Changing the schedule will make the code non
c  conforming.
c
      implicit none
      integer, intent(in) :: n,a(n),nd
      integer, intent(out) :: b(n),d(nd)

      integer i,isum,offset
      integer nth, id, inext
c$    integer omp_get_thread_num, omp_get_num_threads      


c$OMP PARALLEL DEFAULT(none) SHARED(a,b,n,d)
c$OMP$ PRIVATE(i,isum,offset,id,inext,nth)

      id = 1
c$    id = omp_get_thread_num()+1

c     compute cumulative sums of portions (parallel)
      
      
      isum = 0
c$OMP DO SCHEDULE(static)       
      do i = 1,n
         isum = isum+a(i)
         b(i) = isum
      enddo
c$OMP END DO no wait      
      d(id) = isum

c$OMP BARRIER

c     accumulate the ends of the partial sums (each thread)
      
      offset = 0
      do i = 1,id-1
         offset = offset + d(i)
      enddo

c$OMP DO SCHEDULE(static)       
      do i = 1,n
         b(i) = b(i) + offset
      enddo
c$OMP END DO no wait      
      
c$OMP END PARALLEL
      
      return
      end
c
      
