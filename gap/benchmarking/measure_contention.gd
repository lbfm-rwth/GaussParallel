#! @Chapter Finding a suitable number of blocks
#! @Section Measure Contention
#! @Arguments numberBlocks, q, A, [showOutput]
#!
#! @Description
#! This function helps with finding a suitable value for <A>numberBlocks</A>,
#! which has a big influence on the performance of the
#! <Ref Func="EchelonMatTransformationBlockwise"/> function.
#! For inputs <A>numberBlocks</A>, a field size <A>q</A>, and a matrix <A>A</A>
#! over <C>GF(q)</C>, this function does the following:
#! - it calls the parallel function
#!   <Ref Func="EchelonMatTransformationBlockwise"/>, 
#! - it calls the sequential function
#!   <C>EchelonMatTransformation</C> from the package <C>Gauss</C>,
#! - it prints the wall time, that is the elapsed real time, each call took,
#! - for the parallel function, it prints an estimate of the lock contenation
#!   ratio, a short explanation of lock contention is given below. A high lock
#!   contention ratio deteriorates the performance of the algorithm. 
#!
#! The input <A>showOutput</A> can be used to suppress printing of messages.
#!
#! The influence of <A>numberBlocks</A> on the performance is as follows:
#! - if <A>numberBlocks</A> is too small, then not enough calculations can
#!   happen in parallel
#! - if <A>numberBlocks</A> is too big then:
#!   - the lock contention ratio increases
#!   - HPC-GAP's task framework generates a big overhead since it is not
#!     optimized.
#!
#! In a parallel computation several threads may try to access an
#! object at the same time. HPC-GAP needs to prevent situations in which one
#! thread writes into an object that another thread is currently reading, since
#! that could lead to corrupted data. To this end, HPC-GAP can synchronize
#! access to objects via locks.
#!
#! If a thread reads an object, that thread can acquire a lock for that object,
#! that is no other thread can write into it until the lock is released. If a
#! thread tries to write into said object before the lock is released we say
#! that the lock was contended. In such a situation the contending thread needs
#! to wait until the lock is released. This leads to waiting times or
#! unnecessary context-switches and deteriorates performance if it happens too
#! often.
#!
#! @BeginExampleSession
#! gap> n := 4000;; numberBlocks := 8;; q := 5;;
#! gap> A := RandomMat(n, n, GF(q));;
#! gap> GAUSS_MeasureContention(numberBlocks, q, A);
#! Make sure you called GAP with sufficient preallocated memory via `-m` if you
#! try non-trivial examples! Otherwise garbage collection will be a big
#! overhead.
#! 
#! Starting the parallel algorithm
#! EchelonMatTransformationBlockwise.
#!
#! Wall time  parallel execution (ms): 33940.
#! CPU  time  parallel execution (ms): 2
#! Lock statistics(estimates):
#! acquired - 181502, contended - 1120, ratio - 1.%
#! Locks acquired and contention counters per thread
#! [ thread, locks acquired, locks contended ]:
#! [ 5, 54093, 248 ]
#! [ 6, 51004, 228 ]
#! [ 7, 37102, 389 ]
#! [ 8, 39303, 255 ]
#! 
#! Starting the sequential algorithm
#! EchelonMatTransformation
#!
#! Wall time Gauss pkg execution (ms): 66778.
#! 
#! Speedup factor (sequential / parallel wall time):
#! 1.96
#! @EndExampleSession
DeclareGlobalFunction("GAUSS_MeasureContention");
