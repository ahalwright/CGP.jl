# Test of why threading runs slower
using Base.Threads
using BenchmarkTools

function time_threading( nreps::Int64, p::Parameters, numcircuits::Int64 )
  @btime test_threading( nreps, p, numcircuits )
end

function test_threading( nreps::Int64, p::Parameters, numcircuits::Int64 )
  println("nthreads: ",nthreads())
  println("nprocs: ",nprocs())
  if nprocs() > 2
    proc_nreps = div(nreps,nprocs()-1)
    println("proc_nreps: ",proc_nreps)  
    pmap(x->test_worker( nreps, p, numcircuits ), collect(1:nprocs()-1))
  elseif Threads.nthreads() > 1
    thread_nreps = div(nreps,Threads.nthreads())
    println("thread_nreps: ",thread_nreps)
    Threads.@threads for i = 1:Threads.nthreads()
      test_worker( thread_nreps, p, numcircuits )
    end
  else
    println("nreps: ",nreps)
    test_worker( nreps, p, numcircuits )
  end
end

function test_worker( nreps::Int64, p::Parameters, numcircuits::Int64 )
  count_outputs( nreps, p, numcircuits, use_lincircuit=false )
  nothing
end
  
