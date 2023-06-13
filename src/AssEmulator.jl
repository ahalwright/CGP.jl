# Assembly language emulator for genotype-phenotype hypothesis for assembly language sorting programs

#=
# Load only once
mutable struct Data
  N::Int64
  registers::Vector{Int64}
  ZF::Bool   # Zero flag result of cmp  # see https://www.aldeid.com/wiki/X86-assembly/Instructions/cmp
  CF::Bool   # Carry flag result of cmp
end  

# Assume that Data dd is defined externally
struct Instruction
  operation::Function
  operand1::Int64
  operand2::Int64
end
=#  

# Defins a random Data struct
# See also function random_data( N::Int64 )
function emulator_struct( N::Int64, init_values::Vector{Int64} )
  push!(init_values,rand(1:N))
  Data( N, init_values, rand([false,true]), rand([false,true]))
end

function set( dd::Data, value::Int64, to::Int64 )
  dd.registers[to] = value
  dd
end

# Do nothing, a null instruction
function pass( dd::Data, dest::Int64, src::Int64 )
  dd
end

function mov( dd::Data, src::Int64, dest::Int64 )
  dd.registers[dest] = dd.registers[src]
  dd
end

function cmp( dd::Data, dest::Int64, src::Int64 )
  dd.ZF = (dd.registers[dest] == dd.registers[src])
  dd.CF = (dd.registers[dest] > dd.registers[src])
  dd
end

function cmovg( dd::Data, src::Int64, dest::Int64 )
  if dd.CF 
    dd.registers[dest] = dd.registers[src]
  end
  dd
end

function cmovl( dd::Data, src::Int64, dest::Int64 )
  if !dd.CF 
    dd.registers[dest] = dd.registers[src]
  end
  dd
end

function ADprog( dd::Data )  # AlphaDev from Figure 3c of paper.  Fails for ddg.
  (P,Q,R,S)=(1,2,3,4)
  A=dd.registers[1]; B=dd.registers[2]; C=dd.registers[3];
  mov( dd, R, S )
  cmp( dd, P, R )
  cmovg( dd, P, R )  # R == max(A,C)
  @assert dd.registers[R] == max(A,C)
  cmovl( dd, P, S )  # S == min(A,C)
  @assert dd.registers[S] == min(A,C)
  cmp( dd, S, Q )
  cmovg( dd, Q, P )  # P == min(A,B)
  #@assert dd.registers[P] == min(A,B)  # Fails when C < A (ddg)
  cmovg( dd, S, Q )  # Q == max(min(A,C),B)
  @assert dd.registers[Q] == max(min(A,C),B)
  dd
end

function ORprog( dd::Data )  # Original from Figure 3b of paper.  Succeeds for ddl and ddg.
  (P,Q,R,S)=(1,2,3,4)
  A=dd.registers[1]; B=dd.registers[2]; C=dd.registers[3];
  mov( dd, R, S )
  cmp( dd, P, R )
  cmovg( dd, P, R )  # R == max(A,C)
  @assert dd.registers[R] == max(A,C)
  cmovl( dd, P, S )  # S == min(A,C)
  @assert dd.registers[S] == min(A,C)
  mov( dd, S, P )
  cmp( dd, S, Q )
  cmovg( dd, Q, P )
  @assert dd.registers[P] == min(A,B,C)
  cmovg( dd, S, Q )
  @assert dd.registers[Q] == max(min(A,C),B)
  dd
end

function ADfunct( dd::Data )
  (P,Q,R,S)=(1,2,3,4)
  A=dd.registers[1]; B=dd.registers[2]; C=dd.registers[3];
  PP = A
  QQ = B
  RR = max(A,C)
  SS = min(A,C)
  if SS > QQ
    PP = QQ
  end
  PP == min(A,B)  # Fails when C < A (ddg)
end

# Assumes dd is defined externally
function random_instruction(N::Int64)
  Instruction( rand([mov,cmp,cmovl,cmovg]), rand(1:N+1), rand(1:N+1) )
end

# Has dd as an argument
function random_instruction(N::Int64, dd::Data)
  Instruction( dd, rand([mov,cmp,cmovl,cmovg]), rand(1:N+1), rand(1:N+1) )
end

# Returns a random Data struct
function random_data( N::Int64 )
  emulator_struct( N, nthperm( collect(1:N),rand(1:factorial(N))) )
end

# Assumes dd is defined externally
function random_prog( N::Int64, len::Int64 )::Vector{Instruction}
  prog = Instruction[]
  for i = 1:len
    push!(prog, random_instruction( N ) )
  end
  prog
end

# Assumes dd is defined externally
function execute_instruction( inst::Instruction )
  inst.operation( dd, inst.operand1, inst.operand2 )
end

# Assumes dd is defined externally
function execute_assprog( prog::Vector{Instruction} )::Data
  for inst in prog
    dd = execute_instruction( inst )
    #println("inst: ",inst,"   dd: ",dd)
  end
  dd
end

# Returns the count of how many instructions can be deleted/replaced without changing the N-tuple output of the program.
# If delete_inst==true, an instruction is deleted by replacing it with the "pass" instruction which does nothing.
# If delete_inst==false, an instruction is replace by replacing it with a random instruction
function robustness( prog::Vector{Instruction}, N::Int64; delete_inst::Bool=true )
  global dd
  dd0 = deepcopy(dd)
  base_result_dd = execute_assprog( prog )
  bvect = base_result_dd.registers[1:N]
  println("i: ",0,"  bvect: ",bvect)
  robust_count = 0
  for i = 1:length(prog)
    dd = deepcopy(dd0)
    iprog = deepcopy( prog )
    iprog[i] = delete_inst ? Instruction(pass, 1, 1 ) : random_instruction( N, dd )
    iresult = execute_assprog( iprog )
    ivect = iresult.registers[1:N]
    println("i: ",i,"  ivect: ",ivect)
    if ivect == bvect
      robust_count += 1
    end
  end
  dd = deepcopy(dd0)
  (length(prog), robust_count)
end

# All length K Vectors with values from 1:N
# all_vects( N, N ) returns all N^N length k vectors with values from 1:N
# Example:  all_vects(2,2) returns [[1,1],[1,2],[2,1],[2,2]]
function all_vects( N::Int64, k::Int64 )
  if k == 1
    return map(x->[x],collect(1:N))
  end
  result = Vector{Int64}[]
  for u in all_vects( N, k-1 ) 
    for v in  map(x->[x],collect(1:N) )
      push!(result, vcat(u,v) )
    end
  end
  result
end

# Assumes dd is defined externally
# returns an N^N vector of the number of occurences each of the N^N possilbe vectors obtained by running nreps random assembler programs of length len
function count_vects( N::Int64, reps::Int64, len::Int64 )
  global dd
  all_vecs = all_vects( N, N );
  rvect = dd.registers[1:N]; ff = findfirst(x->x==rvect,all_vecs)
  println("count_vects(): dd: ",dd,"  ff: ",ff)
  counts = zeros( Int64, length( all_vecs ) );
  dd0 = deepcopy(dd);
  for r = 1:reps
    dd = deepcopy(dd0)  # restore previous dd
    rprog = random_prog( N, len )
    #print("r: ",r,"  dd: ",dd)
    ddresult = execute_assprog( rprog )
    #println("  rprog: ",rprog,"  ddresult: ",ddresult)
    rvect = ddresult.registers[1:N]
    ff = findfirst(x->x==rvect,all_vecs)
    if ff == nothing
      println( "ddresult: ",ddresult )
    end
    counts[ff] += 1
  end
  dd = deepcopy(dd0)  # restore previous dd
  counts
end 
