# Sets up to test random_neutral_walk() which is in src/Evolve.jl.

#include("../test/miller_examples.jl")

# tries to find a small numactive 3-bit multiplier circuit starting with
#   saved 3-bit multiplier circuits
g3m = goal_for_three_bit_multiplier()
funcs = [AND, XOR ]
c3m1 = build_chromosome( # First saved 3-bit multiplier circuit
       (1,2,3,4,5,6),
       ((XOR,Integer[3, 6]),(AND,Integer[6, 3]),(AND,Integer[5, 2]),(AND,Integer[3, 9]),(AND,Integer[4, 7]),
        (AND,Integer[1, 5]),(AND,Integer[4, 2]),(AND,Integer[1, 4]),(AND,Integer[5, 3]),(AND,Integer[8, 14]),
        (AND,Integer[6, 2]),(XOR,Integer[12, 16]),(XOR,Integer[4, 1]),(XOR,Integer[19, 10]),(XOR,Integer[11, 9]),
        (AND,Integer[6, 20]),(AND,Integer[11, 16]),(AND,Integer[22, 21]),(AND,Integer[13, 18]),(XOR,Integer[24, 13]),
        (AND,Integer[15, 16]),(XOR,Integer[18, 13]),(AND,Integer[10, 26]),(AND,Integer[17, 12]),
        (XOR,Integer[29, 30]),(AND,Integer[8, 8]),(AND,Integer[28, 31]),(XOR,Integer[14, 33]),(XOR,Integer[12, 16]),
        (XOR,Integer[25, 27]),(XOR,Integer[31, 28]),(XOR,Integer[15, 17]),(AND,Integer[32, 32]),
        (XOR,Integer[21, 22]),(XOR,Integer[34, 36])),
       (36,37,38,39,40,41)) 
c3m2 = build_chromosome( # Second saved 3-bit multiplier circuit 
         (1,2,3,4,5,6),
         ((AND,Integer[3, 3]),(AND,Integer[7, 5]),(AND,Integer[2, 4]),(AND,Integer[1, 4]),
          (AND,Integer[8, 9]),(AND,Integer[5, 2]),(AND,Integer[1, 5]),(AND,Integer[2, 6]),
          (XOR,Integer[12, 11]),(AND,Integer[9, 13]),(AND,Integer[14, 8]),(AND,Integer[7, 5]),
          (AND,Integer[6, 1]),(AND,Integer[19, 12]),(XOR,Integer[4, 15]),(XOR,Integer[9, 13]),
          (XOR,Integer[15, 10]),(XOR,Integer[17, 4]),(XOR,Integer[11, 22]),(AND,Integer[6, 7]),
          (XOR,Integer[20, 23]),(XOR,Integer[27, 12]),(XOR,Integer[12, 19]),(AND,Integer[26, 27]),
          (AND,Integer[7, 24]),(AND,Integer[25, 30]),(XOR,Integer[30, 20]),(AND,Integer[16, 12]),
          (XOR,Integer[32, 21]),(XOR,Integer[14, 18]),(XOR,Integer[34, 32]),(AND,Integer[26, 26]),
          (XOR,Integer[25, 33]),(AND,Integer[28, 35]),(XOR,Integer[29, 31])),
         (36,37,38,39,40,41))
rnw1 = random_neutral_walk( c3m1, [g3m], funcs, 20000, num_mutations=2 )
rnw1 = random_neutral_walk( c3m1, [g3m], funcs, 20000, num_mutations=2, reduce_numactive_reps=500 )
