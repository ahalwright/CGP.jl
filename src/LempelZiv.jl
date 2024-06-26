#!/usr/bin/env julia
# -*- coding: utf-8 -*-
"""Lempel-Ziv complexity for a binary sequence, in naive Julia code.

- How to use it? From Julia, it's easy:

julia> using LempelZiv
julia> s = "1001111011000010"
julia> LempelZiv.lempel_ziv_complexity(s)    # 1 / 0 / 01 / 11 / 10 / 110 / 00 / 010
8


- Run this .jl file with an argument "test" to run a small benchmark (1000 strings of size 10000).
- Note: there is also a Python version, if you need.

- MIT Licensed, (C) 2017 Lilian Besson (Naereen)
  https://GitHub.com/Naereen/Lempel-Ziv_Complexity
"""

__author__ = "Lilian Besson (Naereen)"
__version__ = "0.2"


module LempelZiv
export lempel_ziv_complexity

"""Lempel-Ziv complexity for a binary sequence, in naive Julia code.

- How to use it? From Julia, it's easy:

>>> using LempelZiv
>>> s = "1001111011000010"
>>> LempelZiv.lempel_ziv_complexity(s)  # 1 / 0 / 01 / 11 / 10 / 110 / 00 / 010
8

- MIT Licensed, (C) 2017 Lilian Besson (Naereen)
  https://GitHub.com/Naereen/Lempel-Ziv_Complexity
"""
function lempel_ziv_complexity(sequence)
    sub_strings = Set()
    n = length(sequence)

    ind = 1
    inc = 0
    while true
        if ind + inc > n
            break
        end
        sub_str = sequence[ind : ind + inc]
        if sub_str in sub_strings
            inc += 1
        else
            push!(sub_strings, sub_str)
            ind += (inc+1)
            inc = 0
        end
    end
    return length(sub_strings)
end

end


if "test" in ARGS
    # import LempelZiv
    lzs = "1001111011000010"
    # LempelZiv.lempel_ziv_complexity(lzs)  # 1 / 0 / 01 / 11 / 10 / 110 / 00 / 010
    println("For s = ", s, " its Lempel-Ziv Complexity is = ", LempelZiv.lempel_ziv_complexity(s))

    M = 1000;
    N = 5000;
    for _ in 1:M
        lzs = join(rand(0:1, rand(N:10*N)));
        println("For a random string lzs of size = ", length(lzs), " its Lempel-Ziv Complexity is = ", LempelZiv.lempel_ziv_complexity(lzs))
        println(@time LempelZiv.lempel_ziv_complexity(lzs))
    end
end
