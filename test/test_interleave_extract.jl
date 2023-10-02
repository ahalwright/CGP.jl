function test_interleave_extract( P::Parameters, funcs::Vector{Func} )
  ph1 = randgoal(P)
  ph2 = randgoal(P)
  phi = [interleave( 2^P.numinputs, ph1[1], ph2[1] )]
  println("ph1: ",ph1,"  ph2: ",ph2,"  phi: ",phi)
  println( extract_odd_even( phi[1], P.numinputs+1 ) )
  (phe1,phe2) = extract_odd_even( phi[1], P.numinputs+1 )
  println("phe1: ",[phe1],"  phe2: ",[phe2])
  phii = [interleave( 2^P.numinputs, phe2[1], phe1[1] )]
  println("phii: ",phii)
end

