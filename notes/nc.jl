function neutral_component( ch::Chromosome, funcs::Vector{Func} )::Set
  p = ch.params
  ov = output_values(ch)[1]
  nc = Set([chromosome_to_int(ch,funcs)])
  stack = [chromosome_to_int(ch,funcs)]
  while length(stack) > 0
    circ = int_to_chromosome( pop!(stack), p, funcs )
    circs = filter(cch->ov==output_values(cch)[1], mutate_all( circ, funcs, output_outputs=false, output_circuits=true ))
    for c in circs
      if !( chromosome_to_int(c) in nc )
        push!(nc,chromosome_to_int(c,funcs))
        push!(stack,chromosome_to_int(c,funcs))
      end
    end
  end
  nc
end
