export execute_chromosome
using Printf

function apply( f::Function,args)
  if length(args) == 0
    f()
  elseif length(args) == 1
    f(args[1])
  elseif length(args) == 2
    f(args[1],args[2])
  elseif length(args) == 3
    f(args[1],args[2],args[3])
  end
end

function evaluate_node(c::Chromosome, node::InputNode, context::Vector)
    #println("eval node  node: ",node)
    #print("context[node.index]:")
    #Printf.@printf("0x0%4x\n",context[node.index])
    return context[node.index]
end

function evaluate_node(c::Chromosome, node::InteriorNode, context::Vector)
    #println("eval node  node: ",node)
    func = node.func
    args = map(node.inputs[1:func.arity]) do position
        index = position
        evaluate_node(c, c[index], context)
    end
    result = apply(func.func, args) & c.params.mask  # Trim unnecessary bits
    #print("func: ",func,"  args: ",args,"  result: ")
    #Printf.@printf("0x0%4x\n",result)
    #return apply(func.func, args) & c.params.mask  # Trim unnecessary bits
    return result
end

function evaluate_node(c::Chromosome, node::OutputNode, context::Vector)
    index = node.input
    return evaluate_node(c, c[index], context)
end

function execute_chromosome(c::Chromosome, context::Vector)
    return [evaluate_node(c, node, context) for node = c.outputs]
end
