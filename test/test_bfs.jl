function test_bfs()
  adj = [
    [3,6],#1    
    [4,5],#2    
    [3,5],#3    
    [5,6],#4    
    [1,4],#5    
    [8,7],#6    
    [5,8],#7    
    [4] #8    
  ]
  prev_dict = Dict{Int64,Tuple{Int64,Int64,Bool}}()  # the predecessor and its plength
  queue = Queue{Int64}()    
  start = 1
  goal = 2
  prev_dict[start] = (start,0,false)
  current_plen = 0
  enqueue!(queue,start)
  println("initial queue: ",queue)
  maxrep = 50
  rep = 1
  while length(queue) > 0 && rep <= maxrep
    current = dequeue!(queue)
    println("new current: ",current)
    (prev,plength,done) = prev_dict[current]
    if done 
      println("done inf loop:  current: ",current)
      return nothing
    end
    println("current: ",current,"  adj[current]: ",adj[current])
    for nxt in adj[current]
      println("(nxt,current): ",(nxt,current))
      if haskey( prev_dict, nxt )
        (pred,plen,done) = prev_dict[nxt]
        if plen > plength+1
          prev_dict[nxt] = (current,plength+1,true)
        else
          prev_dict[nxt] = (pred,plength+1,true)
        end
      else
        prev_dict[nxt] = (current,plength+1,false)
      end
      println("prev_dict: ",prev_dict)
      if nxt == goal
        reverse_path = [nxt]
        while nxt != start
          (prev,plength,done) = prev_dict[nxt]
          println("done (prev,plength): ",(prev,plength))
          push!(reverse_path,prev)
          nxt = prev
        end
        return reverse_path
      end
      enqueue!(queue,nxt)
      println("new queue: ",queue)
    end  # for
    rep += 1
  end  # while
end
