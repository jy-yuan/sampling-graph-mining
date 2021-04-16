from z3 import *

for k in range(10):
    edges_a = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]
    edges_b = [[1, 4], [4, 2], [2, 5], [5, 3], [3, 0], [0, 1]]
    
    E_a = 6
    E_b = 6
    V_a = 6
    V_b = 6
    
    assert E_a == E_b
    assert V_a == V_b
    
    s = Solver()
    f = Function('f', IntSort(), IntSort())
    
    u = Int('u')
    v = Int('v')
    
    s.add(ForAll(u, Implies(And(u >= 0, u < V_a), And(f(u) >= 0, f(u) < V_b))))
    s.add(ForAll((u, v), Implies(And(u != v, u >= 0, u < V_a, v >= 0, v < V_a), f(u) != f(v))))
    
    for i in range(E_a):
        predicate = And(f(edges_a[i][0]) == edges_b[0][0], f(edges_a[i][1]) == edges_b[0][1])
        for j in range(1, E_b):
            predicate = Or(predicate, And(f(edges_a[i][0]) == edges_b[j][0], f(edges_a[i][1]) == edges_b[j][1]))
        s.add(predicate)
    
    print s.check()
    #print s.model()
    
    #for i in range(V_a):
    #    print i, "->", s.model().eval(f(i))
