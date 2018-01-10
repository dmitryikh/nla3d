#encoding: utf-8
import nla3d as nla

nodeTable = [ [0.0,  0.0,  0.0]
            , [16.0, 0.0,  0.0]
            , [16.0, 12.0, 0.0]
            , [32.0, 0.0,  0.0]
            , [32.0, 12.0, 0.0]
            ]
nnodes = len(nodeTable)

elementTable = [ [1, 3]
               , [1, 2]
               , [2, 3]
               , [3, 5]
               , [3, 4]
               , [2, 5]
               , [2, 4]
               , [4, 5]
               ];
nelements = len(elementTable)

storage = nla.FEStorage()
solver = nla.LinearFESolver()

for p in nodeTable:
    node = nla.Node()
    node.pos[0], node.pos[1], node.pos[2] = p
    node.pos *= 12.0
    storage.addNode(node)

for nodes in elementTable:
  el = nla.ElementTRUSS3();
  el.E = 3.0e4;
  el.A = 10.0;
  el.setNodeNumber(0, nodes[0])
  el.setNodeNumber(1, nodes[1])
  storage.addElement(el);

for i in range(1, nnodes + 1):
    solver.addFix(i, nla.Dof.UZ)

solver.addFix(1, nla.Dof.UX);
solver.addFix(1, nla.Dof.UY);
solver.addFix(4, nla.Dof.UX);
solver.addFix(4, nla.Dof.UY);

solver.addLoad(2, nla.Dof.UY, -100.0);
solver.addLoad(5, nla.Dof.UX, 50.0);

solver.attachFEStorage(storage);

solver.solve()

print "DoF solution:"
for i in range(1, nnodes + 1):
    print("%d: UX = %f" % (i, storage.getNodeDofSolution(i, nla.Dof.UX)))
    print("%d: UY = %f" % (i, storage.getNodeDofSolution(i, nla.Dof.UY)))


def check(v1, v2, eps=1.0e-7):
    assert(abs(v1 - v2) < eps)

# check for correct result
check(storage.getNodeDofSolution(2, nla.Dof.UX), 0.0146067)
check(storage.getNodeDofSolution(2, nla.Dof.UY), -0.1046405)

check(storage.getNodeDofSolution(3, nla.Dof.UX), 0.0027214)
check(storage.getNodeDofSolution(3, nla.Dof.UY), -0.0730729)

check(storage.getNodeDofSolution(5, nla.Dof.UX), 0.0055080)
check(storage.getNodeDofSolution(5, nla.Dof.UY), -0.0164325)
