#encoding: utf-8
import nla3d as nla
import sys


def readTempData(filename):
  temps = []
  with open(filename, 'r') as fin:
    fin.readline()
    for line in fin:
      temps.append(float(line.strip().split()[1]))
  return temps


def check(v1, v2, eps=1.0e-7):
  assert(abs(v1 - v2) < eps)


def main():
  cdb_filename = sys.argv[1];
  res_filename = sys.argv[2] if len(sys.argv) > 2 else None;

  md = nla.MeshData();
  if not nla.readCdbFile(cdb_filename, md):
    raise RuntimeError("Can't read FE info from %s file" % cdb_filename)

  md.compressNumbers();

  storage = nla.FEStorage()
  solver = nla.LinearFESolver()
  solver.attachFEStorage(storage);

  # add nodes
  sind = storage.createNodes(md.nodesNumbers.size())
  for i, j in enumerate(sind):
    storage.getNode(j).pos = md.nodesPos[i] # NOTE: we need copy by value here

  # add elements
  ind = md.getCellsByAttribute("TYPE", 1)
  sind = storage.createElements(ind.size(), nla.ElementType_QUADTH)
  for i, j in zip(sind, ind):
    el = storage.getElementQUADTH(i)
    for k in range(el.getNNodes()):
      el.setNodeNumber(k, md.cellNodes[j][k])
    el.k = 0.018

  # add surface elements
  ind = md.getCellsByAttribute("TYPE", 2)
  sind = storage.createElements(ind.size(), nla.ElementType_SurfaceLINETH)
  for i, j in zip(sind, ind):
    el = storage.getElementSurfaceLINETH(i)
    for k in range(el.getNNodes()):
      el.setNodeNumber(k, md.cellNodes[j][k])
    el.htc = 0.0034
    el.etemp[0] = -6.0
    el.etemp[1] = -6.0

  #add loadBc
  for v in md.loadBcs:
    solver.addLoad(v.node, v.node_dof, v.value)

  # add fixBc
  for v in md.fixBcs:
    solver.addFix(v.node, v.node_dof, v.value)

  vtk = nla.VtkProcessor(storage, "QUADTHpy");
  vtk.writeAllResults();
  solver.addPostProcessor(vtk);

  solver.solve();

  # check temperature results with Ansys data
  if res_filename:
    print('Checking temperatures..')
    ans_temps = readTempData(res_filename)
    for i, v in zip(range(1, storage.nNodes() + 1), ans_temps):
      check(storage.getNodeDofSolution(i, nla.Dof.TEMP), v, 1.0e-3)
    print('OK!')


if __name__ == '__main__':
    main()
