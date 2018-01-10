#encoding: utf-8
import nla3d as nla
import sys


def readDispData(filename):
  disps = []
  with open(filename, 'r') as fin:
    fin.readline()
    for line in fin:
      v = nla.Vec3()
      v[0], v[1], v[2] = map(float, line.strip().split())[1:]
      disps.append(v)
  return disps


def readStressData(filename):
  stresses = []
  with open(filename, 'r') as fin:
    fin.readline()
    for line in fin:
      m = nla.MatSym3()
      m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[1, 2], m[0, 2] = map(float, line.strip().split())[1:]
      stresses.append(m)
  return stresses


def check(v1, v2, eps=1.0e-7):
  assert(abs(v1 - v2) < eps)


def main():
  cdb_filename = sys.argv[1];
  res_disp_filename = sys.argv[2] if len(sys.argv) > 2 else None;
  res_stress_filename = sys.argv[3] if len(sys.argv) > 3 else None;

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
  sind = storage.createElements(ind.size(), nla.ElementType_TETRA0)
  for i, j in zip(sind, ind):
    el = storage.getElementTETRA0(i)
    el.setNodeNumber(0, md.cellNodes[j][0])
    el.setNodeNumber(1, md.cellNodes[j][1])
    el.setNodeNumber(2, md.cellNodes[j][2])
    el.setNodeNumber(3, md.cellNodes[j][4])
    el.E = 1.0e8
    el.my = 0.3

  #add loadBc
  for v in md.loadBcs:
    solver.addLoad(v.node, v.node_dof, v.value)

  # add fixBc
  for v in md.fixBcs:
    solver.addFix(v.node, v.node_dof, v.value)

  # TODO:
  # VtkProcessor* vtk = new VtkProcessor (&storage, "tetra");
  # solver.addPostProcessor(vtk);
  # vtk->writeAllResults();

  solver.solve();

  # check displacements results with Ansys data
  if res_disp_filename:
    print('Checking displacements..')
    ans_disps = readDispData(res_disp_filename)
    for i, v in zip(range(1, storage.nNodes() + 1), ans_disps):
      for j, d in enumerate([nla.Dof.UX, nla.Dof.UY, nla.Dof.UZ]):
        check(storage.getNodeDofSolution(i, d), v[j], 1.0e-10)
    print('OK!')

  # check stress results with Ansys data
  if res_stress_filename:
    print('Checking stresses..')
    ans_stresses = readStressData(res_stress_filename)
    mat = nla.MatSym3()
    for i, v in zip(range(1, storage.nElements() + 1), ans_stresses):
      mat.zero()
      storage.getElement(i).getTensor(mat, nla.tensorQuery_E)
      assert(mat.compare(v, 1e-6))
    print('OK!')


if __name__ == '__main__':
    main()
