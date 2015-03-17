import numpy as np
pdi = self.GetInput()
pdo = self.GetOutput()

nn = pdi.GetNumberOfPoints()
en = pdi.GetNumberOfPoints()

disp = vtk.vtkFloatArray()
disp.SetName('u')
disp.SetNumberOfComponents(3)
disp.SetNumberOfTuples(nn)
for i in range(nn):
  ux = pdi.GetPointData().GetArray('ux').GetTuple1(i)
  uy = pdi.GetPointData().GetArray('uy').GetTuple1(i)
  uz = pdi.GetPointData().GetArray('uz').GetTuple1(i)
  disp.SetTuple(i, (ux, uy, uz))

pdo.GetPointData().AddArray(disp)

Tdir1 = vtk.vtkFloatArray()
Tdir1.SetName('Tdir1')
Tdir1.SetNumberOfComponents(3)
Tdir1.SetNumberOfTuples(en)

Tdir2 = vtk.vtkFloatArray()
Tdir2.SetName('Tdir2')
Tdir2.SetNumberOfComponents(3)
Tdir2.SetNumberOfTuples(en)

Tdir3 = vtk.vtkFloatArray()
Tdir3.SetName('Tdir3')
Tdir3.SetNumberOfComponents(3)
Tdir3.SetNumberOfTuples(en)


T1 = vtk.vtkFloatArray()
T1.SetName('T1')
T1.SetNumberOfComponents(1)
T1.SetNumberOfTuples(en)

T2 = vtk.vtkFloatArray()
T2.SetName('T2')
T2.SetNumberOfComponents(1)
T2.SetNumberOfTuples(en)

T3 = vtk.vtkFloatArray()
T3.SetName('T3')
T3.SetNumberOfComponents(1)
T3.SetNumberOfTuples(en)

for i in range(en):
    txy = pdi.GetCellData().GetArray('TXY').GetTuple1(i)
    txx = pdi.GetCellData().GetArray('TXX').GetTuple1(i)
    txz = pdi.GetCellData().GetArray('TXZ').GetTuple1(i)
    tyy = pdi.GetCellData().GetArray('TYY').GetTuple1(i)
    tyz = pdi.GetCellData().GetArray('TYZ').GetTuple1(i)
    tzz = pdi.GetCellData().GetArray('TZZ').GetTuple1(i)
    T = np.array([[txx,txy,txz],[txy,tyy,tyz],[txz,tyz,tzz]])
    try:
        vals, dirs = np.linalg.eigh(T,'U')
    except:
        T1.SetTuple1(i, 0.0)
        T2.SetTuple1(i, 0.0)
        T3.SetTuple1(i, 0.0)
        Tdir1.SetTuple(i, [0.0, 0.0, 0.0])
        Tdir2.SetTuple(i, [0.0, 0.0, 0.0])
        Tdir3.SetTuple(i, [0.0, 0.0, 0.0])
        continue

    if (vals[0] >= vals[1] and vals[0] >= vals[2]):
        T1.SetTuple1(i, vals[0])
        Tdir1.SetTuple(i, dirs[:,0])
        if (vals[1] >= vals[2]):
            T2.SetTuple1(i, vals[1])
            Tdir2.SetTuple(i, dirs[:,1])
            T3.SetTuple1(i, vals[2])
            Tdir3.SetTuple(i, dirs[:,2])
        else:
            T2.SetTuple1(i, vals[2])
            Tdir2.SetTuple(i, dirs[:,2])
            T3.SetTuple1(i, vals[1])
            Tdir3.SetTuple(i, dirs[:,1])
    elif (vals[1] >= vals[0] and vals[1] >= vals[2]):
        T1.SetTuple1(i, vals[1])
        Tdir1.SetTuple(i, dirs[:,1])
        if (vals[0] >= vals[2]):
            T2.SetTuple1(i, vals[0])
            Tdir2.SetTuple(i, dirs[:,0])
            T3.SetTuple1(i, vals[2])
            Tdir3.SetTuple(i, dirs[:,2])
        else:
            T2.SetTuple1(i, vals[2])
            Tdir2.SetTuple(i, dirs[:,2])
            T3.SetTuple1(i, vals[0])
            Tdir3.SetTuple(i, dirs[:,0])
    else:
        T1.SetTuple1(i, vals[2])
        Tdir1.SetTuple(i, dirs[:,2])
        if (vals[0] >= vals[1]):
            T2.SetTuple1(i, vals[0])
            Tdir2.SetTuple(i, dirs[:,0])
            T3.SetTuple1(i, vals[1])
            Tdir3.SetTuple(i, dirs[:,1])
        else:
            T2.SetTuple1(i, vals[1])
            Tdir2.SetTuple(i, dirs[:,1])
            T3.SetTuple1(i, vals[0])
            Tdir3.SetTuple(i, dirs[:,0])


pdo.GetCellData().AddArray(T1)
pdo.GetCellData().AddArray(T2)
pdo.GetCellData().AddArray(T3)
pdo.GetCellData().AddArray(Tdir1)
pdo.GetCellData().AddArray(Tdir2)
pdo.GetCellData().AddArray(Tdir3)




lambda1 = vtk.vtkFloatArray()
lambda1.SetName('lambda1')
lambda1.SetNumberOfComponents(1)
lambda1.SetNumberOfTuples(en)

lambda2 = vtk.vtkFloatArray()
lambda2.SetName('lambda2')
lambda2.SetNumberOfComponents(1)
lambda2.SetNumberOfTuples(en)

lambda3 = vtk.vtkFloatArray()
lambda3.SetName('lambda3')
lambda3.SetNumberOfComponents(1)
lambda3.SetNumberOfTuples(en)

for i in range(en):
    txy = pdi.GetCellData().GetArray('EXY').GetTuple1(i)
    txx = pdi.GetCellData().GetArray('EXX').GetTuple1(i)
    txz = pdi.GetCellData().GetArray('EXZ').GetTuple1(i)
    tyy = pdi.GetCellData().GetArray('EYY').GetTuple1(i)
    tyz = pdi.GetCellData().GetArray('EYZ').GetTuple1(i)
    tzz = pdi.GetCellData().GetArray('EZZ').GetTuple1(i)
    T = np.array([[txx,txy,txz],[txy,tyy,tyz],[txz,tyz,tzz]])
    try:
        vals, dirs = np.linalg.eigh(T,'U')
    except:
        lambda1.SetTuple1(i, 0.0)
        lambda2.SetTuple1(i, 0.0)
        lambda3.SetTuple1(i, 0.0)
        continue

    if (vals[0] >= vals[1] and vals[0] >= vals[2]):
        lambda1.SetTuple1(i, vals[0])
        if (vals[1] >= vals[2]):
            lambda2.SetTuple1(i, vals[1])
            lambda3.SetTuple1(i, vals[2])
        else:
            lambda2.SetTuple1(i, vals[2])
            lambda3.SetTuple1(i, vals[1])
    elif (vals[1] >= vals[0] and vals[1] >= vals[2]):
        lambda1.SetTuple1(i, vals[1])
        if (vals[0] >= vals[2]):
            lambda2.SetTuple1(i, vals[0])
            lambda3.SetTuple1(i, vals[2])
        else:
            lambda2.SetTuple1(i, vals[2])
            lambda3.SetTuple1(i, vals[0])
    else:
        lambda1.SetTuple1(i, vals[2])
        if (vals[0] >= vals[1]):
            lambda2.SetTuple1(i, vals[0])
            lambda3.SetTuple1(i, vals[1])
        else:
            lambda2.SetTuple1(i, vals[1])
            lambda3.SetTuple1(i, vals[0])


pdo.GetCellData().AddArray(lambda1)
pdo.GetCellData().AddArray(lambda2)
pdo.GetCellData().AddArray(lambda3)
