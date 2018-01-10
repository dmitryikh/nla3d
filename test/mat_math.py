#/bin/env python2.7
#encoding: utf-8

import nla3d as nla

# test1
mat = nla.Mat33()
mat.from_list([ [1.0, 2.0, 3.0]
              , [4.0, 5.0, 6.0]
              , [7.0, 8.0, 9.0]
              ])
print mat
mat2 = mat.copy()
assert(mat2.compare(mat))


# test2
mat2.zero()
zero = nla.Mat33()
assert(zero.compare(mat2))


# test3
mat = nla.Mat33()
mat.from_list([ [1.0, 2.0, 3.0]
              , [43.0, 5.0, 6.0]
              , [7.0, 8.0, 20.0]
              ])
assert(mat.det() == -657)

mat_inv = nla.Mat33()
mat_inv.from_list([ [-0.07914764,  0.02435312,  0.00456621]
                  , [ 1.24505327,  0.00152207, -0.18721461]
                  , [-0.47031963, -0.00913242,  0.12328767]
                  ])
assert(mat_inv.compare(mat.inv(mat.det())))


# test4
sym3 = nla.MatSym3()
sym3.zero()
sym3[0, 0] = 5.0
sym3[0, 1] = -35.0
sym3[0, 2] = 12.0
sym3[1, 1] = 15.0
sym3[1, 2] = -3.0
sym3[2, 2] = 14.0

print sym3
print sym3.toMat()
