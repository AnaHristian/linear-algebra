from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        temp_1 = s[row1]
        temp_2 = s[row2]
        s[row1] = temp_2
        s[row2] = temp_1
        return LinearSystem(s)


    def multiply_coefficient_and_row(self, coefficient, row):
        pass # add your code here


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        pass # add your code here


    def indices_of_first_nonzero_terms_in_each_row(self):
        """
        which term of each equation is the first non-zero term
        usefull for finding the pivot var in each eq
        """
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p0 = Plane(normal_vector=['1','1','1'], constant_term='1')
p1 = Plane(normal_vector=['0','1','0'], constant_term='2')
p2 = Plane(normal_vector=['1','1','-1'], constant_term='3')
p3 = Plane(normal_vector=['1','0','-2'], constant_term='2')

s = LinearSystem([p0,p1,p2,p3])

#print('{},{},{},{}'.format(s[0],s[1],s[2],s[3]))
#print(len(s))
#print(s)
#print(s.swap_rows(1,3))
#print(s.swap_rows(1,2))
#s[0] = p0
#print(s[0] == p1)
#p0 = s[1]
#p1 = s[0]
#s[0] = p1
#s[1] = p0
print(s[1])
s = LinearSystem([p1,p0,p2,p3])
print(s[1])
#print(s)
#print(s[0] == p1)
#print(s[0])

#print(p1)

#s.swap_rows(0,1)
#print(s[0])
#print(type(s.planes[0]))
#print(type(s[0]))
#print(s[0] == p1)
#print(s.planes[0].normal_vector)
#print(s.planes[0].constant_term)

#if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#    print('test case 1 failed')


#print(p0.normal_vector)
#print(s.indices_of_first_nonzero_terms_in_each_row())
#s[0] = p1
#print(s)
#print(MyDecimal('1e-9').is_near_zero())
#print(MyDecimal('1e-11').is_near_zero())