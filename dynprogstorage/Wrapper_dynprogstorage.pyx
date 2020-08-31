# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string

# follow example here http://docs.cython.org/en/latest/src/userguide/wrapping_CPlusPlus.html#add-public-attributes
# Declare the class with cdef
cdef extern from "cplfunction.hpp" :
    cdef cppclass cplfunction:
        cplfunction() except +
        cplfunction(cplfunction)
        cplfunction(double,double) except +
        #cplfunction(cplfunction const) except +
        cplfunction(vector[double],vector[double],double) except +
        #cplfunction operator=(cplfunction)
        double evalf(double )
        void AddSimple(double , double , double , double )
        void print_()
        map[double,double] get_BreakPoints()
        double FirstBreakVal_
        double FirstSlopeVal_
        map[double, double] Breakpoints_



#Create Cython wrapper class http://docs.cython.org/en/latest/src/userguide/wrapping_CPlusPlus.html#create-cython-wrapper-class
# le guide cython ne fonctionne pas
# another inspiration https://stackoverflow.com/questions/33677231/how-to-expose-a-function-returning-a-c-object-to-python-without-copying-the-ob


#from cplfunction cimport cplfunction

cdef class Pycplfunction:
    cdef cplfunction thisptr ## * is required here see https://stackoverflow.com/questions/33573038/how-to-expose-a-function-returning-a-c-object-to-python-using-cython
    # typiquement le genre de bazar que l'on n'aurait pas si l'on pouvait utiliser boost-python ;)
    # would it be possible to use a switch case ?
    def __cinit__(self,vector[double] Slopes,vector[double] BreakPoints, double FirstBreakVal):
        self.thisptr = cplfunction(Slopes,BreakPoints,FirstBreakVal)

    cdef copy(self, cplfunction other):
        self.thisptr = cplfunction(other)

    def evalf(self,double x):
        return self.thisptr.evalf(x)

    def getBreakPoints(self):
        return self.thisptr.get_BreakPoints()

    def AddSimple(self,double val1,double val2,double val3,double val4):
        return self.thisptr.AddSimple(val1,val2,val3,val4)

    def Myprint(self):
        return self.thisptr.print_()

    def __repr__(self):
        return "<cplfunction: Breakpoints_={}, FirstBreakVal_={0}, FirstSlopeVal_={0}>".format(self.Breakpoints_,self.FirstBreakVal_, self.FirstSlopeVal_)

    property Breakpoints_:
        def __get__(self): return self.thisptr.Breakpoints_
        def __set__(self, map[double,double] Breakpoints_): self.thisptr.Breakpoints_ = Breakpoints_

    property FirstBreakVal_:
        def __get__(self): return self.thisptr.FirstBreakVal_
        def __set__(self, FirstBreakVal_): self.thisptr.FirstBreakVal_ = FirstBreakVal_

    property FirstSlopeVal_:
        def __get__(self): return self.thisptr.FirstSlopeVal_
        def __set__(self, FirstSlopeVal_): self.thisptr.FirstSlopeVal_ = FirstSlopeVal_

#    @property
#    def x1(self):
#        return self.c_rect.x1
#    @x1.setter
#    def x1(self, x1):
#        self.c_rect.x1 = x1

## a regarder http://nicolas-hug.com/blog/cython_notes


# distutils: language = c++
# distutils: sources = cplfunction.cpp

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string


# follow example here http://docs.cython.org/en/latest/src/userguide/wrapping_CPlusPlus.html#add-public-attributes
# Declare the class with cdef
cdef extern from "cplfunction.cpp" :
    cdef cppclass cplfunctionvec:
        cplfunctionvec() except +
        cplfunctionvec(int) except +
        cplfunctionvec(cplfunctionvec &) except +
        #cplfunctionvec(cplfunctionvec const) except +
        cplfunctionvec(vector[double],vector[double],vector[double])  except +
        cplfunctionvec(vector[double],vector[double],vector[double],vector[double],vector[double])  except +
        cplfunction vec_get(int)
        int size()
        cplfunctionvec Maxf(cplfunctionvec &)
        void Max_(cplfunctionvec &,cplfunctionvec &)
        vector[cplfunction] MycplfunctionList_
        void Maxf_1Breaks_withO(vector[double],vector[double],vector[double])
        void Maxf_2Breaks_withO(vector[double],vector[double],vector[double],vector[double],vector[double])
        vector[double] EvalDeltaf2(vector[double] )
        vector[double] EvalDeltafMoins(vector[double] )
        vector[double] EvalDeltafPlus(vector[double] )
        vector[double] Evalf(vector[double])
        vector[double] OptimMargInt(vector[double],vector[double],vector[double],vector[double])

#from cplfunctionvec cimport cplfunctionvec

#https://stackoverflow.com/questions/33677231/how-to-expose-a-function-returning-a-c-object-to-python-without-copying-the-ob/33677627#33677627
cdef extern from "<utility>":
    vector[cplfunction]&& move(vector[cplfunction]&&) # just define for peak rather than anything else

cdef class Pycplfunctionvec:
    cdef cplfunctionvec thisptr

#https://stackoverflow.com/questions/18260095/cant-override-init-of-class-from-cython-extension
# Trick https://stackoverflow.com/questions/13201886/cython-and-constructors-of-classes
    #def __cinit__(self, vector[double] S1, vector[double] B1, vector[double] f0):
    def __cinit__(self, S1=None,B1=None,f0=None,S2=None,B2=None):
        if S2 is not None and B2 is not None:
            self.thisptr = cplfunctionvec(S1,S2,B1,B2,f0)
        else:
            self.thisptr = cplfunctionvec(S1,B1,f0)

    #
    # def __dealloc__(self):
    #     del self.thisptr

    cdef copy(self, cplfunctionvec &other):
        self.thisptr = cplfunctionvec(other)

    def size(self):
        return(self.thisptr.size())

    def __repr__(self):
        return "<cplfunctionvec: MycplfunctionList_={}>".format(self.MycplfunctionList_)

#    property MycplfunctionList_:
#        def __get__(self): return self.thisptr.MycplfunctionList_
#        def __set__(self, vector[cplfunction] MycplfunctionList_): self.thisptr.MycplfunctionList_ = MycplfunctionList_

# https://stackoverflow.com/questions/62828472/wrapping-a-class-that-contains-a-vector-of-complex-class-with-cython
# impossible here
    def Maxf_1Breaks_withO(self,vector[double] S1, vector[double] B1, vector[double] f0):
        self.thisptr.Maxf_1Breaks_withO( S1,  B1, f0)

    def Maxf_2Breaks_withO(self,vector[double] S1,vector[double] S2, vector[double] B1, vector[double] B2,vector[double] f0):
        self.thisptr.Maxf_2Breaks_withO( S1, S2, B1,B2, f0)

    def OptimMargInt(self,vector[double]  Pmoins,vector[double]  Pplus,vector[double]  Cmoins,vector[double] Cplus):
        return(self.thisptr.OptimMargInt( Pmoins, Pplus, Cmoins, Cplus))

    def Evalf(self, vector[double] x):
        return(self.thisptr.Evalf(x))

    def EvalDeltafMoins(self,vector[double] x):
        return(self.thisptr.EvalDeltafMoins(x))

    def EvalDeltafPlus(self,vector[double] x):
        return(self.thisptr.EvalDeltafPlus(x))



   # def Max_(self, Pycplfunctionvec x, Pycplfunctionvec y):
   #     self.thisptr.Max_(x.thisptr,y.thisptr)

    #def
 #   def vec_get(self,i):
 #       cdef Pycplfunction Pyres
 #       Pyres = Pycplfunction()
 #       Pyres.copy(self.thisptr.vec_get(i))
 #       return(Pyres)
