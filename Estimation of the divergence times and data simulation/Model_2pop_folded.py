import numpy
import moments
from moments import Numerics, Manips, Integration
from moments.Spectrum_mod import Spectrum


######################
####  SI  and SI+ ####
######################

def SI(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	T1= time of split
    """
    nu1,nu2,T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    return fs
    
def SI_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
    nu1b=pop size for pop1 at T2
	nu2b=pop size for pop2 at T2
	T1= time of split
	T2= time of change in NE
	T2= time of change in NE
    """
    nu1,nu2,nu1b,nu2b,T1,T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    return fs
    


######################
####  IM  and IM+ ####
######################
    
def IM(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	T1= time of split
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,T1,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs


    
def IM_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
    nu1b=pop size for pop1 at T2
	nu2b=pop size for pop2 at T2
	T1= time of split
	T2= time of change in NE
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs
######################
####  SC and SC + ####
######################

def SC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	T1= time of population split
	T2= time of secondary contact
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs

def SC_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	nu1b=pop size for pop1 at T2
	nu2b=pop size for pop2 at T2
	T1= time of population split
	T2= time of secondary contact
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs


######################
####   AM and AM +####
######################

def AM(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	T1= time of population split with migration
	T2= time of speciaton
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21,0 ]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs


def AM_NeC(params, ns):
    """
    nu1= pop size for pop1
	nu2=pop size for pop2 
	nu1b=pop size for pop1 at T2
	nu2b=pop size for pop2 at T2
	T1= time of population split with migration
	T2= time of speciaton
	m12= proportion of pop1 made up of migrants from pop 2
	m21= proportion of pop2 made up of migrants from pop 1
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12], [m21,0 ]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs


######################
####   2EP + 2EP+ ####
######################


def two_ep(params, ns):
    """
    nu1= pop size for pop0
	nu2=pop size for pop1 
	T1= time of population split
	T2= time of second migration epoch
	m12_0= migration rate from pop0 to pop2
	m21_0= migration rate from pop1 to pop0
	m12= migration rate from pop0 to pop2
	m21= migration rate from pop1 to pop0
    """
    nu1,nu2,T1,T2,m12_0,m21_0,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs

def two_ep_NeC(params, ns):
    """
    nu1= pop size for pop0
	nu2=pop size for pop1 
 nu1b=pop size for pop1 at T2
	nu2b=pop size for pop2 at T2
	T1= time of population split
	T2= time of second migration epoch
	m12_0= migration rate from pop0 to pop2
	m21_0= migration rate from pop1 to pop0
	m12= migration rate from pop0 to pop2
	m21= migration rate from pop1 to pop0
    """
    nu1,nu2,nu1b,nu2b,T1,T2,m12_0,m21_0,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m = numpy.array([[0, m12_0], [m21_0, 0]]))
    fs.integrate([nu1b, nu2b], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs
