import numpy as np
from collections.abc import Iterable

## Python functions for the sigmoidal collapse model
## Written by Alex Holehouse (alex.holehouse@wustl.edu)
##
## --> From  Moses, D., Yu, F., Ginell, G., Shamoon, N.M., Koenig, P.S., 
##     Holehouse, A.S., and Sukenik, S. (2020). Probing the Hidden 
##     Sensitivity of Intrinsically Disordered Proteins to their Chemical Environment.
##     
##

# RG prefactor defines some arbitrary lengthscale. We make the  
#RG_PREFACTOR = 0.5
#RE_PREFACTOR = RG_PREFACTOR*np.sqrt(6)

#PREFACTOR = RE_PREFACTOR

RG_PREFACTOR = 0.5
RE_PREFACTOR = RG_PREFACTOR*np.sqrt(6)




# ..........................................................................................
#

def nu2rg(nu, N, PREFACTOR):
    """
    Function that converts an arbitray nu value to an Rg. Note here absolute Rg depends on the
    provided prefactor value.

    Computes Rg with the classic polymeric scaling equation
    

    Rg [or Re] = PREFACTOR * N^{nu}

    Parameters
    ------------
    nu : float
        Apparent scaling exponent. Should be between 0.33 and 0.60 in theory, although no bounds
        checking is done

    N : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Default value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    Returns
    --------
    float 
        Returns the radius of gyration in units determined by the PREFACTOR

    """
    return PREFACTOR*np.power(N, nu)



# ..........................................................................................
#
def nu2chi(nu, N):    
    """
    Function that converts an arbitray nu value to a chi values. 


    Parameters
    ------------
    nu : float
        Apparent scaling exponent. Should be between 0.33 and 0.60 in theory, although no bounds
        checking is done

    N : int
        Chain length

    Returns
    --------
    float 
        Returns the chi value

    """


    rg_theta = np.power(N, 0.5)
    rg = np.power(N, nu)
    
    return (rg/rg_theta)-1



# ..........................................................................................
#

def rg2nu(rg, N, PREFACTOR):
    """
    Function that calculates the radius of gyration (or end-to-end distance) from the apparent
    scaling exponent nu using the classic polymer scaling relationship

    Rg [or Re] = PREFACTOR*N^{nu}

    Parameters
    ------------
    rg : float
        Radius of gyration (in units that correspond to the prefactor)

    N : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Suggested value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    Returns
    --------
    float 
        Returns the apparent scaling exponent


    """
    
    logrg = np.log(rg)
    log_prefactor = np.log(PREFACTOR)
    logN = np.log(N)
    
    return (logrg - log_prefactor)/logN


# ..........................................................................................
#

def chi2nu(chi, N, PREFACTOR):
    """
    Function that calculates the apparent scaling exponeny (nu) from the chi value

    Parameters
    ------------
    chi : float
        Chi value

    N : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Suggested value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    Returns
    --------
    float 
        Returns the apparent scaling exponent


    """

    rg = chi2rg(chi, N, PREFACTOR)
    return rg2nu(rg, N, PREFACTOR)



# ..........................................................................................
#
def chi2rg(chi, N, PREFACTOR):
    """
    Function that converts an arbitray chi value to an Rg. Note here absolute Rg depends on the
    provided prefactor value.

    Computes Rg with the classic polymeric scaling equation
    
    Rg [or Re] = PREFACTOR * N^{nu}

    Where nu is then used to calculate chi

    Parameters
    ------------
    chi : float
        Chi value 
        
    N : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Default value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    Returns
    --------
    float 
        Returns the radius of gyration in units determined by the PREFACTOR
    """



    return (chi+1)*PREFACTOR*np.power(N, 0.50)



# ..........................................................................................
#

def energytochi(energy, L, PREFACTOR, theta=None, midpoint=None):
    """
    Function that converts an inter-bead energy into a chi value in a length-dependent
    manner. Note that if theta and midpoint are not provided the function uses values parameterized
    from a homopolymer cubic lattice model which should be pretty general. Alternatively an explicit
    theta (cooperativity) and midpoint can be provided.

    Parameters
    -------------
    energy : float or iterable
         Float value where 0 is no interaction any anything less than zero is an attractive interaction.
         Anything above zero is treated as the same as zero (i.e. no repulsive interaction). This is a feature
         that could be updated in the future....

    L : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Default value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    theta : float
       Parameter that determines cooperativity of the transition. Bigger the number more highly
       cooperative the transition. Default = None, which means a length-dependent value of Theta
       is calculated based on parameters obtained by fitting against PIMMS simulation data

    midpoint : float
       Parameter that deterimines the midpoint of the transition. Should be between 0 and -infinity.
       Default = None, which means a length-dependent value of Theta
       is calculated based on parameters obtained by fitting against PIMMS simulation data


    Returns
    --------------
    float 
       Returns a float  that corresponds to the chi value for the polymer


    """
    

    def inner(e, L, theta, midpoint):
        base  = nu2chi(0.33, L, PREFACTOR)        
        upper  = nu2chi(0.59, L, PREFACTOR) - base

    
        
        if e >= 0:        
            return base + upper

        # parameters derived from fitting to simulations
        a = 0.6077
        b = -0.1800
        c = -0.274
    
        if theta is None:
            theta = -np.log(L)*a + b
        if midpoint is None:        
            midpoint = c*np.power(L,c)
        

        # goes to 0 and energy is large and negative (compact) and 1 when energy
        # is negligable (expanded)
        energy_scalar = ( 1/(1 + np.power(midpoint/e,theta)))
    
        return base + upper*energy_scalar
    
    
    if isinstance(energy, Iterable):
        return_list = []
        for x in energy:
            return_list.append(inner(x, L, theta, midpoint))
        return return_list
    else:
        return inner(energy, L, theta, midpoint)



# ..........................................................................................
#

def chitoenergy(chi,L, PREFACTOR, theta=None, midpoint=None):
    """
    Function that converts a chi value to an  inter-bead energy in a length-dependent
    manner. Note that if theta and midpoint are not provided the function uses values parameterized
    from a homopolymer cubic lattice model which should be pretty general. Alternatively an explicit
    theta (cooperativity) and midpoint can be provided.

    Parameters
    -------------
    chi : float or iterable of floats
         Float value corresponds to the chi value of the sequence

    L : int
        Chain length

    PREFACTOR : float
        Prefactor used with polymeric calculations. Default value is 0.5 which corresponds to
        working with Rg values for a PIMMS simulations and returns Rg in terms of lattice sites.

    theta : float
       Parameter that determines cooperativity of the transition. Bigger the number more highly
       cooperative the transition. Default = None, which means a length-dependent value of Theta
       is calculated based on parameters obtained by fitting against PIMMS simulation data

    midpoint : float
       Parameter that deterimines the midpoint of the transition. Should be between 0 and -infinity.
       Default = None, which means a length-dependent value of Theta
       is calculated based on parameters obtained by fitting against PIMMS simulation data


    Returns
    --------------
    float 
       Returns a float value that corresponds to the inter-monomer energy for the polymer


    """
    
    a = 0.6077
    b = -0.18
    c = -0.274
    
    if theta is None:
        theta = -np.log(L)*a + b

    if midpoint is None:
        midpoint = c*np.power(L,c)

    base  = nu2chi(0.33,L, PREFACTOR)
    upper  = nu2chi(0.59,L, PREFACTOR) - base
        
    energy_scalar = (chi - base)/upper

    if energy_scalar >= 1.0:
        return 0.0
    
    e = midpoint/np.exp((np.log(1/energy_scalar - 1 ))/theta)
    return e
