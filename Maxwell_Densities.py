#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 10:59:41 2023

@author: luiz
"""

import matplotlib.pyplot as plt
import numpy as np
from Miscellaneous import CS_EOS, Integration

# We give this name because is the code that compute the densities
class DensityFunc:
    
    # Variables
    P_0 = 1 # Equilibrium Pressure
    T = 1 # Equilibrium Temperature
    
    # Maximum and minimum
    rho_max = 1
    rho_min = 0
    rho_MAX = 1
    rho_MIN = 0
    rho_ROOT = 1
    
    # Pointer to the desire EOS
    pEOS = 0
    
    # Function we want to integrate
    def DensFunc( self, rho ):
        # Defining necessary parameters
        T = self.T
        P_0 = self.P_0
        P_eos = self.pEOS.P_from_rho( rho, T )
        # Computing the function
        with np.errstate(all='raise'): # Avoid division by zero
            try:
                f = ( P_eos - P_0 ) / rho / rho
            except:
                f = ( P_eos - P_0 ) 
        return f
    
    # Derivative of the function we want to integrate
    def dDensFunc( self, rho ):
        # Defining necessary parameters
        T = self.T
        P_0 = self.P_0
        P_eos = self.pEOS.P_from_rho( rho, T )
        dP_eos = self.pEOS.dP_from_rho( rho, T ) 
        # Computing the function derivative
        with np.errstate(all='raise'): # Avoid division by zero
            try:
                df = dP_eos / rho / rho - 2 * ( P_eos - P_0 ) / rho ** 3
            except:
                df = dP_eos 
        return df
    
    """
    # Find the maximum point of the base function
    def MAX_func( self ):
        # Defining necessary parameters   
        b = self.pEOS.b
        rho_c = self.pEOS.rho_c
        # Interval of search
        rho_i = 0
        rho_f = 3 / b
        # Incremental Search
        drho = rho_c / 100
        rho = rho_i
        verifier = 0
        while verifier == 0:
            df = self.dDensFunc( rho )
            rho_next = rho + drho
            df_next = self.dDensFunc( rho_next )
            if df > 0 and df_next < 0:
                rho_min = rho
                rho_max = rho_next
                verifier = verifier + 1
            rho = rho_next
            if rho > rho_f:
                print("Problems to find function maximum")
                break
        # Bissection Method
        tol = 1e-6
        while ( rho_max - rho_min ) > tol:
            rho_med = ( rho_max + rho_min ) / 2
            df_med = self.dDensFunc( rho_med )
            if df_med > 0:
                rho_min = rho_med
            else:
                rho_max = rho_med
        return rho_med
    
    # Find the minimum point of the base function
    def MIN_func( self ):
        # Defining necessary parameters   
        b = self.pEOS.b
        rho_c = self.pEOS.rho_c
        # Interval of search
        rho_i = 0
        rho_f = 3 / b
        # Incremental Search
        drho = rho_c / 100
        rho = rho_i
        verifier = 0
        while verifier == 0:
            df = self.dDensFunc( rho )
            rho_next = rho + drho
            df_next = self.dDensFunc( rho_next )
            if df < 0 and df_next > 0:
                rho_min = rho
                rho_max = rho_next
                verifier = verifier + 1
            rho = rho_next
            if rho > rho_f:
                print("Problems to find function minimum")
                break
        # Bissection Method
        tol = 1e-6
        while ( rho_max - rho_min ) > tol:
            rho_med = ( rho_max + rho_min ) / 2
            df_med = self.dDensFunc( rho_med )
            if df_med < 0:
                rho_min = rho_med
            else:
                rho_max = rho_med
        return rho_med    
    
    # Find the Zero Point of the Base Function
    def ROOT( self ):
        # Importing the necessary classes
        # Defining necessary parameters  
        rho_left = self.rho_MAX
        rho_right = self.rho_MIN
        # Bissection Method
        tol = 1e-6
        while ( rho_right - rho_left ) > tol:
            rho_med = ( rho_left + rho_right ) / 2
            f_med = self.DensFunc( rho_med )
            if f_med > 0:
                rho_left = rho_med
            else:
                rho_right = rho_med
        return rho_med   
    """
#****************************************************************************

# Find the Condition of the Maxwell Rule
def MAXWELL( EOS, Dens, INT ):
    # Defining range of search
    rho_max = Dens.rho_max
    P_MAX = EOS.P_from_rho( rho_max, T )
    P_MIN = 0
    # Bissection Method
    tol = 1e-8
    Int = 1
    while ( P_MAX - P_MIN > tol ) | ( abs(Int) > tol ):
        P_0 = ( P_MAX + P_MIN ) / 2
        # Updating the function pressure
        Dens.P_0 = P_0
        # Updating the integration limits
        INT.x_i = EOS.rhov_from_P( P_0, Dens.T )
        INT.x_f = EOS.rhol_from_P( P_0, Dens.T )
        # Integrating the base function
        Int = INT.Simpson( Dens.DensFunc )
        # Evaluating the results
        if Int > 0:
            P_MIN = P_0
        else:
            P_MAX = P_0
    rho_v = EOS.rhov_from_P( P_0, Dens.T )
    rho_l = EOS.rhol_from_P( P_0, Dens.T )
    return P_0, rho_v, rho_l
    
#****************************************************************************
# Code Start

# Initialization

# Storing EOS information
EOS = CS_EOS()
# Equation of state parameters
EOS.a = 0.25
EOS.b = 4
EOS.R = 1
# Computing reduced properties
T_c = EOS.a / EOS.b / EOS.R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * EOS.a / EOS.b ** 2
rho_c = 0.5218 / EOS.b
EOS.T_c = T_c
EOS.P_c = P_c
EOS.rho_c = rho_c
# Defining the EOS temperature
EOS.T = 0.6 * T_c

# Storing function information
Dens = DensityFunc()
Dens.T = 0.6 * T_c
Dens.pEOS = EOS
#
T = Dens.T
rho_i = 0.08*rho_c
rho_f = 3.2*rho_c

# Defining max and min for the EOS
Dens.rho_max = EOS.MAX_EOS( T ) # point of maximum for EOS
print("rho_max: ", Dens.rho_max)
Dens.rho_min = EOS.MIN_EOS( T ) # point of minimum for EOS
print("rho_min: ", Dens.rho_min)

# Defining some dynamics properties
rho = np.linspace( rho_i, rho_f, 100 )
P_eos = EOS.P_from_rho( rho, T )
Dens.P_0 = EOS.P_from_rho( rho_i, T )
P_0 = Dens.P_0

# Defining integration properties
INT = Integration()
INT.N = 20000

# Obtaining Maxwell densities
P_0, rho_v, rho_l = MAXWELL( EOS, Dens, INT )
print("P_0: ", P_0)
print("rho_v: ", rho_v)
print("rho_l: ", rho_l)

# Verifying if the integral is zero with Retangular rule
Int_REC = INT.Rectangular( Dens.DensFunc )
print("Integral O1: ", Int_REC)

# Verifying if the integral is zero with Simpson's rule
Int_SMP = INT.Simpson( Dens.DensFunc )
print("Integral O3: ", Int_SMP)

# Ploting the Pressure-density profile for the selected temperature
plt.plot( rho, P_eos )
plt.ylabel('Pressure')
plt.xlabel('density')
# specifying horizontal line type
plt.axhline(y = P_0, color = 'r', linestyle = '-')
# mark in liquid density
plt.plot( rho_v, P_0, marker='o', color = 'k' )
plt.plot( rho_l, P_0, marker='o', color = 'k' )
plt.show()

