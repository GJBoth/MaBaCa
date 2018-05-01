import numpy as np

from MaBaCa.Utilities.Utilities import *
from MaBaCa.General.Device_1D import *
from MaBaCa.General.Calculator_1D import *

class Magnonic_Crystal_1D(Device_1D, Calculator_1D, Utilities):
    '''
    Constructs a 1D Magnonic Crystal type of device.
    
    Note: This class doesn't do much yet, but in the future this class will be used to prepare the input for the 
    Device_1D class, i.e. Fourier transforming a magnetization profile.
    '''
    
    def __init__(self, Configuration, Number_of_k_points=20, modes=10):
        '''
        Initiates the object.
        
        Parameters
        ----------
        Configuration: Dictionary
            Dictionary containing the values of several parameters for the device. See documentation for details.
        
        Number_of_k_points: Integer
            Number of k points at which the bandstructure is calculated.
    
        Returns
        -------
    
        Raises
        ------
        
        '''
        k_points=np.linspace(0, 2*np.pi/Configuration['Periodicity'], Number_of_k_points)
        Satmag=self._Saturation_Magnetization(Configuration['Plane waves'],Configuration['Periodicity'],Configuration['Stripe width'], Configuration['Saturation magnetization stripe'],Configuration['Saturation magnetization matrix'])
        Exchange=self._Exchange(Configuration['Plane waves'],Configuration['Periodicity'],Configuration['Stripe width'],Configuration['Exchange constant stripe'],Configuration['Exchange constant matrix'], Configuration['Saturation magnetization stripe'],Configuration['Saturation magnetization matrix'])
        Device_1D.__init__(self, Configuration, Satmag, Exchange, k_points, modes)     
    
        
    def _Saturation_Magnetization(self,PlaneWaves, Periodicity, Stripe_width, SatMag_stripe, SatMag_matrix):
        Indices=np.arange(-PlaneWaves,PlaneWaves+1)
        I_Grid, J_Grid=np.meshgrid(Indices,Indices,indexing='ij')
        Vector=(I_Grid-J_Grid)*(2*np.pi/Periodicity)
        
        np.fill_diagonal(Vector,np.nan)
        Sat_Mag_factor=2/(Periodicity*Vector)*(SatMag_stripe-SatMag_matrix)*np.sin(Vector*Stripe_width/2)
        np.fill_diagonal(Sat_Mag_factor,(Stripe_width/Periodicity)*(SatMag_stripe-SatMag_matrix)+SatMag_matrix)
        return Sat_Mag_factor
    
    def _Exchange(self,PlaneWaves, Periodicity, Stripe_width, Exchange_constant_stripe, Exchange_constant_matrix, SatMag_stripe, SatMag_matrix):
        mu0=4*np.pi*10**-7
        Indices=np.arange(-PlaneWaves,PlaneWaves+1)
        I_Grid, J_Grid=np.meshgrid(Indices,Indices,indexing='ij')
        Vector=(I_Grid-J_Grid)*(2*np.pi/Periodicity)
        
        Exchange_length_stripe=(2*Exchange_constant_stripe)/(mu0*SatMag_stripe)
        Exchange_length_matrix=(2*Exchange_constant_matrix)/(mu0*SatMag_matrix)
        
        np.fill_diagonal(Vector,np.nan)
        Exchange_factor=2/(Periodicity*Vector)*(Exchange_length_stripe-Exchange_length_matrix)*np.sin(Vector*Stripe_width/2)
        np.fill_diagonal(Exchange_factor,(Stripe_width/Periodicity)*(Exchange_length_stripe-Exchange_length_matrix)+Exchange_length_matrix)
        return Exchange_factor