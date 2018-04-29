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
    
    def __init__(self, Configuration, Number_of_k_points):
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
        Device_1D.__init__(self, Configuration, k_points)     