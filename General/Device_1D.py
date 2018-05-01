import numpy as np
import pickle


class Device_1D():
    #this class represents the device
    def __init__(self, parameters, SatMag, Exchange, k_points, modes):
        #this function initializes the device
        self._External_field=parameters['External Field']
        
        self._Periodicity=parameters['Periodicity']
        self._Stripe_width=parameters['Stripe width']
        self._Film_thickness=parameters['Film thickness']
        self._PlaneWaves=parameters['Plane waves']
        
        self._Gyromagnetic_ratio=parameters['Gyromagnetic ratio']
        self._Vacuum_permeability=4*np.pi*10**-7
        
        self._Reciprocal_lattice_vector=2*np.pi/self._Periodicity
        
        self._k_points=k_points
        self._modes=modes
        
        self._bandstructure=np.zeros((self._modes*len(self._k_points),4),dtype=object)
        
        self._Saturation_magnetization_matrix=SatMag
        self._Exchange_matrix=Exchange

    def _save(self,filename):
        #pickle output, can be very big!
        with open(filename+'.pickle', 'wb') as output_file:
            pickle.dump(self, output_file)
    
    def export(self,filename):
        #save bandstructure and parameters in separate files
        
        #bandstructure
        np.savetxt(filename+'.csv', self._bandstructure[:,:-1], delimiter=",")
        
        #parameters
        skiplist=['_Vacuum_permeability','_k_points','_bandstructure','_profile','_Saturation_magnetization_matrix','_Exchange_matrix']
        parameterdict={}
        for item in vars(self):
            if item not in skiplist:
                value=eval('self.'+str(item))
                item=item.replace('_',' ')
                item=item.strip()
                parameterdict[str(item)]=value
                
        devicetype=type(self).__name__
        devicetype=devicetype.replace('_',' ')
        with open(filename+'_parameters.txt', 'w') as file:
            file.write('Device type: '+devicetype+'\n')
            [file.write('{0}: {1}\n'.format(key, value)) for key, value in parameterdict.items()]

        
    def return_mode(self,modenumber):
        return self._bandstructure[self._bandstructure[:,1]==modenumber,:]
        
    def set_field(self,field):
        #update the external field
        self._External_field=field
        
    def set_planewaves(self,planewaves):
        self._PlaneWaves=planewaves
        
        