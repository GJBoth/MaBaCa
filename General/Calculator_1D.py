import numpy as np

class Calculator_1D():
    #this class contains the actual methods which do the calculation    
    def _Demagnetization(self,Vector):
        Demag_factor=np.exp(-np.absolute(Vector)*self._Film_thickness/2)*np.cosh(np.absolute(Vector)*self._Film_thickness/2)
        return Demag_factor
            
    def Bandstructure_calculation(self):
        Mxx=np.zeros((2*self._PlaneWaves+1,2*self._PlaneWaves+1))
        Myy=np.zeros((2*self._PlaneWaves+1,2*self._PlaneWaves+1))
        
        Indices=np.arange(-self._PlaneWaves,self._PlaneWaves+1)
        I_Grid, J_Grid=np.meshgrid(Indices,Indices,indexing='ij')
        I_Grid=I_Grid*self._Reciprocal_lattice_vector
        J_Grid=J_Grid*self._Reciprocal_lattice_vector
        
        SatMag=self._Saturation_magnetization_matrix/(self._External_field/self._Vacuum_permeability)
        Exchange=self._Exchange_matrix/(self._External_field/self._Vacuum_permeability)
        
        for k_index, k in enumerate(self._k_points):
            Mxy=np.identity(2*self._PlaneWaves+1)+(1-self._Demagnetization(k+J_Grid))*SatMag+Exchange*(k+J_Grid)**2
            Myx=-np.identity(2*self._PlaneWaves+1)-self._Demagnetization(k+J_Grid)*SatMag-Exchange*(k+J_Grid)**2
        
            M=np.concatenate((np.concatenate((Mxx,Myx),axis=0),np.concatenate((Mxy,Myy),axis=0)),axis=1)
        
            eigenvalues, eigenvectors=np.linalg.eig(M)
            
            eigenvalues, unique_list=np.unique(np.absolute(np.imag(eigenvalues))*self._Gyromagnetic_ratio*self._External_field/(2*np.pi),return_index=True)
            eigenvectors=eigenvectors[:,unique_list]
            
            k_values=k*np.ones(len(eigenvalues))
            eigenvectors=eigenvectors.transpose()
            eigenvectors=eigenvectors.tolist()
            
            self._bandstructure[k_index*self._modes:(k_index+1)*self._modes,0]=k_values[:self._modes]
            self._bandstructure[k_index*self._modes:(k_index+1)*self._modes,1]=np.arange(len(eigenvalues[:self._modes]))
            self._bandstructure[k_index*self._modes:(k_index+1)*self._modes,2]=eigenvalues[:self._modes]
            self._bandstructure[k_index*self._modes:(k_index+1)*self._modes,3]=eigenvectors[:self._modes]
