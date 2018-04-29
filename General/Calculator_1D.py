import numpy as np

class Calculator_1D():
    #this class contains the actual methods which do the calculation    
    def _Demagnetization(self,Vector):
        Demag_factor=np.exp(-np.absolute(Vector)*self._Film_thickness/2)*np.cosh(np.absolute(Vector)*self._Film_thickness/2)
        return Demag_factor
    
    def _Saturation_Magnetization(self,Vector):
        np.fill_diagonal(Vector,np.nan)
        Sat_Mag_factor=2/(self._Periodicity*Vector)*(self._SatMag_stripe-self._SatMag_matrix)*np.sin(Vector*self._Stripe_width/2)
        np.fill_diagonal(Sat_Mag_factor,(self._Stripe_width/self._Periodicity)*(self._SatMag_stripe-self._SatMag_matrix)+self._SatMag_matrix)
        return Sat_Mag_factor
    
    def _Exchange(self,Vector):
        np.fill_diagonal(Vector,np.nan)
        Exchange_factor=2/(self._Periodicity*Vector)*(self._Exchange_stripe-self._Exchange_matrix)*np.sin(Vector*self._Stripe_width/2)
        np.fill_diagonal(Exchange_factor,(self._Stripe_width/self._Periodicity)*(self._Exchange_stripe-self._Exchange_matrix)+self._Exchange_matrix)
        return Exchange_factor
            
    def Bandstructure_calculation(self):
        Mxx=np.zeros((2*self._PlaneWaves+1,2*self._PlaneWaves+1))
        Myy=np.zeros((2*self._PlaneWaves+1,2*self._PlaneWaves+1))
        
        Indices=np.arange(-self._PlaneWaves,self._PlaneWaves+1)
        I_Grid, J_Grid=np.meshgrid(Indices,Indices,indexing='ij')
        I_Grid=I_Grid*self._Reciprocal_lattice_vector
        J_Grid=J_Grid*self._Reciprocal_lattice_vector
        
        data_k=[]
        data_eigenvalues=[]
        data_eigenvectors=[]
        data_mode=[]
        
        SatMag_matrix=self._Saturation_Magnetization(I_Grid-J_Grid)/(self._External_field/self._Vacuum_permeability)
        Exchange_matrix=self._Exchange(I_Grid-J_Grid)/(self._External_field/self._Vacuum_permeability)
        
        for k in self._k_points:
            Mxy=np.identity(2*self._PlaneWaves+1)+(1-self._Demagnetization(k+J_Grid))*SatMag_matrix+Exchange_matrix*(k+J_Grid)**2
            Myx=-np.identity(2*self._PlaneWaves+1)-self._Demagnetization(k+J_Grid)*SatMag_matrix-Exchange_matrix*(k+J_Grid)**2
        
            M=np.concatenate((np.concatenate((Mxx,Myx),axis=0),np.concatenate((Mxy,Myy),axis=0)),axis=1)
        
            eigenvalues, eigenvectors=np.linalg.eig(M)
            
            eigenvalues, unique_list=np.unique(np.absolute(np.imag(eigenvalues))*self._Gyromagnetic_ratio*self._External_field/(2*np.pi),return_index=True)
            eigenvectors=eigenvectors[:,unique_list]
            
            k_values=k*np.ones(len(eigenvalues))
            eigenvectors=eigenvectors.transpose()
            eigenvectors=eigenvectors.tolist()
            
            data_k.extend(k_values[:self._modes])
            data_eigenvalues.extend(eigenvalues[:self._modes])
            data_eigenvectors.extend(eigenvectors[:self._modes])
            data_mode.extend(np.arange(len(eigenvalues[:self._modes])))
            
        self._bandstructure[:,0]=data_k
        self._bandstructure[:,1]=data_mode
        self._bandstructure[:,2]=data_eigenvalues
        self._bandstructure[:,3]=data_eigenvectors