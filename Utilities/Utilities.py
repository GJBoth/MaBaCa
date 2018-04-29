import numpy as np

from plotly.offline import *
from plotly.graph_objs import *

class Utilities():
    def Bandstructure_plot(self,number_of_modes=3,width=4,height=8):
        #this function plots the bandstructure; we plot all lines but zoom in on the first three (standard)
        trace=[]
        for modenumber in np.arange(self._modes):
            mode=self.return_mode(modenumber)
            x=mode[:,0]*self._Periodicity/(2*np.pi)
            y=mode[:,2]
            trace.append(Scatter(x=x,y=y,name='Mode '+str(modenumber)))
            if modenumber==0:
                y_plot_min=min(y)
            if modenumber==number_of_modes-1:
                y_plot_max=max(y)
        
        #Setting plot range
        plotrange=y_plot_max-y_plot_min
        y_plot_min=y_plot_min-plotrange/10 
        if y_plot_min<0:
            y_plot_min=0
        y_plot_max=y_plot_max+plotrange/10
        
        layout=Layout(xaxis=dict(title='k (2pi/a)',range=[0,1],fixedrange=True),
                      yaxis=dict(title='Frequency (GHz)',range=[y_plot_min,y_plot_max],rangemode='nonnegative'),
                      showlegend=False,
                      width=width,
                      height=height,
                      hovermode = 'closest')
        figure=Figure(data=trace,layout=layout)
        config = {'scrollZoom': True}

        plot(figure,config=config)
        
        
    def Spinwave_profile(self,plottype,modenumber, k_point,xlimits):
        #this function plots the spinwave profile of modenumber between xlimits at k
        
        self._profile=np.zeros(int(self._PlaneWaves/2),6)
        #finding eigenvector with closest k-value
        mode=self.return_mode(modenumber)
        eigenvector=mode[np.argmin(np.absolute(mode[:,0]*self._Periodicity/(2*np.pi)-k_point)),3]
        
        eigenvector_x=np.array(eigenvector[0:(2*self._PlaneWaves+1)])
        eigenvector_y=np.array(eigenvector[(2*self._PlaneWaves+1):])
        
        x=np.linspace(xlimits[0]*self._Periodicity,xlimits[1]*self._Periodicity,int(self._PlaneWaves/2))
        G=np.arange(-self._PlaneWaves,self._PlaneWaves+1,1)*self._Reciprocal_lattice_vector
        G,x=np.meshgrid(G,x)
        
        factor=np.exp(1j*(k_point+G)*x)
        
        Sw_x=np.matmul(factor,eigenvector_x)
        Sw_y=np.matmul(factor,eigenvector_y)
        
        self._profile[:,0]=np.absolute(Sw_x)
        self._profile[:,1]=np.angle(Sw_x)
        self._profile[:,2]=np.real(np.absolute(Sw_x)*np.exp(-1j*np.angle(Sw_x)))
        self._profile[:,3]=np.absolute(Sw_y)
        self._profile[:,4]=np.angle(Sw_y)
        self._profile[:,5]=np.real(np.absolute(Sw_y)*np.exp(-1j*np.angle(Sw_y)))
        xplot=np.linspace(xlimits[0],xlimits[1],int(self._PlaneWaves/2))
        if plottype=='Phase':
            trace1=Scatter(x=xplot,y=self._profile[:,1]/np.pi,name='Out-of-plane')
            trace2=Scatter(x=xplot,y=self._profile[:,4]/np.pi,name='In-plane')
            layout=Layout(xaxis=dict(title='y (Periodicity)'),yaxis=dict(title='Phase (pi)'))
                     
        elif plottype=='Amplitude':
            normalization_factor=np.amax(np.absolute(self._profile[:,[0,3]]))
            trace1=Scatter(x=xplot,y=self._profile[:,0]/normalization_factor,name='Out-of-plane')
            trace2=Scatter(x=xplot,y=self._profile[:,3]/normalization_factor,name='In-plane')
            layout=Layout(xaxis=dict(title='y (Periodicity)'),yaxis=dict(title='Amplitude (norm., arb.)'))
                     
        elif plottype=='Total':
            normalization_factor=np.amax(np.absolute(self._profile[:,[2,5]]))
            trace1=Scatter(x=xplot,y=self._profile[:,2]/normalization_factor,name='Out-of-plane')
            trace2=Scatter(x=xplot,y=self._profile[:,5]/normalization_factor,name='In-plane')
            layout=Layout(xaxis=dict(title='y (Periodicity)'),yaxis=dict(title='M (norm., arb.)'))
        
        elif plottype=='M2':
            normalization_factor=np.amax(np.absolute(self._profile[:,[2,5]]))
            trace1=Scatter(x=xplot,y=(self._profile[:,2]/normalization_factor)**2,name='Out-of-plane')
            trace2=Scatter(x=xplot,y=(self._profile[:,5]/normalization_factor)**2,name='In-plane')
            layout=Layout(xaxis=dict(title='y (Periodicity)'),yaxis=dict(title='M^2 (norm., arb.)'))
                     
        data=[trace1,trace2]
        
        figure=Figure(data=data,layout=layout)
        plot(figure)
    
        