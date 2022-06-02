print("Code Ready")
import numpy as np
import matplotlib.pylab as plt
import scipy.stats
import pandas as pd
print("Imports ready")

########################
#Simulation conditions #
########################
time_unit=0.0001
steps=15000
nbodies=2
#In case of not using X(Omega) set xbase
xbase=0.15
#Steps from xbase or X(Omega)
xi=0.002

#For figures
lim=4
Save=False

###################### 
# Initial conditions #
######################
# masita is a placeholder for the mass as this script can be used to simulate other systems
Energy=-1.0
pyi=1.0
masita=1
epsilon=0.0005
#To use golden ration set winding=phi
phi=(np.sqrt(5.0)-1)/2
#Insert winding number of choice
winding=1.0/2.0
yi=0 


################################
# Hamilton equations of motion #
################################
def PXpunto(x,y,px,py,t,m):
    p1=np.sqrt(((x+epsilon)**2)+(y**2))
    p2=np.sqrt(((x+epsilon-1)**2)+(y**2))
    Ppuntox=-(((1-epsilon)*(x+epsilon))/(p1**3))+py-(((epsilon)*(x+epsilon-1))/(p2**3))
    return Ppuntox

def Xpunto(x,y,px,py,t,m):
    return px+y

def PYpunto(x,y,px,py,t,m):
    p1=np.sqrt(((x+epsilon)**2)+(y**2))
    p2=np.sqrt(((x+epsilon-1)**2)+(y**2))
    Ppuntoy=-(((1-epsilon)*(y))/(p1**3))-px-(((epsilon)*(y))/(p2**3))
    return  Ppuntoy

def Ypunto(x,y,px,py,t,m):
    return py-x
################################
   
def px_inicial(E,x,y,py):
    """Calculates the initial X momentum given the other initial conditions

    Args:
        E (Scalar): Energy
        x (Scalar): Initial position that may be calculated from Windingpos
        y (Scalar): Initial position on y axis, Normally set to 0
        py (Scalar): Initial momentum in y axis recommended to use 1 in initial conditions

    Returns:
        scalar: Initial momentum in x axis
    """    
    p1=np.sqrt(((x+epsilon)**2)+(y**2))
    p2=np.sqrt(((x+epsilon-1)**2)+(y**2))
    px=np.sqrt(2*(E+((1-epsilon)/(p1))+((epsilon)/(p2))-(((py-x)**2)/2)+((x**2+y**2)/2)))-y
    print(fr"Px inicial={px}")
    return px


def Windingpos(w):
    """Calculates the initial position in terms of R

    Args:
        w (Fraction or Scalar): Winding number set at the initial conditions phase

    Returns:
        Position: Initial position corresponding to w rounded to 6 decimal places
    """    
    Pos=(-0.5*((1/(1-w))**(2/3)))-Energy
    return round(Pos,6)

def simulate(cuerpo,New_time_units=time_unit,steps=steps):
    """Simulates a body with Hamilton equations of motion and initial conditions set by the user

    Args:
        cuerpo (Body): A body with initial conditons set
        New_time_units (Scalar, optional): In case other time units must be used. Defaults to time_unit.
        steps (Int, optional): Number of steps to simulate. Defaults to steps.
    """    
    for i in range(steps):
        time=i*New_time_units
        cuerpo.updateMom(time,time_units=New_time_units)
        cuerpo.updatePos(time,time_units=New_time_units)
    return(cuerpo.historyPos,cuerpo.historyMom)

class Body:
    def __init__(self,x0,y0,px0,py0,mass0,time_units=time_unit):
        self.x=x0
        self.y=y0
        self.px=px0
        self.py=py0
        self.historyPos=[(x0,y0)]
        self.historyMom=[(px0,py0)]
        self.mass=mass0
        self.time_units=time_units
    
    def updatePos(self,t,time_units=time_unit):
        self.x=self.x+(Xpunto(self.x,self.y,self.px,self.py,t,self.mass)*time_units)
        self.y=self.y+(Ypunto(self.x,self.y,self.px,self.py,t,self.mass)*time_units)
        self.historyPos.append((self.x,self.y))
    
    def updateMom(self,t,time_units=time_unit):
        self.px=self.px+(PXpunto(self.x,self.y,self.px,self.py,t,self.mass)*time_units)
        self.py=self.py+(PYpunto(self.x,self.y,self.px,self.py,t,self.mass)*time_units)
        self.historyMom.append((self.px,self.py))



print("Definitions Ready")      

x_inicial=np.array([Windingpos(winding)+(i*xi) for i in range(nbodies)])
WindingNumbers=[((1+np.sqrt(5))/2)-1]
WindingNumbersStr=[winding for i in range(nbodies)]
Energias=[Energy for i in range(nbodies)]
EnergiasStr=[str(Energias[i]) for i in range(nbodies)]

# exit()
planets={}
Posiciones={}
Posiciones_x={}
Posiciones_y={}

Momentos={}
Momento_y={}
Momento_x={}

DataFrames={}
CruzY_0={}
Ypunto_posi={}
for aBody in range(nbodies):
    planets[aBody]=Body(x_inicial[aBody],0,px_inicial(Energias[aBody],x_inicial[aBody],0,pyi),pyi,masita,time_units=time_unit)
    # planets[aBody]=Body(x_inicial[aBody],yi,pxi,py_inicial(Energias[aBody],x_inicial[aBody],yi,pxi),masita,time_units=time_unit)
    Posiciones[aBody],Momentos[aBody]=simulate(planets[aBody],New_time_units=time_unit,steps=steps)
    Posiciones_x[aBody]=[tuplee[0] for tuplee in Posiciones[aBody]]
    Posiciones_y[aBody]=[tuplee[1] for tuplee in Posiciones[aBody]]
    CruzY_0[aBody]=np.abs(np.diff(0.5*np.sign(Posiciones_y[aBody]+[0])))
    
    Momento_x[aBody]=[tuplee[0] for tuplee in Momentos[aBody]]
    Momento_y[aBody]=[tuplee[1] for tuplee in Momentos[aBody]]
    Ypunto_posi[aBody]=(np.sign(Momento_y[aBody])+1)
    
    DataFrames[aBody]=pd.DataFrame({f"PosX{aBody}":Posiciones_x[aBody],f"PosY{aBody}":Posiciones_y[aBody],f"MomX{aBody}":Momento_x[aBody],f"MomY{aBody}":Momento_y[aBody],f"CruzY{aBody}":CruzY_0[aBody],f"Ypunto_posi{aBody}":Ypunto_posi[aBody]})
      
    print(f"Body number {aBody+1} simulation done")
    


def PlotTrayectory(save=False):
    """Plots the trajectory for the nbodies, if initial conditions are not correct the system quickly diverges

    Args:
        save (bool, optional): User chooses if figure is shown or saved to pc. Defaults to False.
    """    
    fig=plt.figure()
    # plt.title(fr"Trayectories of restricted" "\n" fr"body problem with $\epsilon$={epsilon}")
    for aBody in range(nbodies):
        plt.plot(Posiciones_x[aBody],Posiciones_y[aBody],label=r"$x_0=$"+f"{x_inicial[aBody]} R" "\n" r"$\Omega=$"+f"{winding}",linestyle=(0, (5, 1)),linewidth=0.5)
    plt.xlabel("X (R)")
    plt.ylabel("Y (R)")
    plt.scatter([-epsilon],[0],s=50, label="Main Attractor",c="k")
    if epsilon>0.0:
        plt.scatter([1-epsilon],[0],s=30,label="Perturbator",c="g")
        plt.xlim((-lim,lim))
        plt.ylim((-lim,lim))
    plt.legend()
    if save:
        #Insert saving location
        plt.savefig(f"/home/bleon/Documents/Maestria/Clases/Mecanica_Analitica/Caos_proyecto/trayec_e{epsilon}_W{WindingNumbersStr[0]}.pdf")
        print(f"Figure has been saved with name:trayec_e{epsilon}_W{WindingNumbersStr[0]}!")
    else:
        plt.show()



def PlotPoincare():
    """Plots poincare section for all bodies used
    """    
    fig=plt.figure()
    for aBody in range(nbodies):
        df=DataFrames[aBody][(DataFrames[aBody][f"Ypunto_posi{aBody}"]>0) & (DataFrames[aBody][f"CruzY{aBody}"]>0)]
        plt.scatter(df[f"PosX{aBody}"],df[f"MomX{aBody}"],s=2,label=fr"body {aBody+1}, $\Omega=$ {WindingNumbersStr[aBody]}, E={EnergiasStr[aBody]} ")
    plt.legend()
    plt.show()

#Uncomment what you want to plot
# PlotTrayectory(save=Save)
# PlotPoincare()

Lya=pd.DataFrame({"x0":Posiciones_x[0],"x1":Posiciones_x[1],"y0":Posiciones_y[0],"y1":Posiciones_y[1]})
Lya["distY"]=Lya["y1"]-Lya["y0"]
Lya["distX"]=Lya["x1"]-Lya["x0"]
Lya["Distance"]=np.sqrt((Lya["distY"]**2)+(Lya["distX"]**2))

#Change the 4000 in case you want to calculate exponent from beggining of time
N=range(steps-4000)
dist=Lya["Distance"].to_numpy()[4000:]

#Calculate Lyapunov exponents for first 2 bodies
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(N, np.log(dist)[:-1])
print("The luyapunov coefficient is:")
print(slope)
print("The correlation for this simulation is")
print(r_value)

#Uncomment if you want to plot 
# fig=plt.figure()
# plt.plot(N,dist[:-1])
# plt.show()
