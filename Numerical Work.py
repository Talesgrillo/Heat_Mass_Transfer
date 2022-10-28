from cProfile import label
from cmath import exp, pi
from mimetypes import suffix_map
from turtle import color, left
import numpy as np
import matplotlib.pyplot as plt
from pandas import array

#Definition of the variables that we will use in the solution of the problem
N = 10
L = 1.0
delta_x = L/N
delta_t = 0.2*delta_x*delta_x
quadrado_x = delta_x*delta_x

#Definition of positions along the bar
x = np.linspace(0.0,L,N+1)

#Definition of the final and initial time that the simulations of the behavior in the bar will take place
tempo_final = 10.0
t = 0.0

#Here we define problem inicial conditions
Tem= np.ones(N+1,float)

#Here we define problem boundary conditions
Tem[0]= 0.0
Tem[N]= 0.0

#Here we get the temperature with the boundary conditions
temp_nova = np.copy(Tem)

#The numerical calculation formula was defined
def transcal(i):
    temp_nova = Tem [i] + (delta_t/(quadrado_x))*(Tem[i+1]-2.0*Tem[i]+Tem[i-1])
    return temp_nova

#Here I define the conditions that will start the simulations at different positions and times.
while t < tempo_final:
    for i in range (1,N):
        temp_nova[i] = transcal(i)
    Tem = np.copy (temp_nova)
    t+=delta_t

#here the exact solution was applied to confirm the values ​​obtained with the numerical method
def somatorio(x,t,inter_count):
    sum = 0
    for n in range(1,inter_count):
        fraction = 4/((2*n-1)*np.pi)
        sin = np.sin((2*n-1)*np.pi*x)
        exp = np.exp(-(2*n-1)**2*(np.pi)**2*(t))
        sum += fraction*sin*exp
    return(sum)

result = somatorio(x, t, 20)

#informations for plot de graph
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot()

#Organization of the visual structure of graphs and their plotting
fig.suptitle ('Time = %.2f (s)'%t, fontsize=18, fontweight='bold')
ax.set_ylabel( '$Temperature(°C)$' , fontsize=18)
ax.set_xlabel( '$Position(m)$' , fontsize=18)
plt.plot (x,Tem, '-r' , lw = 4,label = 'Computed Solution')
plt.plot (x,result, '-r' , lw = 4,color ='green',label = 'Exact Solution')
plt.legend(title='Legend',bbox_to_anchor=(0.65,1.14),loc='upper left',ncol=2)
plt.grid()
plt.show ()

#Exercicio 2
#Definition of the variables that we will use in the solution of the problem
N = 10
L = 1.0
delta_x = L/N
delta_t = 0.2*delta_x*delta_x
quadrado_x = delta_x*delta_x

#Definition of positions along the bar
x = np.linspace(0.0,L,N+1)

#Definition of the final and initial time that the simulations of the behavior in the bar will take place
tempo_final = 10.0
t = 0.0

#Here we define problem inicial conditions
Tem= np.zeros(N+1,float)

#Here we define problem boundary conditions
Tem[0]= 1.0

#Here we get the temperature with the boundary conditions
temp_nova = np.copy (Tem)

#The numerical calculation formula was defined
def transcal(i):
    temp_nova = Tem [i] + (delta_t/(quadrado_x))*(Tem[i+1]-2.0*Tem[i]+Tem[i-1])
    return temp_nova

#Here I define the conditions that will start the simulations at different positions and times.
while t < tempo_final:
    for i in range (1,N):
        temp_nova [i] = transcal(i)
    Tem = np.copy (temp_nova)
    t+=delta_t

#here the exact solution was applied to confirm the values ​​obtained with the numerical method
def somatorio(x,t,inter_count):
    sum = 0
    for n in range(1,inter_count):
        fraction = 2/(n*np.pi)
        sin = np.sin((n)*np.pi*x)
        exp = np.exp((-n**2)*(np.pi**2)*(t))
        sum += fraction*sin*exp
        temp_essa = 1 - x - sum
    return(temp_essa)

result = somatorio(x, t, 20)

#informations for plot de graph
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot()

#Organization of the visual structure of graphs and their plotting
fig.suptitle ('Time = %.2f (s)'%t, fontsize=18, fontweight='bold')
ax.set_ylabel( '$Temperature(°C)$' , fontsize=18)
ax.set_xlabel( '$Position(m)$' , fontsize=18)
plt.plot (x,Tem, '-r' , lw = 4,label = 'Computed Solution')
plt.plot (x,result, '-r' , lw = 4,color = 'Green',label = 'Exact Solution')
plt.legend(title='Legend',bbox_to_anchor=(0.65,1.14),loc='upper left',ncol=2)
plt.grid()
plt.show ()

# #Exercicio 3
#Definition of the variables that we will use in the solution of the problem
N = 10
L = 2.0
delta_x = L/N
delta_t = 0.2*delta_x*delta_x
quadrado_x = delta_x*delta_x

#Definition of positions along the bar
x = np.linspace(0.0,L,N+1)

#Definition of the final and initial time that the simulations of the behavior in the bar will take place
tempo_final = 10.0
t = 0.0


#Here we define problem inicial conditions
Tem=np.sin((np.pi/2)*x)

#Here we define problem boundary conditions
Tem[0]= 0.0
Tem[N]= 0.0

#Here we get the temperature with the boundary conditions
temp_nova = np.copy (Tem)

#The numerical calculation formula was defined
def transcal(i):
    temp_nova = Tem [i] + (delta_t/(quadrado_x))*(Tem[i+1]-2.0*Tem[i]+Tem[i-1])
    return temp_nova

#Here I define the conditions that will start the simulations at different positions and times.
while t < tempo_final:
    for i in range (1,N):
        temp_nova [i] = transcal(i)
    Tem = np.copy (temp_nova)
    t+=delta_t

#here the exact solution was applied to confirm the values ​​obtained with the numerical method
def somatorio(x,t,inter_count):
    for n in range(1,N):
        sin = np.sin((np.pi/2)*(x))
        expo = np.exp(((-np.pi**2)*(t))/4)
        sum = sin*expo
    return (sum)

result = somatorio(x, t, 20)

#informations for plot de graph
fig = plt.figure(figsize=(11,7))
ax = fig.add_subplot()

#Organization of the visual structure of graphs and their plotting
fig.suptitle ('Time = %.1f (s)'%t, fontsize=18, fontweight='bold')
ax.set_ylabel( '$Temperature(°C) $' , fontsize=18)
ax.set_xlabel( '$Position(m)$' , fontsize=18)
plt.plot (x,Tem, '-r' , lw = 4,label = 'Computed Solution')
plt.plot (x,result, '-r' , lw = 4,color = 'Green',label = 'Computed Solution')
plt.legend(title='Legend',bbox_to_anchor=(0.65,1.14),loc='upper left',ncol=2)
plt.grid()
plt.show ()