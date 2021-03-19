# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 08:49:05 2021

@author: Dofi
"""
import numpy as np
import matplotlib.pyplot as plt
#posisi source
sx=[0,0,1000]
sy=[-250,-633,0]
rx=[2000]
ry=[-1000]
tobs=[3.018,2.894,1.141]
v0=883.27
v1=v0
v2=v0

k=0
rms=[]
while k<5:
    if k==0:
        v1=v0
        v2=v0
    else:
        v1=v1[0]
        v2=v2[0]

#ray tracing menggunakan hukum snell
#teta=np.arange(0,np.pi/2,np.pi/360)
    l1=[]
    l2=[]
    tcall=[]
    deltat=[]
    y11=[]
    y12=[]
    y21=[]
    y22=[]
    for i in range(len(sx)):
        if i==2:
            l1.append(0)
            l2.append(1414.2)
            tcall.append(l1[i]*1/v1+l2[i]*1/v2)
            deltat.append(tobs[i]-tcall[i])
        else:
            y=7
            teta1=0
            while y>0.05:
                teta1=teta1+np.pi/3600000
                teta2=np.arctan(np.sin(teta1)*v2/v1)
                dy1=np.tan(teta1)*1000
                a=sy[i]-dy1
                dy2=np.tan(teta2)*1000
                aa=a-dy2
                y=abs(1000+aa)
                if i==0 and k==4:
                    y11.append(a)
                    y12.append(aa)
                elif i==1 and k==4:
                    y21.append(a)
                    y22.append(aa)
            l1.append(1000/np.cos(teta1))
            l2.append(1000/np.cos(teta2))
            tcall.append(l1[i]*1/v1+l2[i]*1/v2)
            deltat.append(tobs[i]-tcall[i])
            
#inversi matrik
    A=np.array(l1+l2).reshape(2,3)
    At=np.matrix.transpose(A)
    delta=np.array(deltat).reshape(3,1)
    rms.append(sum(delta**2)**0.5)
    Ax=np.dot(A,At)
    Axinv=np.linalg.inv(Ax)
    S=np.dot(Axinv,A)
    S=np.dot(S,delta)
    dv1=-S[0]*v1**2/(1+S[0]*v1)
    dv2=-S[1]*v2**2/(1+S[1]*v2)
    v1=v1+dv1
    v2=v2+dv2
    v1=v1.tolist()
    v2=v2.tolist()
    k=k+1
    
y1=y21[1:len(y21):1000]
y1=y1[-20:]
y2=y22[1:len(y22):1000]
y2=y2[-20:]
y3=y11[1:len(y11):2000]
y3=y3[-20:]
y4=y12[1:len(y12):2000]
y4=y4[-20:]
#plot ray tracy using snell's law
for i,j,k,l in zip(y1,y2,y3,y4):
    plt.plot([0,1000,2000],[-633,i,j],'-r',[0,1000,2000],[-250,k,l],'-r')
    plt.plot([0,1000,2000],[-633,y1[-1],y2[-1]],'-g',[0,1000,2000],[-250,y3[-1],y4[-1]],'-g')
    plt.plot([1000,2000],[0,-1000],'-g')
    plt.plot([1000,1000],[0,-1000],'-b')
    plt.annotate('*',xy=(0,-300),ha='center',fontsize=25,color='yellow')
    plt.annotate('*',xy=(0,-653),ha='center',fontsize=25,color='yellow')
    plt.annotate('*',xy=(1000,-5),ha='center',fontsize=25,color='yellow')
    plt.annotate('V',xy=(2000,-1000),ha='center',fontsize=25,color='green')
    plt.ylim(-1000,0)
    plt.xlim(0,2000)
    
plt.show()

#plot rms
plt.plot(rms)
plt.xlabel('iter')
plt.ylabel('rms')
plt.show()
#Plot velocity and ray tracing
xgrid=np.arange(0,2100,100)
ygrid=np.arange(0,-1100,-100)
X,Y=np.meshgrid(xgrid,ygrid)
z=np.zeros((11,21))
for i in range(11):
    for j in range(21):
        if j<10:
            z[i,j]=z[i,j]+v1
        else:
            z[i,j]=z[i,j]+v2
plt.pcolor(X,Y,z)
plt.plot([0,1000,2000],[-633,y1[-1],y2[-1]],'-g',[0,1000,2000],[-250,y3[-1],y4[-1]],'-k')
plt.plot([1000,2000],[0,-1000],'-b')
plt.plot([1000,1000],[0,-1000],'-r')
plt.colorbar()

