#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 09:13:00 2020

@author: christopherwilson
"""
import numpy as np

import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.linalg
from scipy import integrate
from scipy import special 
from scipy import misc
from scipy import special 
from mpmath import *
from scipy import interpolate
from numpy import polynomial 


"Methods involving basic void catalog analysis"

def buildDensityProfile(cat,voidType=False,rMin=5,rMax=15,deltaMin=-1,deltaMax=1000,weight=0,outTo=50):
    densVals=[[] for i in range(outTo)]
    xpoints=cat[0][3][0][0:outTo]

    if isinstance(voidType, list):
        voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] in voidType or voidType == "all"))]
    else:
        voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] >=deltaMin and v[1] <= deltaMax))]

    for void in voidsToUse:
        for i in range(outTo):
            densVals[i]=np.append(densVals[i],void[3][1][i]*(void[2]**weight))

    dataToSave=np.array([xpoints,[np.average(densVals[i]) for i in range(outTo)],[np.std(densVals[i])/np.sqrt(len(densVals[i])) for i in range(outTo)]])
    return dataToSave



def buildIntegratedDensityProfile(cat,voidType=["R"],rMin=5,rMax=15,catType='halo',weight=0,outTo=30):
    densityVals=[[] for i in range(30)]
    if (catType=="part"):
        nBar=0.05 
        
    xpoints=np.array(cat[0][4][0])+.05
    voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] in voidType or voidType == "all"))]
    
    for void in voidsToUse:
        if (catType=="halo"):
            nBar=void[0][2][0]
        for i in range(outTo):
            numPartTot=0
            density=0
            for j in range(i+1):
                delta=void[3][1][j]
                inRad=0
                
                "BIN SEPARATION IS HARD CODED IN!!!"
                inRad=(void[3][0][j]-0.05)*void[2]
                outRad=(void[3][0][j]+0.05)*void[2]
                shellVol=4/3*np.pi*(outRad**3-inRad**3)
                numParts=round(shellVol*(nBar*delta+nBar))
                totalVol=4/3*np.pi*(outRad**3)
                numPartTot+=numParts
                density=((numPartTot/totalVol)-nBar)/nBar
        
            densityVals[i]=np.append(densityVals[i],density*(void[2]**weight))

    dataToSave=np.array([xpoints,[np.average(densityVals[i]) for i in range(outTo)],[np.std(densityVals[i])/np.sqrt(len(densityVals[i])) for i in range(outTo)]])
    return dataToSave



def buildVelocityProfile(cat,voidType=["R"],rMin=10,rMax=20,weight=0,catType='part',corFactor=1,minFromCenter=0,outTo=30):
    
    velocityVals=[[] for i in range(outTo)]
    xpoints=cat[0][4][0]
    voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] in voidType or voidType == "all"))]
   
    for void in voidsToUse:
        for i in range(outTo):
            if (void[4][1][i]!=0 and void[3][0][i]*void[2] >= minFromCenter):
                velocityVals[i]=np.append(velocityVals[i],void[4][1][i]*(void[2]**weight))

    dataToSave=np.array([xpoints,[np.average(velocityVals[i])*corFactor for i in range(outTo)],[np.std(velocityVals[i])/np.sqrt(len(velocityVals[i]))*corFactor for i in range(outTo)]])
    
    return dataToSave



def buildIntegratedVelocityProfile(cat,voidType=["R"],rMin=5,rMax=15,catType='part',weight=0,corFactor=1,startIndex=0,minFromCenter=0,outTo=30):

    velocityVals=[[] for i in range(outTo)]
    if (catType=="part"):
        nBar=0.05 

    xpoints=np.array(cat[0][4][0])+.05
    voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] in voidType or voidType == "all"))]

    for void in voidsToUse:
        if (catType=="halo"):
            nBar=void[0][2][0]
        numPartTot=0
        vel=0
        for i in range(outTo):
            if void[3][0][i]*void[2] >= minFromCenter:
                
                delta=void[3][1][i]
                inRad=0
                
                "BIN SEPARATION IS HARD CODED IN!!!"
                inRad=(void[3][0][i]-0.05)*void[2]
                outRad=(void[3][0][i]+0.05)*void[2]
                shellVol=4/3*np.pi*(outRad**3-inRad**3)
                numParts=round(shellVol*(nBar*delta+nBar))
                numPartTot+=numParts
                vel+=numParts*void[4][1][i]
            
                if (numPartTot!=0):
                    velocityVals[i]=np.append(velocityVals[i],vel/numPartTot*(void[2]**weight))

    dataToSave=np.array([xpoints,[np.average(velocityVals[i])*corFactor for i in range(outTo)],[np.std(velocityVals[i])/np.sqrt(len(velocityVals[i]))*corFactor for i in range(outTo)]])
    return dataToSave



def reclassifyVoidsRS(cat,Rcut=1.0,densCut=0,fromParts=False):
    
    for void in cat:
        "here we classify with halos or particles"
        if fromParts==True:
            nBar=0.05
        else:
            nBar=void[0][2][0]
        xpoints=void[3][0]
        Reff=void[2]
        i=0
        numParts=0
        "Check S"
        while (xpoints[i]+.0499<=Rcut):
            inRad=(xpoints[i]-0.05)*Reff
            outRad=(xpoints[i]+0.05)*Reff
            delta=void[3][1][i]
            shellVol=4/3*np.pi*(outRad**3-inRad**3)
            numParts+=shellVol*(nBar*delta+nBar)
            i+=1
        
        if (((numParts/(4/3*np.pi*(Rcut*Reff)**3)-nBar)/nBar)>densCut):
            void[1]="S"
        else:
            void[1]="R"



def matchClassificatons(cat1,cat2): 
    "the first cat has already been resorted, the second one will have its classlifications matched to the first one"
    "To have an apple to apples comparison, this should be run anytime after reclasifyVoidsRS is run"            
    if (len(cat1) != len(cat2)):
        print("These catalogs are not the same, should check which catalogs you are resorting")
        return 0
    for i in range(len(cat1)):
        cat2[i][1]=cat1[i][1]
    return cat2


def getAvgReff(cat,voidType=["R"],rMin=5,rMax=15):     
    voidsToUse=[v for v in cat if (v[2] > rMin and v[2] < rMax and (v[1] in voidType or voidType == "all"))]
    avgReff=0
    for void in voidsToUse:
        avgReff+=void[2]
    avgReff=avgReff/len(voidsToUse)
    print("average R_eff = " + str(avgReff))
    return avgReff

"Methods involving calculating forces in voids"

"these methods get called incredibly frequently during other methods, want them to be as bare as possible. These are all greens functions"
def Fminus(r,R=1.2,R_eff=15,a=1,mu=0.41908449189639996,rho=77985400000.0,G=4.301*10**-9):
    return 4*np.pi*G*rho*R/(3*mu)*np.exp(-a*R*R_eff*mu)*((a*R_eff*mu)*(np.cosh(a*r*R_eff*mu))/r-(np.sinh(a*r*R_eff*mu))/(r**2))

def Fplus(r,R=1.2,R_eff=15,a=1,mu=0.41908449189639996,rho=77985400000.0,G=4.301*10**-9):
    return 4*np.pi*G*rho*R/(3*mu)*np.sinh(a*R*R_eff*mu)*(-(a*R_eff*mu)*(np.exp(-a*r*R_eff*mu))/r-(np.exp(-a*r*R_eff*mu))/(r**2))

def Gminus(r,R=1.2,R_eff=15,a=1,mu=0.41908449189639996,rho=77985400000.0,G=4.301*10**-9):
    return (8*np.pi*a*G*rho/(3*mu))*R/r*R_eff*np.exp(-a*R*R_eff*mu)*np.sinh(a*r*R_eff*mu)

def Gplus(r,R=1.2,R_eff=15,a=1,mu=0.41908449189639996,rho=77985400000.0,G=4.301*10**-9):
    return (8*np.pi*a*G*rho/(3*mu))*R/r*R_eff*np.exp(-a*r*R_eff*mu)*np.sinh(a*R*R_eff*mu)

def FplusNewton(r,R=1.2,R_eff=15,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -4*a*G*np.pi*R_eff*rho*R**2/(r**2)

def GplusNewton(r,R=1.2,R_eff=15,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -(4*G*np.pi*rho*(R*R_eff*a)**2)/r

def GminusNewton(r,R=1.2,R_eff=15,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -4*((a*R_eff)**2)*G*np.pi*R*rho

"Should explain these methods"
def FplusNewtonR1(r,R=1.2,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -4*a*G*pi*rho*R**2/(r**2)

def GplusNewtonR2(r,R=1.2,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -(4*G*pi*rho*(R*a)**2)/r

def GminusNewtonR2(r,R=1.2,a=1,rho=77985400000.0,G=4.301*10**-9):
    return -4*((a)**2)*G*pi*R*rho


"Below methods calculates the linear fifth force and the newtonian force using greens functions"


def calcGreensFifthForce(densProf,R_eff,F0=1e-6,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    a=1/(1+z)
    mu=(1/2997.92458*(1/(2*F0))**(1/2))*((omega_M*a**(-3)+4*omega_L)**(3/2))/(omega_M+4*omega_L)
    "mu has units of h/Mpc"
    
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    "using H0=100h"
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    
    forcePoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        "Out integrates over shells exterior to the location r"
        In=sp.integrate.quad(lambda R: (max(deltaInterp(R)+0,-1))*Fplus(r,float(R),R_eff,a,mu,rho,G),dx,r-dx,limit=Limit)[0]
        Out=sp.integrate.quad(lambda R: (max(deltaInterp(R)+0,-1))*Fminus(r,float(R),R_eff,a,mu,rho,G),r+dx,outTo,limit=Limit)[0]
        forcePoints[i]=In+Out
    
    convFac=29979245.8 
    "dividing by convFac changes units from h/Mpc (km/s)^2 to H0*c"
    
    "Force is returned with units of H0*c"
    return [xpoints,forcePoints/convFac]


def calcGreensFifthField(densProf,R_eff,F0=1e-6,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    "field that is output has units of (km/s)^2, to convert back to a dimensionless field, divide by (2.997e5)**2"
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    a=1/(1+z)
    mu=(1/2997.92458*(1/(2*F0))**(1/2))*((omega_M*a**(-3)+4*omega_L)**(3/2))/(omega_M+4*omega_L)
    "mu has units of h/Mpc"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    "using H0=100h"
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    fieldPoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        "Out integrates over shells exterior to the location r"
        In=sp.integrate.quad(lambda R: (max(deltaInterp(R)+0,-1))*Gplus(r,float(R),R_eff,a,mu,rho),dx,r-dx,limit=Limit)[0]
        Out=sp.integrate.quad(lambda R: (max(deltaInterp(R)+0,-1))*Gminus(r,float(R),R_eff,a,mu,rho),r+dx,5,limit=Limit)[0]
        fieldPoints[i]=In+Out
    fieldPoints[0]=fieldPoints[1]
    "fieldPoints has units of (km/s)^2, to convert back to a dimensionless field, divide by (2.997e5)**2"
    "If this is changed, things will also have to be changed within the convergeForce method, which depends on the output of this calculation"
    
    "fieldPoints is the field's deviation from its background value, which is what gives the fifth force"
    
    return [xpoints,fieldPoints]



def calcGreensNewtonForce(densProf,R_eff,F0=1e-6,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    a=1/(1+z)
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    #xpoints=np.linspace(0,R_eff_final,400)
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    
    forcePoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*FplusNewton(r,float(R),R_eff,a,rho),dx,r,limit=Limit)[0]
        forcePoints[i]=In
    
    convFac=29979245.8 
    "dividing by convFac changes units from h/Mpc (km/s)^2 to H0*c"
    "Force is returned with units of H0*c"
    return [xpoints,forcePoints/convFac]


def calcGreensNewtonField(densProf,R_eff,F0=1e-6,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    a=1/(1+z)
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    #xpoints=np.linspace(0,R_eff_final,400)
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    fieldPoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        "Out integrates over shells exterior to the location r"
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*GplusNewton(r,float(R),R_eff,a,rho),dx,r-dx,limit=Limit)[0]
        Out=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*GminusNewton(r,float(R),R_eff,a,rho),r+dx,5,limit=Limit)[0]
        fieldPoints[i]=In+Out
    fieldPoints[0]=fieldPoints[1]
    "fieldPoints has units of (km/s)^2, to convert back to a dimensionless field, divide by (2.997e5)**2"
    return [xpoints,fieldPoints]


def calcGreensNewtonForceR1(densProf,F0=1e-6,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    "densProf must be weighted by R_eff^1"
    a=1/(1+z)
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    #xpoints=np.linspace(0,R_eff_final,400)
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    
    forcePoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*FplusNewtonR1(r,float(R),a,rho),dx,r-dx,limit=Limit)[0]
        forcePoints[i]=In
    
    convFac=29979245.8 
    "dividing by convFac changes units from h/Mpc (km/s)^2 to H0*c"
    "Force is returned with units of H0*c"
    return [xpoints,forcePoints/convFac]

def calcGreensNewtonFieldR2(densProf,z=0,Limit=10000,dx=0.00001,outTo=5,omega_M=0.281,omega_L=0.719):
    "densProf must be weighted by R_eff^2"
    a=1/(1+z)
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    #xpoints=np.linspace(0,R_eff_final,400)
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    
    fieldPoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        "In integrates over shells interior to the location r"
        "Out integrates over shells exterior to the location r"
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*GplusNewtonR2(r,float(R),a,rho),dx,r-dx,limit=Limit)[0]
        Out=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*GminusNewtonR2(r,float(R),a,rho),r+dx,outTo,limit=Limit)[0]
        fieldPoints[i]=In+Out
    "fieldPoints has units of (km/s)^2, to convert back to a dimensionless field, divide by (2.997e5)**2"
    fieldPoints[0]=fieldPoints[1]
    return [xpoints,fieldPoints]


"Below methods are used in the iterative algorithm which calulates the full fifth force"

def fullEqnCheckInterp(forcePoints,fieldPoints,R_eff,F0=1e-6,z=0,outTo=5,omega_M=0.281,omega_L=0.719): 
    "when entered, forcePoints should have units of H0*c and fieldpoint should have units of (km/s)^2"
    
    a=1/(1+z)
    H0=100
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    
    if (len(forcePoints) !=2):
        Q=np.array([np.arange(0.02,outTo,.02),forcePoints[1:]])
        forcePoints=Q
    else:
        Q=np.array([forcePoints[0][1:],forcePoints[1][1:]])
        forcePoints=Q
    convFac=29979245.8 
    Force=convFac*forcePoints[1]
    Force[0]=0
    "modFac turns the fields average value today into the fields average value at redshift z"
    modFac=((1+4*omega_L/omega_M)/(a**(-3)+4*omega_L/omega_M))**2

    "forcePoints now has units of h/Mpc (km/s)^2"
    forceInterp = sp.interpolate.CubicSpline(forcePoints[0],Force,bc_type = 'not-a-knot')
    fieldInterp = sp.interpolate.CubicSpline(fieldPoints[0],fieldPoints[1],bc_type = 'not-a-knot')

    "field still has units of (km/s)^2"
    
    def lap(x,dx=0.0001):
        return a/(R_eff*(x**2))*(((2*forceInterp(x+dx)+0)*((x+dx)**2))-((2*forceInterp(x-dx)+0)*((x-dx)**2)))/(2*dx)
    def term2(x):
        return omega_M*((a*H0)**2)*((1+4*omega_L/omega_M)*(-F0/(-F0*modFac+(fieldInterp(x)/((2.997e5)**2))))**(1/2)-(a**(-3)+4*omega_L/omega_M))
    
    delta=[-(lap(x)-term2(x))/(8/3*np.pi*G*(a**2)*rho) for x in Q[0]]
    
    return [Q[0],delta]


def forceFromField(fieldPoints,R_eff,z=0,dx=0.00001):
    "fieldPoints should have units of (km/s)^2 when entered here"
    a=1/(1+z)
    convFac=29979245.8
    xpoints=fieldPoints[0]
    forcePoints=np.array([xpoints,1/(2*R_eff*a)*np.gradient(fieldPoints[1],fieldPoints[0])/convFac])
    "forcePoints now has units of Ho*c"
    return forcePoints
    

def fifthFieldCorrection(field,R_eff,z,Delta,F0=1e-6,Limit=10000,dx=0.00001,outTo=5,fracBeforeMod=0.8,omega_M=0.281,omega_L=0.719):
    "Delta should be equal to the real density minus reconstructed density --> delta0-delta1, equal to epsilon from the paper"
    "delta1 is reconstructed from fullEqnCheckInterp above"
    "field should have units of (km/s)^2 here"
    
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    a=1/(1+z)
    
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)

    
    modFac=((1+4*omega_L/omega_M)/(a**(-3)+4*omega_L/omega_M))**2
    "modFac turns the fields average value today into the fields average value at redshift z"
    
    fieldInterp=sp.interpolate.CubicSpline(field[0],field[1],bc_type = 'not-a-knot')
    

    DeltaInterp = sp.interpolate.CubicSpline(Delta[0],Delta[1],bc_type = 'not-a-knot')
    dimFac=((2.997e5)**2)
    "multiplying by dimfac changes a dimensionless field into one with units of (km/s)^2"
    fieldMax=dimFac*F0*modFac
    x=fracBeforeMod
    
    m=omega_M**.5*(100)
    R_bar0=3*(m**2)*(1.+4.*(omega_L/omega_M))
    
    def MU(r):
        f=fieldInterp(r)+0
        if (f>x*fieldMax):
            "Check to see if field values are physical. if not we smoothly modify the field to keep it physical"
            s=f-x*fieldMax
            A=(1-x)*fieldMax
            g=A*s/(A+s)
            f=x*fieldMax+g
        
        return np.sqrt((dimFac*F0/(dimFac*F0*modFac-f))**(3/2)*R_bar0/(6*F0))/299792.458
        
    "to speed up algorithm, change the number of points in xpoints below to a smaller number, making sure to still span the range 0 to 'outTo'. be careful not use too few"
    xpoints=np.arange(0,outTo,.05)
    fieldCorrectionPoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        In=sp.integrate.quad(lambda R: (DeltaInterp(R)+0)*Gplus(r,float(R),R_eff,a,MU(float(R)),rho),dx,r-dx,limit=Limit)[0]
        Out=sp.integrate.quad(lambda R: (DeltaInterp(R)+0)*Gminus(r,float(R),R_eff,a,MU(float(R)),rho),r+dx,xpoints[-1],limit=Limit)[0]

        fieldCorrectionPoints[i]=In+Out
    print("field changed by a maximum of "+ str(np.max(np.abs(fieldCorrectionPoints))) + " (km/s)^2")    
    
    fieldCorrectionPoints[0]=fieldCorrectionPoints[1]
    correctInterp=sp.interpolate.CubicSpline(xpoints,fieldCorrectionPoints,bc_type = 'not-a-knot')
    xpoints=np.arange(0,outTo,.02)
    return np.array([xpoints,[correctInterp(x) for x in xpoints]])
    


def fieldCheck(field,fieldMax,fracBeforeMod=0.8):
    "makes sure in a smooth way that all of field's entries are below fieldMax"
    x=fracBeforeMod
    j=0
    for i in range(len(field[1])):
        f=field[1][i]
        if (f>x*fieldMax):
            if (j==0):
                print("field was modified by fieldCheck")
                j=1
            s=f-x*fieldMax
            A=(1-x)*fieldMax
            g=A*s/(A+s)
            field[1][i]=x*fieldMax+g
    return field 
        

def getInitialFieldGuess(densProf,z,F0=1e-6,outTo=5,omega_M=0.281,omega_L=0.719):
    "used to guess a 'fully screened' field when the initial linear guess is too far off to be used. Only implemented in rare cases"
    "needs to return a field value with units of (km/s)^2"
    
    G=4.301*10**-9 
    "G=newtons constant with units of (km/s)^2 (Mpc/MSun)"
    a=1/(1+z)
    rho_bar0=3*(100)**2/(8*np.pi*G)*omega_M
    rho=rho_bar0/(a**3)
    
    delta=densProf[1]
    xPoints=densProf[0]
    delta=densProf[1]
    xpoints=np.arange(0,outTo,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    
    dimFac=(2.997e5)**2
    "multiply a dimensionless field by dimFac to give it units of (km/s)**2"
    
    m=omega_M**.5*(100)
    R_bar=3*(m**2)*(a**(-3)+4.*(omega_L/omega_M)) # units of h**2(km/s)^2/(Mpc)^2
    R_bar0=3*(m**2)*(1.+4.*(omega_L/omega_M))
    fRbar=F0*((R_bar0/R_bar)**2)
    fR=np.array([-dimFac*(F0*((R_bar0/(8*np.pi*G*rho*deltaInterp(r)+R_bar))**2)-fRbar) for r in xpoints])
    field=np.array([xpoints,fR])
    return field 




"The below methodis used to calculate the newtonian force along with the standard error. the density profile must be weighted by the effective radius"
"still need to go through this "
def calcGreensNewtonForceError(densProf,F0=1e-6,z=0,Limit=10000,dx=0.00001):
    "densProf must be weighted by R_eff^1, and include the standard error in the mean "
    a=1/(1+z)
    mu=(1/2997.92458*(1/(2*F0))**(1/2))*((omega_M*a**(-3)+4*omega_L)**(3/2))/(omega_M+4*omega_L)

    rho=rho_bar0/(a**3)
    xPoints=densProf[0]
    delta=densProf[1]
    
    xpoints=np.arange(0,3,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    forcePoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*FplusNewtonR1(r,float(R),a,mu,rho),dx,r-dx,limit=Limit)[0]
        forcePoints[i]=In
    
    
    delta=densProf[2]
    #xpoints=np.linspace(0,R_eff_final,400)
    xpoints=np.arange(0,3,.02)
    deltaInterp = sp.interpolate.CubicSpline(xPoints,delta,bc_type = 'natural')
    errorPoints=np.zeros(len(xpoints))
    for i in range(1,len(xpoints)):
        r=xpoints[i]
        In=sp.integrate.quad(lambda R: (deltaInterp(R)+0)*FplusNewtonR1(r,float(R),a,mu,rho),dx,r-dx,limit=Limit)[0]
        errorPoints[i]=-In
    
    return [xpoints,forcePoints/convFac,errorPoints/convFac]








