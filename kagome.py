import math
from math import pi
from math import sin as sin
from math import cos as cos
mol = {}
hyd= {}
hxlat = {}
hxlat[0,0,0] = [0,0,0]
coord = {}
hydro = {}
a = 15
CC=1.42
CH = 1.08
horde = 71
const1 = 2.88
pd = math.sqrt(3)*a
rpd = math.sqrt(3)*CC
twist = -75
for cx in range(0,8):
    for cy in range(0,8):
        hxlat[cx,cy,0] = [pd*(cx+0.5*cy),1.5*a*cy,0]
        hxlat[cx,cy,1] = [pd*(cx+0.5*cy)+pd/4,1.5*a*cy+0.866*pd/2,0]
        hxlat[cx,cy,2] = [pd*(cx+0.5*cy)+pd/2,1.5*a*cy,0]    
def turn(i,alpha,upperend):
    for z in range(1,upperend):
        mol[i,z] = [cos(alpha/180*pi)*mol[i,z][0]+sin(alpha/180*pi)*mol[i,z][1],-sin(alpha/180*pi)*mol[i,z][0]+cos(alpha/180*pi)*mol[i,z][1],0]
def turnhydro(i,alpha,upperend):
    for z in range(1,upperend):
        hyd[i,z] = [cos(alpha/180*pi)*hyd[i,z][0]+sin(alpha/180*pi)*hyd[i,z][1],-sin(alpha/180*pi)*hyd[i,z][0]+cos(alpha/180*pi)*hyd[i,z][1],0]    
def shift(i,xshift,yshift,upperend):        
    for z in range(1,upperend):    
        mol[i,z] = [xshift+mol[i,z][0],yshift+mol[i,z][1],0]    
def shifthydro(i,xshift,yshift,upperend):        
    for z in range(1,upperend):    
        hyd[i,z] = [xshift+hyd[i,z][0],yshift+hyd[i,z][1],0] 
def singletetracene(j):
    mol[j,1] = [0,+0.5* CC,0]
    mol[j,2] = [0,-0.5* CC,0]
    mol[j,3] = [-rpd,0.5* CC,0]
    mol[j,4] = [+rpd,-0.5* CC,0]
    mol[j,5] = [-rpd,-0.5* CC,0]
    mol[j,6] = [rpd,0.5* CC,0]
    mol[j,7] = [-2*rpd,0.5* CC,0]
    mol[j,8] = [-2*rpd,-0.5* CC,0]
    mol[j,9] = [2*rpd,-0.5* CC,0]
    mol[j,10] = [2*rpd,0.5* CC,0]
    mol[j,11] = [-0.5*rpd,1* CC,0]
    mol[j,12] = [-0.5*rpd,-1* CC,0]
    mol[j,13] = [0.5*rpd,-1* CC,0]
    mol[j,14] = [0.5*rpd,1* CC,0]
    mol[j,15] = [-1.5*rpd,1* CC,0]
    mol[j,16] = [-1.5*rpd,-1* CC,0]
    mol[j,17] = [1.5*rpd,-1* CC,0]
    mol[j,18] = [1.5*rpd, 1* CC,0]    
def hydrogenedge(j):
    hyd[j,1] = [0.5*rpd,CC+CH,0]
    hyd[j,2] = [-0.5*rpd,CC+CH,0]
    hyd[j,3] = [-0.5*rpd,-CC-CH,0]
    hyd[j,4] = [0.5*rpd,-CC-CH,0]
    hyd[j,5] = [1.5*rpd,CC+CH,0]
    hyd[j,6] = [-1.5*rpd,CC+CH,0]
    hyd[j,7] = [-1.5*rpd,-CC-CH,0]
    hyd[j,8] = [1.5*rpd,-CC-CH,0]
    hyd[j,9] = [2*rpd+0.866*CH,0.5*CC+0.5*CH,0]
    hyd[j,10] = [-2*rpd-0.866*CH,0.5*CC+0.5*CH,0]
    hyd[j,11] = [-2*rpd-0.866*CH,-0.5*CC-0.5*CH,0]
    hyd[j,12] = [2*rpd+0.866*CH,-0.5*CC-0.5*CH,0]    
def maketetracene(xcoord,ycoord,i):
    singletetracene(i)
    turn(i,twist,19)
    shift(i,xcoord,ycoord,19)
def makehydrotetracene(xcoord,ycoord,i):
    hydrogenedge(i)
    turnhydro(i,twist,13)
    shifthydro(i,xcoord,ycoord,13)  
def makeROTtetracene(xcoord,ycoord,i):
    singletetracene(i)
    turn(i,120,19)
    turn(i,twist,19)
    shift(i,xcoord,ycoord,19)    
def makeROThydrotetracene(xcoord,ycoord,i):
    hydrogenedge(i)
    turnhydro(i,120,13)
    turnhydro(i,twist,13)
    shifthydro(i,xcoord,ycoord,13)        
def makeTORtetracene(xcoord,ycoord,i):
    singletetracene(i)
    turn(i,-120,19)
    turn(i,twist,19)
    shift(i,xcoord,ycoord,19)           
def makeTORhydrotetracene(xcoord,ycoord,i):
    hydrogenedge(i)
    turnhydro(i,-120,13)
    turnhydro(i,twist,13)
    shifthydro(i,xcoord,ycoord,13)     
enum = 0
for key in hxlat:
        if key[2] == 0:
            enum = enum +1
            maketetracene(hxlat[key][0],hxlat[key][1],enum)
            makehydrotetracene(hxlat[key][0],hxlat[key][1],enum)
        if key[2] == 1:
            enum = enum +1
            makeROTtetracene(hxlat[key][0],hxlat[key][1],enum) 
            makeROThydrotetracene(hxlat[key][0],hxlat[key][1],enum)
        if key[2] == 2:
            enum = enum +1
            makeTORtetracene(hxlat[key][0],hxlat[key][1],enum) 
            makeTORhydrotetracene(hxlat[key][0],hxlat[key][1],enum)        
#string = str(len(mol)+len(hyd)+horde*horde)+"\n\n"        
string = str(len(mol)+len(hyd))+"\n\n"  
"""
for key in hxlat:
        string += "Au\t"
        string += str(hxlat[key][0])
        string += "\t"
        string += str(hxlat[key][1])
        string += "\t"
        string += str(hxlat[key][2])
        string += "\n"
"""        
for key in mol:
        string += "C\t"
        string += str(mol[key][0])
        string += "\t"
        string += str(mol[key][1])
        string += "\t"
        string += str(mol[key][2])
        string += "\n"   
        
for key in hyd:
        string += "H\t"
        string += str(hyd[key][0])
        string += "\t"
        string += str(hyd[key][1])
        string += "\t"
        string += str(hyd[key][2])
        string += "\n"     
"""
for a in range(0,horde):
    for b in range(-horde,0):
        string += "Au\t"
        string += str(a*const1/2-b*const1/2)
        string += "\t"
        string += str(math.sqrt(3)/2*a*const1+math.sqrt(3)/2*b*const1)
        string += "\t"
        string += str(-3.2)
        string += "\n"   
"""
filename = "kagomedense.xyz"
with open(filename,"w") as f:
    f.write(string)
