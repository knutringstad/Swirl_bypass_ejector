
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import diff, sort
from matplotlib import cm
import pickle

file_name="tangential_vel_yplusminus"
# file_name="mixing_chamber_BL_zdir"
# file_name="mixing_chamber_BL_bypass_ydir"

ydir = True
if(file_name.endswith("ydir")) :
    ydir = True

with open(file_name) as f:
    lines = f.readlines()

pvec=[]
uvec=[]
p_=[]
u_=[]
pos_x_vec=[]
dir_vec=[]

for line in lines[2:]:

    if line.startswith("(("):
        BL=line.split('\"')[1]
        posx = BL.split("_")[2]
        pos_x_vec.append(float(posx))
        if ydir==True:
            dir  = BL.split("_")[3]
        else:
            dir  = BL.split("_")[4]

        dir_vec.append(dir)
        if not p_==[]:
            p_, u_ = zip(*sorted(zip(p_, u_)))
            pvec.append(p_)
            uvec.append(u_)
        p_=[]
        u_=[]
        
    elif not line.startswith("\n") and not line.startswith("(") and not line.startswith(")"):
        p,u = [float(x) for x in line.split("\t")]
        p_.append(p)
        u_.append(u)

pvec.append(p_)
uvec.append(u_)


delta_w_vec_zplus=[]
delta_w_vec_yplus=[]
delta_w_vec_zminus=[]
delta_w_vec_yminus=[]


for i,position in enumerate(pvec):
    u=np.array(uvec[i])
    p=np.array(pvec[i])
    direction = dir_vec[i]
    umax =max(u)
    arg =np.where( np.absolute(p)  == np.max(np.absolute(p)) )[0][0]
    usuc =  u[ arg ]

    if any(diff(np.array(p))==0):
        u=u[1:]
        u = u[    np.absolute( diff(p) )  >1e-6]
        p_l=p[1:]
        p = p_l[    np.absolute( diff(p) )  >1e-6]

        
    dudr =  diff(np.array(u)) /diff(np.array(p))


    maxdudr= np.max(np.absolute(dudr))

    delta_w = (umax-usuc)/maxdudr

    nondimcolorval = (pos_x_vec[i]-min(pos_x_vec))/(max(pos_x_vec)-min(pos_x_vec))

    if dir_vec[i] == "z+":
        delta_w_vec_zplus.append(delta_w)
        plt.plot(p,u,c=cm.winter(nondimcolorval),label='_nolegend_')
    elif dir_vec[i] == "z-":
        delta_w_vec_zminus.append(delta_w)
        plt.plot(p,u,c=cm.winter(nondimcolorval))
    elif dir_vec[i] == "y+":
        delta_w_vec_yplus.append(delta_w)
        plt.plot(p,u,c=cm.winter(nondimcolorval),label='_nolegend_')
    elif dir_vec[i] == "y-":
        delta_w_vec_yminus.append(delta_w)
        plt.plot(p,u,c=cm.winter(nondimcolorval))
    

plt.xlabel("Radial coordinate - r [m]")
plt.ylabel("Tangential velocity - u [m/s]")

legend_string=[]
for i,n in enumerate(pos_x_vec[::2]):
    legend_string.append( "%.2f" %((n-min(pos_x_vec))/0.0031) )

plt.legend(legend_string,ncol=2,title='$x/D_{mix}$')


pos_x_vec_2 = pos_x_vec[::2]
posx_vec_nondim = np.divide(np.array(pos_x_vec_2)-min(np.array(pos_x_vec_2)),0.0031) 



# f = open('store1.pckl', 'wb')
# pickle.dump(delta_w_vec_zplus, f)
# f.close()
# f = open('store2.pckl', 'wb')
# pickle.dump(delta_w_vec_zminus, f)
# # f.close()
# f = open('store3.pckl', 'wb')
# pickle.dump(delta_w_vec_yminus, f)
# f.close()
# f = open('store4.pckl', 'wb')
# pickle.dump(delta_w_vec_yplus, f)
# f.close()

# f = open('store1.pckl', 'rb')
# delta_w_vec_zplus = pickle.load(f)
# f.close()
# f = open('store2.pckl', 'rb')
# delta_w_vec_zminus = pickle.load(f)
# f.close()
# f = open('store3.pckl', 'rb')
# delta_w_vec_yminus = pickle.load(f)
# f.close()
# f = open('store4.pckl', 'rb')
# delta_w_vec_yplus = pickle.load(f)
# f.close()


plt.figure()

# plt.plot(posx_vec_nondim,delta_w_vec_zminus,marker="*")
# plt.plot(posx_vec_nondim,delta_w_vec_yminus,marker="*")
# plt.plot(posx_vec_nondim,delta_w_vec_zplus,marker="*")
# plt.plot(posx_vec_nondim,delta_w_vec_yplus,marker="*")

plt.xlabel("Axial coordinate - $x/D_{mix}$")
plt.ylabel("Mixing layer thickness, Eqn. (6) - $\delta_w$ [m]")

plt.legend([r"Negative z - direction, $\theta_{bp}=0^\circ$",r"Negative y - direction, $\theta_{bp}=90^\circ$",r"Positive z - direction, $\theta_{bp}=180^\circ$" ,r"Positive y - direction, $\theta_{bp}=270^\circ$"])




plt.show()
