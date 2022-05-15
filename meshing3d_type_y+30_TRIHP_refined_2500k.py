


from math import sqrt
from numpy.core.fromnumeric import compress
from numpy.lib.function_base import angle


def meshing3D_test():
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import subprocess

    SampleName = 'geometry_parameters_TRIHP.dat' 
    A=[]
    f=open(SampleName, "r")
    for line in f:
        A.append(float(line))



    write_mesh=False

    f.close()
    Lmotive = A[0]
    DmotiveIn = A[1]
    DmotiveOut = A[2]
    Dthroat = A[3]
    alphaMotiveDiff = A[4]
    alphaMotiveConv = A[5]
    ThicknessNozzle =A[6]
    Lmch = A[7]
    alphaSuction = A[8]
    Dmix = A[9]
    Lmix = A[10]
    alphadiff = A[11]
    DdiffOut = A[12]
    Loutlet = A[13]

    ThicknessNozzleOuter = ThicknessNozzle*2  #THIS COULD BE CHANGED UP
    Ldiff = (DdiffOut/2-Dmix/2)/(math.tan(math.radians(alphadiff/2)))
    Lthroat = Lmotive - (DmotiveOut/2-Dthroat/2)/(math.tan(math.radians(alphaMotiveDiff/2)))
    L1 = Lthroat - (-Dthroat/2+DmotiveIn/2)/(math.tan(math.radians(alphaMotiveConv/2)))

    if L1 <0:
        print("in meshing")
        print(Lmotive)
        print(Lthroat)
        print(L1)
        raise Exception("error in geom setup: L1<0")


    x   = np.ndarray(10,float)
    y1  = np.ndarray(10,float)
    y2  = np.ndarray(10,float)
    y3  = np.ndarray(10,float)
    yi  = np.ndarray(10,float)

    alphaMotiveOuter = alphaSuction
    xmOutConv = Lmotive-(DmotiveIn/2-DmotiveOut/2+ThicknessNozzleOuter-ThicknessNozzle)/math.tan(math.radians((alphaMotiveOuter/2)))

    xsConv = xmOutConv*0.9
    Dsuc = 2*(math.tan(math.radians(alphaSuction/2))*(Lmotive+Lmch-xsConv) + Dmix/2)
    
    x[0]= 0
    x[1]= L1
    x[2]= xsConv
    x[3]= xmOutConv
    x[4]= Lthroat
    x[5]= Lmotive
    x[6]= Lmotive+Lmch
    x[7]= Lmotive+Lmch+Lmix
    x[8]= Lmotive+Lmch+Lmix+Ldiff
    x[9]= Lmotive+Lmch+Lmix+Ldiff+Loutlet

    outerthinningfactor = 1.2

    y1[0]= DmotiveIn/2
    y1[1]= DmotiveIn/2
    y1[2]= DmotiveIn/2 +  (x[1]-x[2])*math.tan(math.radians((alphaMotiveConv/2)))
    y1[3]= DmotiveIn/2 +  (x[1]-x[3])*math.tan(math.radians((alphaMotiveConv/2)))
    y1[4]= Dthroat/2
    y1[5]= DmotiveOut/2 
    y1[6]= DmotiveOut/2 * outerthinningfactor
    y1[7]= DmotiveOut/2 * outerthinningfactor
    y1[8]= y1[7] * DdiffOut /Dmix 
    y1[9]= y1[7] * DdiffOut /Dmix 
    
    suctionTipThinningFactor = 1.4

    y2[0]= DmotiveIn/2 +ThicknessNozzleOuter
    y2[1]= DmotiveIn/2 +ThicknessNozzleOuter
    y2[2]= DmotiveIn/2 +ThicknessNozzleOuter
    y2[3]= DmotiveIn/2 +ThicknessNozzleOuter
    y2[4]= DmotiveIn/2 +ThicknessNozzleOuter -  (x[4]-x[3])*math.tan(math.radians((alphaMotiveOuter/2)))
    y2[5]= DmotiveOut/2 + ThicknessNozzle
    y2[6]= (DmotiveOut/2 + ThicknessNozzle)*suctionTipThinningFactor
    y2[7]= (DmotiveOut/2 + ThicknessNozzle)*suctionTipThinningFactor
    y2[8]= y2[7] * DdiffOut /Dmix
    y2[9]= y2[7] * DdiffOut /Dmix

    y3[0]= Dsuc/2 
    y3[1]= Dsuc/2 
    y3[2]= Dsuc/2 
    y3[3]= Dsuc/2 -  (x[3]-x[2])*math.tan(math.radians((alphaSuction/2)))
    y3[4]= Dsuc/2 -  (x[4]-x[2])*math.tan(math.radians((alphaSuction/2)))
    y3[5]= Dsuc/2 -  (x[5]-x[2])*math.tan(math.radians((alphaSuction/2)))
    y3[6]= Dmix/2
    y3[7]= Dmix/2
    y3[8]= DdiffOut/2
    y3[9]= DdiffOut/2

    superinnerthinningfactor = 0.8

    yi[0]= y1[0]*0.6
    yi[1]= y1[1]*0.6
    yi[2]= y1[2]*0.6
    yi[3]= y1[3]*0.6
    yi[4]= y1[4]*0.6
    yi[5]= y1[5]*0.6
    yi[6]= y1[6]*0.6*superinnerthinningfactor
    yi[7]= y1[7]*0.6*superinnerthinningfactor
    yi[8]= y1[8]*0.6*superinnerthinningfactor
    yi[9]= y1[9]*0.6*superinnerthinningfactor


    xb=4*Dmix
    L=0.5*Dmix
    
    x_swirl_1 =Lmotive+xb
    x_swirl_2 =Lmotive+xb+L
    R = Dmix/2
    H = Dmix*2/3


    
    theta = 30

    alpha_i = 0
    alpha_t = 40 #+ theta/2




    
    ScriptName = 'EjectorBypassSwirl_meshing_3D_y+30_F.rpl' 

    fid=open(ScriptName,'w')

    # POINTS
    fid.write('ic_set_global geo_cad 0 toptol_userset\n') 
    fid.write('ic_set_global geo_cad 0.0 toler\n') 
    fid.write('ic_geo_new_family GEOM\n') 
    fid.write('ic_boco_set_part_color GEOM\n') 
    fid.write('ic_empty_tetin\n') 

    j=1  #point index
    k = 1 # curve index
    l = 1

    #curves and points
    for i in range(10):
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],y1[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],-y1[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,-y1[i] ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,y1[i] ))
        j = j +1
        fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,j-5,j-4,j-1,y1[i])  )
        k = k +1

        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],y2[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],-y2[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,-y2[i] ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,y2[i] ))
        j = j +1
        fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,j-9,j-4,j-1,y2[i])  )
        k = k +1

        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],y3[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],-y3[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,-y3[i] ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,y3[i] ))
        j = j +1
        fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,j-13,j-4,j-1,y3[i])  )
        k = k +1

        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],yi[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],-yi[i],0 ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,-yi[i] ))
        j = j +1
        fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x[i],0,yi[i] ))
        j = j +1


    krevolve=k
    for i in range(9):
        fid.write("ic_curve point GEOM crv.%d {pnt.%d pnt.%d}\n" % (k, i*17+12,(i+1)*17+12))
        k=k+1




    y_1 = R*math.sin(math.radians(0))
    y_2 = R*math.sin(math.radians(theta))

    z_1= R*math.cos(math.radians(0))
    z_2= R*math.cos(math.radians(theta))




    
    
    mode = "no offset"
    if mode =="offset":
        y_1 = R*math.sin(math.radians(-theta/2))
        y_2 = R*math.sin(math.radians(theta/2))

        z_1= R*math.cos(math.radians(-theta/2))
        z_2= R*math.cos(math.radians(theta/2))
    
    
    # swirl box

    pi = 3.14159
    compression_ratio_base=1.2
    angle_factor=(theta/180*pi)/(1-math.cos(math.radians(theta)))
    
    compression_ratio=compression_ratio_base
    # L*sqrt(cr) *W*sqrt(cr) 

    cr_displacement_x = (x_swirl_2-x_swirl_1)*(compression_ratio)

    cr_displacement_z = (z_1-z_2)/angle_factor
    # cr_displacement_x=0
    # cr_displacement_z=0

    H_1= H -y_1
    H_2= H -y_2

    delta_x_l_alpha_i = (y_2-y_1)*math.tan(math.radians(alpha_i))
    delta_x_u_alpha_i = H_2*math.tan(math.radians(alpha_i))

    delta_z_1_alpha_t = H
    delta_z_2_alpha_t = H
    if alpha_t >5:
        delta_z_1_alpha_t = H_1*math.tan(math.radians(90-alpha_t))
        delta_z_2_alpha_t = H_2*math.tan(math.radians(90-alpha_t))

    y_outer_1= y_1 + H*math.sin(math.radians(alpha_t)) + (y_2-y_1)*math.sin(math.radians(alpha_t))
    y_outer_2= y_2 + H*math.sin(math.radians(alpha_t))


    # inner
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1 +delta_x_l_alpha_i ,y_1, z_1  )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1 ,y_2, z_2  )    )
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2 +delta_x_l_alpha_i,y_1, z_1  )    )
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2 ,y_2, z_2  )    )
    j = j +1

    
    # outer
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1-delta_x_u_alpha_i -cr_displacement_x/2 ,y_outer_1, z_1 +delta_z_1_alpha_t + cr_displacement_z  )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1-delta_x_u_alpha_i-cr_displacement_x/2 ,y_outer_2, z_2 +delta_z_2_alpha_t - cr_displacement_z)    )
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2-delta_x_u_alpha_i+cr_displacement_x/2,y_outer_1, z_1 +delta_z_1_alpha_t + cr_displacement_z)    )
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2-delta_x_u_alpha_i+cr_displacement_x/2,y_outer_2, z_2 +delta_z_2_alpha_t - cr_displacement_z )    )
    j = j +1

    # center
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1,0, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2,0, 0 )  )    
    j = j +1

    #circle support (for slanted inlets)
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1,Dmix/2, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1,0, Dmix/2 )  )    
    j = j +1

    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2,Dmix/2, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2,0, Dmix/2 )  )    
    j = j +1




    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1+delta_x_l_alpha_i,0, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2+delta_x_l_alpha_i,0, 0)  )    
    j = j +1

    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1+delta_x_l_alpha_i,Dmix/2, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2+delta_x_l_alpha_i,Dmix/2, 0 )  )    
    j = j +1

    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_1,-Dmix/2, 0 )  )    
    j = j +1
    fid.write("ic_point {{}} GEOM pnt.%d {%f %f %f}\n" % (j, x_swirl_2,-Dmix/2, 0 )  )     
    j = j +1
    


    fid.write("ic_curve point GEOM crv.%d {pnt.172 pnt.174}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.171 pnt.173}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.172 pnt.171}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.174 pnt.173}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.174 pnt.178}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.178 pnt.177}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.177 pnt.173}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.178 pnt.176}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.177 pnt.175}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.175 pnt.176}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.176 pnt.172}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.171 pnt.175}\n" %k)
    k= k+1
    

    fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,179,181,182,Dmix/2)  )
    k= k+1
    fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,180,183,184,Dmix/2)  )
    k= k+1





    fid.write("ic_curve point GEOM crv.%d {pnt.6 pnt.23}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.23 pnt.40}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.40 pnt.57}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.57 pnt.74}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.74 pnt.91}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.2 pnt.19}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.19 pnt.36}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.36 pnt.53}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.53 pnt.70}\n" %k)
    k= k+1
    fid.write("ic_curve point GEOM crv.%d {pnt.70 pnt.87}\n" %k)
    k= k+1


    fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,185,171,187,Dmix/2)  )
    k= k+1
    fid.write("ic_curve arc_ctr_rad GEOM crv.%d {pnt.%d pnt.%d pnt.%d %f 0 360}\n" %(k,186,173,188,Dmix/2)  )
    k= k+1

    fid.write("ic_curve point GEOM crv.%d {pnt.171 pnt.189}\n" %k)
    k= k+1

    
    fid.write("ic_curve point GEOM crv.%d {pnt.173 pnt.190}\n" %k)
    k= k+1


    


    

    fid.write("""
    ic_undo_group_end
    ic_undo_group_begin 
    ic_geo_new_family FLUID
    ic_boco_set_part_color FLUID
    ic_hex_initialize_blocking {} FLUID 0 101
    ic_hex_unblank_blocks 
    ic_hex_multi_grid_level 0
    ic_hex_projection_limit 0
    ic_hex_default_bunching_law default 2.0
    ic_hex_floating_grid off
    ic_hex_transfinite_degree 1
    ic_hex_unstruct_face_type one_tri
    ic_hex_set_unstruct_face_method uniform_quad
    ic_hex_set_n_tetra_smoothing_steps 20
    ic_hex_error_messages off_minor
    ic_undo_group_end 
    ic_hex_mark_blocks unmark
    ic_undo_group_begin 
    ic_hex_mark_blocks superblock 13
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks face_neighbors corners { 21 25 22 26 } { 37 41 38 42 }
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_ogrid 1 m GEOM FLUID -version 50
    ic_hex_mark_blocks unmark
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks superblock 13
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks face_neighbors corners { 68 70 69 71 } { 72 74 73 75 }
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_ogrid 1 m GEOM FLUID -version 50
    ic_hex_mark_blocks unmark
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks superblock 13
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks face_neighbors corners { 84 86 85 87 } { 88 90 89 91 }
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_ogrid 1 m GEOM FLUID -version 50
    ic_hex_mark_blocks unmark
    ic_undo_group_end 
    ic_hex_mark_blocks unmark\n""")

    fid.write("""ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 25 41 pnt.27 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 121 41 pnt.44 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 149 41 pnt.61 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 177 41 pnt.95 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 205 41 pnt.112 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 233 41 pnt.127 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 261 41 pnt.146 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 177 205 pnt.78 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end
     \n""")




    fid.write("""
    ic_hex_move_node 70 pnt.8
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 86 pnt.4
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 102 pnt.16
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 68 pnt.7
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 84 pnt.3
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 100 pnt.15
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 101 pnt.17
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 85 pnt.5
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 69 pnt.9
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 71 pnt.6
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 87 pnt.2
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 103 pnt.14
    ic_undo_group_end 
    ic_undo_group_begin
    ic_hex_move_node 26 pnt.10
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 25 pnt.12
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 21 pnt.11
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 22 pnt.13
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_align_vertices 0 1 1 0 0 -csys global GEOM FLUID
    ic_undo_group_end  \n""")


    list_of_points = [11, 13, 12, 10, 7, 3, 15, 9, 5, 17, 8, 4, 16, 6 , 2, 14 ]

    for ii in range(2,11):
        counter=0

        for jj in range(1,3):
            for kk in range(1,3):
                point = 17*(ii-1) + list_of_points[counter]
                fid.write("ic_hex_move_node [ic_hex_vertex_number { %d %d %d }] pnt.%d\n" %(ii,jj,kk,point))
                counter=counter+1

        for jj in range(1,3):
            for kk in range(1,3):
                for ll in range(1,4):
                    point = 17*(ii-1) + list_of_points[counter]
                    counter=counter+1
                    fid.write("ic_hex_move_node [ic_hex_vertex_number { %d %d %d 3:%d }] pnt.%d\n" %(ii,jj,kk,ll,point))

    

    # Swirl part

    fid.write("""ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 243 271 pnt.179 m GEOM FLUID 
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 355 271 pnt.180 m GEOM FLUID
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_mark_blocks unmark
    ic_hex_mark_blocks superblock 178
    ic_hex_undo_major_start split_grid
    ic_hex_split_grid 351 355 pnt.172 m GEOM FLUID marked
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    \n""")

    fid.write("""
    ic_geo_check_family FLUID
    ic_geo_check_family GEOM
    ic_geo_check_family ORFN
    ic_geo_check_family VORFN
    ic_hex_extrude_faces 1 { 351 379 494 498 } FLUID 0.000754359 -nsub 1 -type original FLUID GEOM -version 50
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 534 pnt.176
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 538 pnt.178
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 532 pnt.175
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 536 pnt.177
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 494 pnt.172
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_hex_move_node 498 pnt.174
    ic_undo_group_end 
    \n""")

    #ic_undo_group_begin 
    # ic_hex_move_node 570 pnt.176
    # ic_undo_group_end 
    # ic_undo_group_begin 
    # ic_hex_move_node 575 pnt.178
    # ic_undo_group_end 
    # ic_undo_group_begin 
    # ic_hex_move_node 568 pnt.175
    # ic_undo_group_end 
    # ic_undo_group_begin 
    # ic_hex_move_node 573 pnt.177
    # ic_undo_group_end 


    
    fid.write("""
    ic_undo_group_begin
    ic_surface bsinterp GEOM srf.00 crv.1
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_begin
    ic_surface 2-4crvs GEOM srf.01 {0.000001 {crv.2 crv.1}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end
    ic_undo_group_begin
    ic_surface 2-4crvs GEOM srf.02 {0.000001 {crv.3 crv.2}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_set_global geo_cad 5e-7 toler
    ic_set_global geo_cad 5e-7 toler
    ic_set_global geo_cad 5e-7 toler
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.03 {0.000001 {crv.3 crv.6}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.04 {0.000001 {crv.6 crv.9}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.05 {0.000001 {crv.12 crv.9}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.06 {0.000001 {crv.12 crv.15}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.07 {0.000001 {crv.15 crv.18}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.08 {0.000001 {crv.18 crv.21}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_geo_configure_one_attribute surface shade wire
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.09 {0.000001 {crv.37 crv.21 crv.52}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_geo_configure_one_attribute surface shade wire
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.10 {1e-6 {crv.21 crv.24}}
    ic_reinit_geom_objects 
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.11 {1e-6 {crv.24 crv.27}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.12 {0.000001 {crv.27 crv.30}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    \n""")


    
    fid.write("""
    ic_undo_group_begin
    ic_surface bsinterp GEOM srf.23 crv.28
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_surface 2-4crvs GEOM srf.14 {0.000000001 {crv.51 crv.41 crv.46 crv.48}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.15 {0.000000001 {crv.50 crv.40 crv.44 crv.47}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.16 {0.000000001 {crv.42 crv.50 crv.49 crv.51}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.17 {0.000000001 {crv.43 crv.44 crv.46 crv.45}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.18 {0.000000001 {crv.47 crv.49 crv.48 crv.45}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.19 {0.000000001 {crv.24 crv.27}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.20 {0.000000001 {crv.27 crv.30}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.21 {0.000000001 {crv.30 crv.29}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.22 {0.000000001 {crv.29 crv.28}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    \n""")

    
    


    fid.write("""
    ic_set_global geo_cad 5e-7 toler
    ic_set_global geo_cad 5e-7 toler
    ic_set_global geo_cad 5e-7 toler
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.23 {0.0001 {crv.2 crv.5}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.24 {0.0001 {crv.5 crv.8}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.25 {0.0001 {crv.8 crv.11}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.26 {0.0001 {crv.11 crv.14}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.27 {0.0001 {crv.14 crv.17}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.28 {0.0001 {crv.1 crv.4}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.29 {0.0001 {crv.4 crv.7}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.30 {0.0001 {crv.7 crv.10}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.31 {0.0001 {crv.10 crv.13}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.32 {0.0001 {crv.13 crv.16}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.33 {0.0001 {crv.17 crv.16}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_surface 2-4crvs GEOM srf.34 {0.0001 {crv.15 crv.18}}
    ic_set_global geo_cad 5e-7 toler
    ic_set_dormant_pickable point 0 {}
    ic_set_dormant_pickable curve 0 {}
    ic_undo_group_end
    ic_geo_intersect_surfaces GEOM { srf.10 srf.17}
    ic_geo_intersect_surfaces GEOM { srf.10 srf.16} 
    \n""")

    
    fid.write("ic_geo_project_curve_to_surface crv.%d srf.10 crv.%d GEOM 0 0\n" %(k-1,k))
    k= k+1
    fid.write("ic_geo_project_curve_to_surface crv.%d srf.10 crv.%d GEOM 0 0\n" %(k-3,k))
    k= k+1





    
    

    fid.write("""
    ic_hex_find_comp_curve crv.4
    ic_undo_group_begin 
    ic_hex_set_edge_projection 70 71 0 1 crv.2
    ic_hex_set_edge_projection 411 70 0 1 crv.2
    ic_hex_set_edge_projection 68 411 0 1 crv.2
    ic_hex_set_edge_projection 68 69 0 1 crv.2
    ic_hex_set_edge_projection 69 467 0 1 crv.2
    ic_hex_set_edge_projection 467 71 0 1 crv.2
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.1
    ic_undo_group_begin 
    ic_hex_set_edge_projection 86 87 0 1 crv.1
    ic_hex_set_edge_projection 468 87 0 1 crv.1
    ic_hex_set_edge_projection 84 85 0 1 crv.1
    ic_hex_set_edge_projection 85 468 0 1 crv.1
    ic_hex_set_edge_projection 84 412 0 1 crv.1
    ic_hex_set_edge_projection 412 86 0 1 crv.1
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.5
    ic_undo_group_begin 
    ic_hex_set_edge_projection 122 132 0 1 crv.5
    ic_hex_set_edge_projection 415 122 0 1 crv.5
    ic_hex_set_edge_projection 471 132 0 1 crv.5
    ic_hex_set_edge_projection 128 471 0 1 crv.5
    ic_hex_set_edge_projection 118 128 0 1 crv.5
    ic_hex_set_edge_projection 118 415 0 1 crv.5
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.2
    ic_hex_find_comp_curve crv.4
    ic_undo_group_begin 
    ic_hex_set_edge_projection 416 123 0 1 crv.4
    ic_hex_set_edge_projection 123 133 0 1 crv.4
    ic_hex_set_edge_projection 472 133 0 1 crv.4
    ic_hex_set_edge_projection 129 472 0 1 crv.4
    ic_hex_set_edge_projection 119 129 0 1 crv.4
    ic_hex_set_edge_projection 119 416 0 1 crv.4
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.8
    ic_undo_group_begin 
    ic_hex_set_edge_projection 419 150 0 1 crv.8
    ic_hex_set_edge_projection 146 419 0 1 crv.8
    ic_hex_set_edge_projection 146 156 0 1 crv.8
    ic_hex_set_edge_projection 150 160 0 1 crv.8
    ic_hex_set_edge_projection 475 160 0 1 crv.8
    ic_hex_set_edge_projection 156 475 0 1 crv.8
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.7
    ic_undo_group_begin 
    ic_hex_set_edge_projection 420 151 0 1 crv.7
    ic_hex_set_edge_projection 147 420 0 1 crv.7
    ic_hex_set_edge_projection 147 157 0 1 crv.7
    ic_hex_set_edge_projection 476 161 0 1 crv.7
    ic_hex_set_edge_projection 151 161 0 1 crv.7
    ic_hex_set_edge_projection 157 476 0 1 crv.7
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.11
    ic_undo_group_begin 
    ic_hex_set_edge_projection 423 178 0 1 crv.11
    ic_hex_set_edge_projection 174 423 0 1 crv.11
    ic_hex_set_edge_projection 174 184 0 1 crv.11
    ic_hex_set_edge_projection 479 188 0 1 crv.11
    ic_hex_set_edge_projection 184 479 0 1 crv.11
    ic_hex_set_edge_projection 178 188 0 1 crv.11
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.10
    ic_undo_group_begin 
    ic_hex_set_edge_projection 424 179 0 1 crv.10
    ic_hex_set_edge_projection 179 189 0 1 crv.10
    ic_hex_set_edge_projection 175 424 0 1 crv.10
    ic_hex_set_edge_projection 175 185 0 1 crv.10
    ic_hex_set_edge_projection 185 480 0 1 crv.10
    ic_hex_set_edge_projection 480 189 0 1 crv.10
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.15
    ic_undo_group_begin 
    ic_hex_set_edge_projection 426 317 0 1 crv.15
    ic_hex_set_edge_projection 313 426 0 1 crv.15
    ic_hex_set_edge_projection 313 323 0 1 crv.15
    ic_hex_set_edge_projection 317 327 0 1 crv.15
    ic_hex_set_edge_projection 323 482 0 1 crv.15
    ic_hex_set_edge_projection 482 327 0 1 crv.15
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.14
    ic_undo_group_begin 
    ic_hex_set_edge_projection 314 427 0 1 crv.14
    ic_hex_set_edge_projection 427 318 0 1 crv.14
    ic_hex_set_edge_projection 324 483 0 1 crv.14
    ic_hex_set_edge_projection 483 328 0 1 crv.14
    ic_hex_set_edge_projection 318 328 0 1 crv.14
    ic_hex_set_edge_projection 314 324 0 1 crv.14
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.13
    ic_undo_group_begin 
    ic_hex_set_edge_projection 428 319 0 1 crv.13
    ic_hex_set_edge_projection 315 428 0 1 crv.13
    ic_hex_set_edge_projection 315 325 0 1 crv.13
    ic_hex_set_edge_projection 484 329 0 1 crv.13
    ic_hex_set_edge_projection 319 329 0 1 crv.13
    ic_hex_set_edge_projection 325 484 0 1 crv.13
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.16
    ic_undo_group_begin 
    ic_hex_set_edge_projection 432 207 0 1 crv.16
    ic_hex_set_edge_projection 203 432 0 1 crv.16
    ic_hex_set_edge_projection 203 213 0 1 crv.16
    ic_hex_set_edge_projection 488 217 0 1 crv.16
    ic_hex_set_edge_projection 213 488 0 1 crv.16
    ic_hex_set_edge_projection 207 217 0 1 crv.16
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.17
    ic_undo_group_begin 
    ic_hex_set_edge_projection 431 206 0 1 crv.17
    ic_hex_set_edge_projection 206 216 0 1 crv.17
    ic_hex_set_edge_projection 212 487 0 1 crv.17
    ic_hex_set_edge_projection 487 216 0 1 crv.17
    ic_hex_set_edge_projection 202 212 0 1 crv.17
    ic_hex_set_edge_projection 202 431 0 1 crv.17
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.20
    ic_undo_group_begin 
    ic_hex_set_edge_projection 435 234 0 1 crv.20
    ic_hex_set_edge_projection 230 435 0 1 crv.20
    ic_hex_set_edge_projection 230 240 0 1 crv.20
    ic_hex_set_edge_projection 491 244 0 1 crv.20
    ic_hex_set_edge_projection 240 491 0 1 crv.20
    ic_hex_set_edge_projection 234 244 0 1 crv.20
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.42
    ic_hex_find_comp_curve GEOM/1
    ic_undo_group_begin 
    ic_hex_set_edge_projection 351 494 0 1 GEOM/1
    ic_undo_group_end 
    ic_hex_move_node 351 pnt.171
    ic_hex_find_comp_curve crv.69
    ic_undo_group_begin 
    ic_hex_set_edge_projection 341 351 0 1 crv.69
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.43
    ic_undo_group_begin 
    ic_hex_set_edge_projection 379 498 0 1 crv.43
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.68
    ic_undo_group_begin 
    ic_hex_set_edge_projection 369 379 0 1 crv.68
    ic_undo_group_end 
    ic_hex_place_node 379 curve:crv.43 1
    ic_hex_set_edge_projection 351 494 0 1 GEOM/1
    ic_undo_group_end
    ic_hex_find_comp_curve crv.51
    ic_undo_group_begin 
    ic_hex_set_edge_projection 351 532 0 1 crv.51
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.50
    ic_undo_group_begin 
    ic_hex_set_edge_projection 494 534 0 1 crv.50
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.50
    ic_hex_find_comp_curve crv.49
    ic_undo_group_begin 
    ic_hex_set_edge_projection 532 534 0 1 crv.49
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.47
    ic_undo_group_begin 
    ic_hex_set_edge_projection 534 538 0 1 crv.47
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.48
    ic_undo_group_begin 
    ic_hex_set_edge_projection 532 536 0 1 crv.48
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.45
    ic_undo_group_begin 
    ic_hex_set_edge_projection 536 538 0 1 crv.45
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.46
    ic_undo_group_begin 
    ic_hex_set_edge_projection 379 536 0 1 crv.46
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.44
    ic_undo_group_begin 
    ic_hex_set_edge_projection 498 538 0 1 crv.44
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.53
    ic_undo_group_begin 
    ic_hex_set_edge_projection 373 383 0 1 crv.53
    ic_hex_set_edge_projection 369 442 0 1 crv.53
    ic_hex_set_edge_projection 379 498 0 1 GEOM/0
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.30
    ic_undo_group_begin 
    ic_hex_set_edge_projection 41 42 0 1 crv.30
    ic_hex_set_edge_projection 37 454 0 1 crv.30
    ic_hex_set_edge_projection 37 38 0 1 crv.30
    ic_hex_set_edge_projection 38 510 0 1 crv.30
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.69
    ic_undo_group_begin 
    ic_hex_set_edge_projection 341 351 0 1 crv.69
    ic_undo_group_end 
    \n""")



    fid.write("""
    ic_hex_find_comp_curve crv.3
    ic_undo_group_begin 
    ic_hex_set_edge_projection 25 26 0 1 crv.3
    ic_hex_set_edge_projection 21 410 0 1 crv.3
    ic_hex_set_edge_projection 22 466 0 1 crv.3
    ic_hex_set_edge_projection 21 22 0 1 crv.3
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.6
    ic_undo_group_begin 
    ic_hex_set_edge_projection 121 131 0 1 crv.6
    ic_hex_set_edge_projection 117 414 0 1 crv.6
    ic_hex_set_edge_projection 117 127 0 1 crv.6
    ic_hex_set_edge_projection 127 470 0 1 crv.6
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.9
    ic_undo_group_begin 
    ic_hex_set_edge_projection 149 159 0 1 crv.9
    ic_hex_set_edge_projection 145 418 0 1 crv.9
    ic_hex_set_edge_projection 155 474 0 1 crv.9
    ic_hex_set_edge_projection 145 155 0 1 crv.9
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.12
    ic_undo_group_begin 
    ic_hex_set_edge_projection 177 187 0 1 crv.12
    ic_hex_set_edge_projection 173 422 0 1 crv.12
    ic_hex_set_edge_projection 183 478 0 1 crv.12
    ic_hex_set_edge_projection 173 183 0 1 crv.12
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.18
    ic_undo_group_begin 
    ic_hex_set_edge_projection 205 215 0 1 crv.18
    ic_hex_set_edge_projection 201 430 0 1 crv.18
    ic_hex_set_edge_projection 211 486 0 1 crv.18
    ic_hex_set_edge_projection 201 211 0 1 crv.18
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.21
    ic_undo_group_begin 
    ic_hex_set_edge_projection 233 243 0 1 crv.21
    ic_hex_set_edge_projection 229 434 0 1 crv.21
    ic_hex_set_edge_projection 229 239 0 1 crv.21
    ic_hex_set_edge_projection 239 490 0 1 crv.21
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.24
    ic_undo_group_begin 
    ic_hex_set_edge_projection 261 271 0 1 crv.24
    ic_hex_set_edge_projection 257 446 0 1 crv.24
    ic_hex_set_edge_projection 257 267 0 1 crv.24
    ic_hex_set_edge_projection 267 502 0 1 crv.24
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.27
    ic_undo_group_begin 
    ic_hex_set_edge_projection 289 299 0 1 crv.27
    ic_hex_set_edge_projection 295 506 0 1 crv.27
    ic_hex_set_edge_projection 285 450 0 1 crv.27
    ic_hex_set_edge_projection 285 295 0 1 crv.27
    ic_undo_group_end
    \n""")



    


    fid.write("""
    ic_geo_rename_family BODY SOLID 0
    ic_geo_rename_family BODY SOLID 1
    ic_geo_rename_family FLUID FLUID_BLK 0
    ic_geo_rename_family FLUID FLUID_BLK 1
    ic_undo_group_begin 
    ic_geo_new_family BODY
    ic_boco_set_part_color BODY
    ic_geo_create_body {srf.00 srf.02 srf.03 srf.04 srf.05 srf.06 srf.07 srf.10 srf.11 srf.12 srf.23 srf.14 srf.15 srf.16 srf.17 srf.18 srf.21 srf.22 srf.24 srf.25 srf.26 srf.27 srf.28 srf.29 srf.30 srf.31 srf.32 srf.33 srf.34} {} BODY
    ic_undo_group_end 
    ic_geo_rename_family BODY FLUID 0
    ic_geo_rename_family BODY FLUID 1
    ic_undo_group_begin 
    ic_geo_new_family SOLID_BLK
    ic_boco_set_part_color SOLID_BLK
    ic_hex_mark_blocks unmark
    ic_hex_mark_blocks superblock 157
    ic_hex_mark_blocks superblock 158
    ic_hex_mark_blocks superblock 163
    ic_hex_mark_blocks superblock 166
    ic_hex_mark_blocks superblock 43
    ic_hex_mark_blocks superblock 44
    ic_hex_mark_blocks superblock 45
    ic_hex_mark_blocks superblock 46
    ic_hex_mark_blocks superblock 67
    ic_hex_mark_blocks superblock 68
    ic_hex_mark_blocks superblock 72
    ic_hex_mark_blocks superblock 74
    ic_hex_mark_blocks superblock 79
    ic_hex_mark_blocks superblock 80
    ic_hex_mark_blocks superblock 85
    ic_hex_mark_blocks superblock 88
    ic_hex_mark_blocks superblock 91
    ic_hex_mark_blocks superblock 92
    ic_hex_mark_blocks superblock 96
    ic_hex_mark_blocks superblock 99
    ic_hex_change_element_id SOLID_BLK
    ic_undo_group_end 
    \n""")

    

    fid.write("""
    ic_undo_group_begin 
    ic_geo_set_part surface srf.00 INLET_M 0
    ic_delete_empty_parts 
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_geo_set_part surface srf.02 INLET_S 0
    ic_delete_empty_parts 
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_geo_set_part surface srf.18 INLET_B 0
    ic_delete_empty_parts 
    ic_undo_group_end 
    ic_undo_group_begin 
    ic_geo_set_part surface {srf.21 srf.22 srf.23} OUTLET 0
    ic_delete_empty_parts 
    ic_undo_group_end 
    \n""")

    
    fid.write("""
    ic_hex_set_mesh 215 243 n 12 h1 0.000286778 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 187 327 n 10 h1 0.002 h2 0.000286778 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 131 159 n 8 h1 0.002 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked 
    ic_hex_set_mesh 233 345 n 20 h1 0.002 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 341 369 n 15 h1 0.002 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 373 261 n 50 h1 0.002 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 25 26 n 10 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 21 410 n 15 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 25 26 n 15 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 21 22 n 15 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 26 71 n 8 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 187 327 n 20 h1 0.002 h2 0.000286778 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 327 215 n 5 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 187 327 n 30 h1 0.002 h2 0.000286778 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 215 243 n 25 h1 0.000286778 h2 0.002 r1 1.3 r2 1.3 lmax 0.002 default unlocked
    ic_hex_set_mesh 243 355 n 30 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked 
    ic_hex_set_mesh 274 302 n 30 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 289 41 n 15 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 261 289 n 50 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 149 177 n 6 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 532 534 n 5 h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked
    ic_hex_set_mesh 87 103 n 14 h1 1e-005 h2 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 86 102 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 84 100 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 85 101 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 123 124 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 133 134 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 119 120 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 129 130 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 161 162 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 151 152 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 147 148 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 157 158 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 87 103 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 319 320 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 315 316 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 203 204 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 207 208 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 329 330 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 217 218 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 325 326 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_hex_set_mesh 213 214 n 14 h1rel 0.0131578947368 h2rel 0.0 r1 1.3 r2 2 lmax 0 default unlocked
    ic_undo_group_end
    \n""")







    N_y_bypass=30
    fid.write("ic_hex_set_mesh 351 532 n %d h1 0.002 h2 0.00160067 r1 1.3 r2 1.3 lmax 0.002 default unlocked\n" %(N_y_bypass))
    
    Nx_tot= 200
    N_x_mixpart_start = math.ceil((xb-Lmch)/Lmix*Nx_tot)
    if (xb-Lmch<delta_x_l_alpha_i):
        N_x_mixpart_start = math.ceil((delta_x_l_alpha_i+xb-Lmch)/Lmix*Nx_tot)

    N_x_mixpart_bypass = math.ceil(L/Lmix*Nx_tot*1.2)
    N_x_mixpart_end = math.ceil( (Lmix-(xb+L-Lmch) )/Lmix *Nx_tot)

    N_theta_tot=80
    N_theta_side_tot=N_theta_tot/4
    weightingBypass = 1.3
    N_theta_bypass = math.ceil(theta/90*N_theta_side_tot*weightingBypass)
    N_theta_neighborBlock = N_theta_side_tot - N_theta_bypass
    if N_theta_neighborBlock < 2:
        N_theta_neighborBlock =2



    zero=0


    d_inner = 1e-4











    # ##################################################################################################################   THETA ##################################################################################################################

    fid.write("ic_hex_set_mesh 351 494 n %d h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %(N_theta_bypass))   
    fid.write("ic_hex_set_mesh 355 494 n %d h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %(N_theta_neighborBlock))  


    # theta dir toward bypass - theta


    # fid.write("ic_hex_set_mesh 379 498 n 10 h1 0.0002 h2 0.0 r1 1.3 r2 1.3 lmax 0 default unlocked\n")
    # fid.write("ic_hex_set_mesh 341 351 n 20 h1 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n")
    # fid.write("ic_hex_set_mesh 369 379 n 20 h1 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n")


    # ##################################################################################################################  R ##################################################################################################################


    # ###################### SUCTION ######################
    
    N_r_suc = 26
    d_s_start = 2.0e-5
    d_s_end = 1e-5
    d_s_mix = 1e-5
    d_s_mix_start_wall = 5e-6
    d_s_mix_end_wall = 5e-6
    d_s_tobypass = 5e-6
    d_s_diff = d_inner
    d_s_diff_wall_start = 2e-5


    # Suction section start - r
    fid.write("ic_undo_group_begin\n")  
    fid.write("ic_hex_set_mesh 22 69 n %d h1 %f h2 %f r1 1.3 r2 2 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 21 68 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 25 70 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 26 71 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 131 132 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 159 160 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 187 188 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 127 128 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 155 156 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 183 184 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 117 118 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 145 146 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 173 174 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 177 178 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 149 150 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_hex_set_mesh 121 122 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_start,d_s_start) )  
    fid.write("ic_undo_group_end\n")  


    # Suction section end -r

    fid.write("ic_undo_group_begin\n")
    fid.write("ic_hex_set_mesh 327 328 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 323 324 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 313 314 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 201 202 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 205 206 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 211 212 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 317 318 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_hex_set_mesh 215 216 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_end,d_s_end) )  
    fid.write("ic_undo_group_end\n")


    # outer suction in mixer -r

    fid.write("ic_hex_set_mesh 243 244 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_start_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 239 240 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_start_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 233 234 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_start_wall,d_s_mix) ) 
    fid.write("ic_hex_set_mesh 229 230 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_start_wall,d_s_mix) )   
    fid.write("ic_hex_set_mesh 271 272 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 261 262 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 267 268 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 257 258 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 323 324 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 313 314 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 201 202 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 205 206 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 211 212 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 317 318 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 355 356 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 341 342 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 345 346 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 369 370 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 383 384 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 373 374 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  

    # outer suction to bypass -r

    fid.write("ic_hex_set_mesh 494 495 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_tobypass,d_s_mix) )  
    fid.write("ic_hex_set_mesh 379 380 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_tobypass,d_s_mix) )  
    fid.write("ic_hex_set_mesh 351 352 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_tobypass,d_s_mix) )  
    fid.write("ic_hex_set_mesh 498 499 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_tobypass,d_s_mix) )  


    # outer suction in diffuser -r

    fid.write("ic_hex_set_mesh 42 75 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 41 74 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 37 72 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 38 73 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 289 290 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 299 300 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 295 296 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  
    fid.write("ic_hex_set_mesh 285 286 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_diff_wall_start,d_s_diff) )  





    # ###################### LIP ######################

    # motive lip -r

    N_lip = 25
    d_lip_inner_motive_outer =1e-6
    d_lip_outer= 2e-5

    fid.write("ic_hex_set_mesh 216 217 n %d h1 %f  h2 %f  r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_lip_outer,d_lip_inner_motive_outer) )   
    fid.write("ic_hex_set_mesh 212 213 n %d h1 %f  h2 %f  r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_lip_outer,d_lip_inner_motive_outer) )  
    fid.write("ic_hex_set_mesh 206 207 n %d h1 %f  h2 %f  r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_lip_outer,d_lip_inner_motive_outer) )  
    fid.write("ic_hex_set_mesh 202 203 n %d h1 %f  h2 %f  r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_lip_outer,d_lip_inner_motive_outer) )  

    # Lip and mixing layer -r
    fid.write("ic_hex_set_mesh 240 241 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 234 235 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 230 231 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 352 353 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 346 347 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 356 357 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 342 343 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 384 385 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 374 375 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 380 381 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 370 371 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 272 273 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 262 263 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 268 269 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 258 259 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  

    # Diffuser central -r

    fid.write("ic_hex_set_mesh 75 91 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 74 90 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 73 89 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  
    fid.write("ic_hex_set_mesh 72 88 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_lip,d_inner,d_inner) )  




    
    # ###################### INNER MOTIVE R ######################

    # central motive mixing chamber blocks -r
    fid.write("ic_hex_set_mesh 273 274 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 263 264 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 259 260 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 269 270 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 385 386 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 375 376 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 381 382 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 371 372 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 347 348 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 357 358 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 353 354 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 235 236 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 231 232 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 241 242 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  

    # central motive diffuser chamber blocks -r
    fid.write("ic_hex_set_mesh 297 298 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 291 292 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 287 288 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 91 107 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 90 106 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 88 104 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  
    fid.write("ic_hex_set_mesh 89 105 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )   
    fid.write("ic_hex_set_mesh 301 302 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_theta_side_tot,zero,zero) )  

    
    # ###################### OUTER MOTIVE R ######################

    # Motive outer start - r

    N_outer_mot_r = 16
    d_motive_wall_start = 2e-5

    fid.write("ic_hex_set_mesh 85 101 n %d h1 %f  h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 86 102 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 84 100 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 87 103 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 119 120 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 147 148 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 175 176 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 129 130 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 157 158 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 185 186 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 133 134 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 161 162 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 189 190 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 179 180 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 151 152 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  
    fid.write("ic_hex_set_mesh 123 124 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_wall_start,d_inner) )  

        # Motive nozzle -r

    d_motive_inner_wall_end = d_lip_inner_motive_outer
    d_motive_inner_centercore_end =2e-5

    fid.write("ic_hex_set_mesh 329 330 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 325 326 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 319 320 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 315 316 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 203 204 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 213 214 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 207 208 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  
    fid.write("ic_hex_set_mesh 217 218 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_outer_mot_r,d_motive_inner_wall_end,d_motive_inner_centercore_end) )  




    # bypass inlet -r

    N_r_bypass = 19
    # Inner wall wall - improve quality
    fid.write("ic_hex_set_mesh 494 534 n %d h1 0.0002 h2 0.0 r1 1.5 r2 2 lmax 0 default unlocked\n" %N_r_bypass)
    fid.write("ic_hex_set_mesh 498 538 n %d h1 0.0002 h2 0.0 r1 1.5 r2 2 lmax 0 default unlocked\n" %N_r_bypass)

    # Outer wall - improve quality
    fid.write("ic_hex_set_mesh 351 532 n %d h1 0.0007 h2 0.0 r1 3 r2 1.3 lmax 0 geo1 unlocked\n" %N_r_bypass)
    fid.write("ic_hex_set_mesh 379 536 n %d h1 0.0007 h2 0.0 r1 3 r2 1.3 lmax 0 geo1 unlocked\n" %N_r_bypass)

    # Along theta wall - improve quality

    fid.write("ic_hex_set_mesh 379 498 n %d h1 0.0 h2 0.0001 r1 2 r2 2 lmax 0 hcosinus2 unlocked\n" %N_theta_bypass)
    fid.write("ic_hex_set_mesh 351 494 n %d h1 0.0 h2 0.0001 r1 2 r2 2 lmax 0 hcosinus2 unlocked\n" %N_theta_bypass)

    N_x_bypass = 21
    # Inner wall wall - improve quality
    fid.write("ic_hex_set_mesh 532 536 n %d h1 0.000 h2 0.00005 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_bypass)
    fid.write("ic_hex_set_mesh 534 538 n %d h1 0.000 h2 0.00005 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_bypass)



    # ##################################################################################################################  X ##################################################################################################################

    #  x dir - diffuser
    N_x_diff= 260
    fid.write("ic_hex_set_mesh 271 299 n %d h1 0.0002 h2 0.0 r1 1.3 r2 1.3 lmax 0 default copy_to_parallel unlocked \n" %(N_x_diff)) 




    # x dir - mch outer

    N_x_mch = 60

    fid.write("ic_hex_set_mesh 205 233 n %d h1 0.0002 h2 0.000 r1 1.3 r2 1.3 lmax 0.002 default unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 215 243 n %d h1 0.0002 h2 0.0 r1 1.3 r2 1.3 lmax 0.002 default unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 211 239 n %d h1 0.0002 h2 0.0 r1 1.3 r2 1.3 lmax 0.002 default unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 201 229 n %d h1 0.0002 h2 0.0 r1 1.3 r2 1.3 lmax 0.002 default unlocked\n" %N_x_mch)

    # x dir - mch to wall (or motive)

    fid.write("ic_hex_set_mesh 216 244 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 217 245 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 206 234 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 207 235 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 213 241 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 212 240 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 203 231 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 202 230 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 218 246 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 208 236 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 214 242 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)
    fid.write("ic_hex_set_mesh 204 232 n %d h1 0.00003 h2 0.0 r1 1.3 r2 2 lmax 0.002 geo1 unlocked\n" %N_x_mch)




    
    #  x dir - mixing chamber
    
    delta_mix_start = Lmch/N_x_mch
    fid.write("ic_hex_set_mesh 244 356 n %d h1 %f h2 0.0 r1 1.3 r2 1.3 lmax 0 default copy_to_parallel\n" %(N_x_mixpart_start,delta_mix_start))
    fid.write("ic_hex_set_mesh 383 271 n %d h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0 default copy_to_parallel \n" %(N_x_mixpart_end))
    fid.write("ic_hex_set_mesh 341 369 n %d h1 0.0 h2 0.0 r1 1.3 r2 1.3 lmax 0.002 default copy_to_parallel\n" %(N_x_mixpart_bypass))

    N_x_inlet = 26  
    N_x_inlet_part2 = 3 
    fid.write("ic_hex_set_mesh 26 131 n %d h1rel 0.0 h2rel 0.0 r1 2 r2 2 lmax 0 default copy_to_parallel unlocked\n" %(N_x_inlet))
    fid.write("ic_hex_set_mesh 131 159 n %d h1rel 0.0 h2rel 0.0 r1 1.3 r2 1.3 lmax 0.002 default copy_to_parallel unlocked\n" %(N_x_inlet_part2))



    N_x_suction_to_mch =11
    # x dir - suction to mch outer

    fid.write("ic_hex_set_mesh 327 215 n %d h1rel 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 317 205 n %d h1rel 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 323 211 n %d h1rel 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 313 201 n %d h1rel 0.0 h2 0.0002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
 

    # x dir - suction to mch inner

    fid.write("ic_hex_set_mesh 328 216 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 318 206 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 324 212 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 314 202 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)

    #  x dir - nozzle out

    fid.write("ic_hex_set_mesh 330 218 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 320 208 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 326 214 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 316 204 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 319 207 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 325 213 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 329 217 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)
    fid.write("ic_hex_set_mesh 315 203 n %d h1rel 0.0 h2 0.00002 r1 1.3 r2 1.3 lmax 0 default unlocked\n" %N_x_suction_to_mch)






    # x dir - suction to suction end
    N_x_convergingnozzle = 30

    fid.write("ic_hex_set_mesh 177 317 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 187 327 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 183 323 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 173 313 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 184 324 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 178 318 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 174 314 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 188 328 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)



    # x dir - motive nozzle conv

    fid.write("ic_hex_set_mesh 175 315 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 179 319 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 185 325 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 190 330 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 180 320 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 186 326 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 176 316 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)
    fid.write("ic_hex_set_mesh 189 329 n %d h1rel 0.0 h2 0.0001 r1 2 r2 1.3 lmax 0 default unlocked\n" %N_x_convergingnozzle)






    #Smoothin cut... before bypass



    fid.write("""
    ic_undo_group_begin 
    ic_hex_mark_blocks unmark
    ic_hex_mark_blocks superblock 115
    ic_hex_mark_blocks superblock 116
    ic_hex_mark_blocks superblock 117
    ic_hex_mark_blocks superblock 118
    ic_hex_mark_blocks superblock 119
    ic_hex_mark_blocks superblock 120
    ic_hex_mark_blocks superblock 121
    ic_hex_mark_blocks superblock 122
    ic_hex_mark_blocks superblock 123
    ic_hex_mark_blocks superblock 124
    ic_hex_mark_blocks superblock 125
    ic_hex_mark_blocks superblock 126
    ic_hex_mark_blocks superblock 127
    ic_hex_undo_major_start split_grid
    \n""")  
    
    fid.write("ic_hex_split_grid 243 355 abs:%f m GEOM FLUID_BLK FLUID SHELL LUMP SOLID_BLK INLET_M INLET_S INLET_B OUTLET marked    \n" %((xb-Lmch)*0.9)  )

    fid.write("""
    ic_hex_undo_major_end split_grid
    ic_undo_group_end 
    \n""")  

    fid.write("ic_hex_get_node_location {  677  } _tempx _tempy _tempz \n")  
    fid.write("ic_hex_set_node_location x {$_tempx} -csys global node_numbers {{  667  }}    \n")  


    L_pre_bypass_block = (xb-Lmch)*0.1
    delta_pre_bypass_block=1e-4
    N_pre_bypass_block = L_pre_bypass_block/delta_pre_bypass_block

    fid.write("ic_hex_set_mesh 667 668 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 677 678 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )  
    fid.write("ic_hex_set_mesh 652 653 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) ) 
    fid.write("ic_hex_set_mesh 662 663 n %d h1 %f h2 %f r1 1.3 r2 1.3 lmax 0 default unlocked\n"%(N_r_suc,d_s_mix_end_wall,d_s_mix) )   
    fid.write("ic_hex_set_mesh 662 345 n %d h1rel 0.0 h2rel 0.0 r1 1.3 r2 1.3 lmax 0 default copy_to_parallel unlocked\n" %N_pre_bypass_block)
    

    fid.write("""
    ic_undo_group_begin 
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.23
    ic_undo_group_begin 
    ic_hex_set_edge_projection 262 272 0 1 crv.23
    ic_hex_set_edge_projection 268 503 0 1 crv.23
    ic_hex_set_edge_projection 258 268 0 1 crv.23
    ic_hex_set_edge_projection 258 447 0 1 crv.23
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.22
    ic_undo_group_begin 
    ic_hex_set_edge_projection 263 273 0 1 crv.22
    ic_hex_set_edge_projection 269 504 0 1 crv.22
    ic_hex_set_edge_projection 259 269 0 1 crv.22
    ic_hex_set_edge_projection 259 448 0 1 crv.22
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.26
    ic_undo_group_begin 
    ic_hex_set_edge_projection 290 300 0 1 crv.26
    ic_hex_set_edge_projection 296 507 0 1 crv.26
    ic_hex_set_edge_projection 286 451 0 1 crv.26
    ic_hex_set_edge_projection 286 296 0 1 crv.26
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.25
    ic_undo_group_begin 
    ic_hex_set_edge_projection 291 301 0 1 crv.25
    ic_hex_set_edge_projection 297 508 0 1 crv.25
    ic_hex_set_edge_projection 287 297 0 1 crv.25
    ic_hex_set_edge_projection 287 452 0 1 crv.25
    ic_undo_group_end 
    ic_hex_find_comp_curve srf.21e147
    ic_undo_group_begin 
    ic_hex_set_edge_projection 74 75 0 1 srf.21e147
    ic_hex_set_edge_projection 73 511 0 1 srf.21e147
    ic_hex_set_edge_projection 72 455 0 1 srf.21e147
    ic_hex_set_edge_projection 72 73 0 1 srf.21e147
    ic_undo_group_end 
    ic_hex_find_comp_curve crv.28
    ic_undo_group_begin 
    ic_hex_set_edge_projection 90 91 0 1 crv.28
    ic_hex_set_edge_projection 89 512 0 1 crv.28
    ic_hex_set_edge_projection 88 89 0 1 crv.28
    ic_hex_set_edge_projection 88 456 0 1 crv.28
    ic_undo_group_end \n""")


    
    return
    fid.write("""
    ic_hex_create_mesh GEOM FLUID_BLK FLUID SOLID_BLK INLET_M INLET_S INLET_B OUTLET proj 2 dim_to_mesh 3 nproc 20
    ic_hex_write_file ./hex.uns GEOM FLUID_BLK FLUID SHELL LUMP SOLID_BLK INLET_M INLET_S INLET_B OUTLET proj 2 dim_to_mesh 3 no_boco
    ic_uns_load ./hex.uns 3 0 {} 1
    ic_uns_update_family_type visible {LUMP GEOM SHELL INLET_M SOLID_BLK FLUID FLUID_BLK OUTLET ORFN INLET_B INLET_S} {!NODE !LINE_2 QUAD_4 !HEXA_8} update 0
    ic_boco_solver 
    ic_boco_clear_icons 
    ic_uns_update_family_type visible {LUMP GEOM SHELL INLET_M SOLID_BLK FLUID FLUID_BLK OUTLET ORFN INLET_B INLET_S} {!NODE !LINE_2 !QUAD_4 !HEXA_8} update 0
    ic_undo_group_begin 
    ic_uns_create_diagnostic_subset 
    ic_uns_create_diagnostic_edgelist 1
    ic_uns_diagnostic subset all diag_type vol_orient fix_fam FIX_VOL_ORIENT disp_subset uns_diag_0 diag_verb {Volume orientations} fams {} busy_off 1 skip 1
    ic_uns_fix_volume_orientation uns_diag_0
    ic_uns_subset_delete uns_diag_0
    ic_uns_create_diagnostic_edgelist 0
    ic_undo_group_end 
    \n""")

    filename="Mesh_3D_y+30.msh"
    fid.write("""
    ic_boco_solver 
    ic_boco_solver {ANSYS Fluent}
    ic_solution_set_solver {ANSYS Fluent} 1
    ic_boco_save C:/Users/knuterin/Documents/PhD/Fluent/Python/3Dejector/ansys.fbc
    ic_boco_save_atr C:/Users/knuterin/Documents/PhD/Fluent/Python/3Dejector/ansys.atr
    ic_delete_empty_parts 
    ic_delete_empty_parts 
    ic_save_tetin project2.tin 0 0 {} {} 0 0 1
    ic_uns_check_duplicate_numbers 
    ic_save_unstruct project2.uns 1 {} {} {}
    ic_uns_set_modified 1
    ic_hex_save_blocking project2.blk
    ic_boco_solver 
    ic_boco_solver {ANSYS Fluent}
    ic_solution_set_solver {ANSYS Fluent} 1
    ic_boco_save project2.fbc
    ic_boco_save_atr project2.atr
    ic_save_project_file C:/Users/knuterin/Documents/PhD/Fluent/Python/3Dejector/project2.prj {array\ set\ file_name\ \{ {    catia_dir .} {    parts_dir .} {    domain_loaded 0} {    cart_file_loaded 0} {    cart_file {}} {    domain_saved project2.uns} {    archive {}} {    med_replay {}} {    topology_dir .} {    ugparts_dir .} {    icons {{$env(ICEM_ACN)/lib/ai_env/icons} {$env(ICEM_ACN)/lib/va/EZCAD/icons} {$env(ICEM_ACN)/lib/icons} {$env(ICEM_ACN)/lib/va/CABIN/icons}}} {    tetin project2.tin} {    family_boco project2.fbc} {    iges_dir .} {    solver_params_loaded 0} {    attributes_loaded 0} {    project_lock {}} {    attributes project2.atr} {    domain project2.uns} {    domains_dir .} {    settings_loaded 0} {    settings project2.prj} {    blocking project2.blk} {    hexa_replay {}} {    transfer_dir .} {    mesh_dir .} {    family_topo {}} {    gemsparts_dir .} {    family_boco_loaded 0} {    tetin_loaded 0} {    project_dir .} {    topo_mulcad_out {}} {    solver_params {}} \} array\ set\ options\ \{ {    expert 1} {    remote_path {}} {    tree_disp_quad 2} {    tree_disp_pyra 0} {    evaluate_diagnostic 0} {    histo_show_default 1} {    select_toggle_corners 0} {    remove_all 0} {    keep_existing_file_names 0} {    record_journal 0} {    edit_wait 0} {    face_mode all} {    select_mode all} {    med_save_emergency_tetin 1} {    user_name knuterin} {    diag_which all} {    uns_warn_if_display 500000} {    bubble_delay 1000} {    external_num 1} {    tree_disp_tri 2} {    apply_all 0} {    temporary_directory C:/Users/knuterin/AppData/Local/Temp} {    flood_select_angle 0} {    home_after_load 1} {    project_active 0} {    histo_color_by_quality_default 1} {    undo_logging 1} {    tree_disp_hexa 0} {    histo_solid_default 1} {    host_name NTNU07708} {    xhidden_full 1} {    replay_internal_editor 1} {    editor notepad} {    mouse_color orange} {    clear_undo 1} {    remote_acn {}} {    remote_sh csh} {    tree_disp_penta 0} {    n_processors 20} {    remote_host {}} {    save_to_new 0} {    quality_info Quality} {    tree_disp_node 0} {    med_save_emergency_mesh 1} {    redtext_color red} {    tree_disp_line 0} {    select_edge_mode 0} {    use_dlremote 0} {    max_mesh_map_size 1024} {    show_tris 1} {    remote_user {}} {    enable_idle 0} {    auto_save_views 1} {    max_cad_map_size 512} {    display_origin 0} {    uns_warn_user_if_display 1000000} {    detail_info 0} {    win_java_help 0} {    show_factor 1} {    boundary_mode all} {    clean_up_tmp_files 1} {    auto_fix_uncovered_faces 1} {    med_save_emergency_blocking 1} {    max_binary_tetin 0} {    tree_disp_tetra 0} \} array\ set\ disp_options\ \{ {    uns_dualmesh 0} {    uns_warn_if_display 500000} {    uns_normals_colored 0} {    uns_icons 0} {    uns_locked_elements 0} {    uns_shrink_npos 0} {    uns_node_type None} {    uns_icons_normals_vol 0} {    uns_bcfield 0} {    backup Wire} {    uns_nodes 0} {    uns_only_edges 0} {    uns_surf_bounds 0} {    uns_wide_lines 0} {    uns_vol_bounds 0} {    uns_displ_orient Triad} {    uns_orientation 0} {    uns_directions 0} {    uns_thickness 0} {    uns_shell_diagnostic 0} {    uns_normals 0} {    uns_couplings 0} {    uns_periodicity 0} {    uns_single_surfaces 0} {    uns_midside_nodes 1} {    uns_shrink 100} {    uns_multiple_surfaces 0} {    uns_no_inner 0} {    uns_enums 0} {    uns_disp Wire} {    uns_bcfield_name {}} {    uns_color_by_quality 0} {    uns_changes 0} {    uns_cut_delay_count 1000} \} {set icon_size1 24} {set icon_size2 35} {set thickness_defined 0} {set solver_type 1} {set solver_setup -1} array\ set\ prism_values\ \{ {    n_triangle_smoothing_steps 5} {    min_smoothing_steps 6} {    first_layer_smoothing_steps 1} {    new_volume {}} {    height {}} {    prism_height_limit {}} {    interpolate_heights 0} {    n_tetra_smoothing_steps 10} {    do_checks {}} {    delete_standalone 1} {    ortho_weight 0.50} {    max_aspect_ratio {}} {    ratio_max {}} {    incremental_write 0} {    total_height {}} {    use_prism_v10 0} {    intermediate_write 1} {    delete_base_triangles {}} {    ratio_multiplier {}} {    verbosity_level 1} {    refine_prism_boundary 1} {    max_size_ratio {}} {    triangle_quality {}} {    max_prism_angle 180} {    tetra_smooth_limit 0.3} {    max_jump_factor 5} {    use_existing_quad_layers 0} {    layers 3} {    fillet 0.10} {    into_orphan 0} {    init_dir_from_prev {}} {    blayer_2d 0} {    do_not_allow_sticking {}} {    top_family {}} {    law exponential} {    min_smoothing_val 0.1} {    auto_reduction 0} {    stop_columns 1} {    stair_step 1} {    smoothing_steps 12} {    side_family {}} {    min_prism_quality 0.01} {    ratio 1.2} \} {set aie_current_flavor {}} array\ set\ vid_options\ \{ {    wb_import_mat_points 0} {    wb_NS_to_subset 0} {    wb_import_cad_att_pre {SDFEA;DDM}} {    wb_import_tritol 0.001} {    wb_import_mix_res -1} {    wb_import_save_pmdb {}} {    composite_tolerance 1.0} {    wb_import_save_partfile 0} {    wb_NS_to_entity_parts 0} {    wb_import_reference_key 0} {    replace 0} {    tdv_axes 1} {    vid_mode 0} {    auxiliary 0} {    wb_import_surface_bodies 1} {    show_name 0} {    wb_import_cad_att_trans 1} {    wb_import_solid_bodies 1} {    default_part GEOM} {    wb_import_mix_res_solid 0} {    new_srf_topo 1} {    DelPerFlag 0} {    wb_import_associativity_model_name {}} {    wb_import_work_points 0} {    wb_import_sel_proc 1} {    wb_NS_only 0} {    wb_import_scale_geo Default} {    wb_import_lcs 0} {    same_pnt_tol 1e-4} {    wb_import_transfer_file_scale 1.0} {    DelBlkPerFlag 0} {    wb_import_mesh 0} {    wb_import_mix_res_surface 0} {    wb_import_analysis_type 3} {    wb_import_geom 1} {    wb_import_refresh_pmdb 0} {    wb_import_load_pmdb {}} {    wb_import_mix_res_line 0} {    wb_import_delete_solids 0} {    inherit 1} {    wb_import_line_bodies 0} {    wb_import_en_sym_proc 1} {    wb_run_mesher tetra} {    wb_import_mix_res_point 0} {    wb_import_pluginname {}} {    wb_import_create_solids 0} {    wb_import_sel_pre {}} {    wb_import_cad_associativity 0} \} {set savedTreeVisibility {geomNode 1 geom_subsetNode 2 geomPointNode 0 geomCurveNode 2 geomSurfNode 2 geomBodyNode 2 meshNode 0 mesh_subsetNode 0 meshPointNode 0 meshLineNode 0 meshShellNode 0 meshQuadNode 0 meshVolumeNode 0 meshHexaNode 0 blockingNode 1 block_subsetNode 2 block_vertNode 0 block_edgeNode 2 block_faceNode 0 block_blockNode 0 block_meshNode 0 topoNode 2 topo-root 2 partNode 2 part-FLUID 2 part-FLUID_BLK 2 part-GEOM 2 part-INLET_B 2 part-INLET_M 2 part-INLET_S 2 part-OUTLET 2 part-SOLID_BLK 2 part-VORFN 0}} {set last_view {rot {-0.145617070209 0.371810537523 0.0251288596233 0.916472112746} scale {31621.1602507 31621.1602507 31621.1602507} center {0.075413968701499995 -4.6566128730800003e-10 -4.6566128730800003e-10} pos {-1156.02650724 94.7608725703 0}}} array\ set\ cut_info\ \{ {    a 0} {    b 0} {    xyz {0.075413968701468548 -4.6566128730773926e-10 -4.6566128730773926e-10}} {    c 1} {    active 0} {    d -4.6566128730773926e-10} {    pt1 {-0.00010205002763541415 -0.0095430510118603706 -4.6566128730773926e-10}} {    nx 0} {    pt2 {-0.00010205002763541415 0.009543050080537796 -4.6566128730773926e-10}} {    ny 0} {    pt3 {0.15092998743057251 -0.0095430510118603706 -4.6566128730773926e-10}} {    nz 1} {    whole 1} \} array\ set\ hex_option\ \{ {    default_bunching_ratio 2.0} {    floating_grid 0} {    project_to_topo 0} {    n_tetra_smoothing_steps 20} {    sketching_mode 0} {    trfDeg 1} {    wr_hexa7 0} {    smooth_ogrid 0} {    find_worst 1-3} {    hexa_verbose_mode 0} {    old_eparams 0} {    uns_face_mesh_method uniform_quad} {    multigrid_level 0} {    uns_face_mesh one_tri} {    check_blck 0} {    proj_limit 0} {    check_inv 0} {    project_bspline 0} {    hexa_update_mode 1} {    default_bunching_law BiGeometric} {    worse_criterion Quality} \} array\ set\ saved_views\ \{ {    views {}} \}} {ICEM CFD}
    ic_exec {C:/Program Files/ANSYS Inc/v192/icemcfd/win64_amd/icemcfd/output-interfaces/fluent6} -dom C:/Users/knuterin/Documents/PhD/Fluent/Python/3Dejector/project2.uns -b project2.fbc ./%s
    ic_uns_num_couplings 
    ic_undo_group_begin 
    ic_uns_create_diagnostic_edgelist 1
    ic_uns_diagnostic subset all diag_type uncovered fix_fam FIX_UNCOVERED diag_verb {Uncovered faces} fams {} busy_off 1 quiet 1
    ic_uns_create_diagnostic_edgelist 0
    ic_undo_group_end 
    ic_uns_min_metric Quality {} {}
    \n"""%filename)



    fid.close()









    if write_mesh == True:
        command ='"C:/Program Files/ANSYS Inc/v192/icemcfd/win64_amd/bin/icemcfd.bat" -batch -script ./EjectorBypassSwirl_meshing_3D.rpl & '
        subprocess.call(command)

    








meshing3D_test()
