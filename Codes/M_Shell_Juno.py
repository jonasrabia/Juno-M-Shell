import numpy as np
from JRM33 import GetJupiterMag 


# ------------------------------------------------------------------
#   Compute the M-Shell of Juno
# 
#   Required input : 
#       - r | Radial distance in Rj (1 Rj = 71,492 km)
#       - t | Colatitude S3RH in deg 
#       - p | Longitude S3RH in deg 
#
#   Output : 
#       - M | /   | M-Shell, radial distance of the intersection between the magnetic field line and the magnetic equator, normalized to 1 Rj
#       - t | deg | Colatitude (S3RH) of the intersection between the magnetic field line and the magnetic equator
#       - p | deg | Longitude (S3RH) of the intersection between the magnetic field line and the magnetic equator 	
#       - S |  Rj | Curvilinear distance along the magnetic field line between Juno and the magnetic equator
#
#   Requirements : 
#       - JRM33 module 
#       - CON2020 module 
# 
#   Written by Jonas RABIA (IRAP-CNRS), contact : jonas.rabia@irap.omp.eu 
#-------------------------------------------------------------------



# -------------- Coordinates conversion -----------------------------

def Sph2CarP(rp, thetap, phip):
    # ==================================================
    # Spherical to cartesian coordinates for positions.
    # ==================================================
    # Function:
    #   input: r, theta, phi (arrays) in spherical coords.
    #   output: x, y, z

    r = np.array(rp)
    theta = np.array(thetap)*np.pi/180
    phi = np.array(phip)*np.pi/180

    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)    
    return np.transpose(np.array([x, y, z]))

def Car2SphP(xp, yp, zp):
    # ==================================================
    # Cartesian to Spherical coordinates for positions.
    # ==================================================
    # Function:
    #   input: x, y, z (arrays) in cartesian coords.
    #   output: r, theta, phi in degrees
    x = np.array(xp)
    y = np.array(yp)
    z = np.array(zp)
    n = len(x)
    r = np.maximum(np.sqrt(x**2 + y**2 + z**2), 1e-15*np.ones(n))

    theta = np.arccos(z/r)*180/np.pi
    phi = np.arctan2(y, x)*180/np.pi

    ir = np.where(phi < 0)[0]
    if len(ir) > 0:
        phi[ir] += 360.0

    return np.transpose(np.array([r, theta, phi]))


def magnetic_lat_CON2020(lon,lat):
    # Magnetic latitude wrt CON2020 coordinates 
    # See Connerney et al.2020
    # Input : longitude & latitude [°] in S3RH
    # Output : magnetic latitude [°]
    m_lat=9.3*np.cos(np.radians(lon-155.8))+lat
    return m_lat


# -------------- M-Shell calculation -----------------------------

def M_Shell(r,t,p):
    # =================================================
    # This program compute the Juno's M-Shell parameter.
    # 
    # Written by Jonas Rabia - 2023-11-08.
    # =================================================    
    # Get inital position 
    
    jno_r = [r]
    jno_t = [t]
    jno_p = [p]

    # Numerical parameter to rebuild the magnetic field line, has to be small enough   
    ds = 1/250  

    # Get Juno State
    jno_xyz = Sph2CarP(jno_r, jno_t, jno_p)

    # Get Juno distance to Jupiter
    r_in = jno_r

    # Get Magnetic field at Juno's position
    Bxyz = GetJupiterMag(jno_r, jno_t, jno_p, xyz_in=False, Bxyz_Out=True, CAN=True)
    Bx   = Bxyz[:,0]
    By   = Bxyz[:,1]
    Bz   = Bxyz[:,2]
    Bmag = Bxyz[:,3]

    

    # rebuilding the magnetic field line to get the M-Shell parameter of Juno,
    # the distance from the field line to Jupiter at the magnetic equator
    # we create a fictive particle which will follows to the magnetic field line
    fictive_particle_position_xyz = jno_xyz
    fictive_particle_position_rtp = Car2SphP([fictive_particle_position_xyz[0][0]], [fictive_particle_position_xyz[0][1]],
                                             [fictive_particle_position_xyz[0][2]])

    mlat = magnetic_lat_CON2020(fictive_particle_position_rtp[0][2], 90-fictive_particle_position_rtp[0][1])

    # direction : can be 1 or -1 and corresponds to the direction in which the "fictive particle" follows the magnetic field line
    # 1 : Juno is above the magnetic equator : we follow the field line to the south
    #-1 : Juno is below the magnetic equator : we follow the field line to the north
    if mlat > 0:
        direction = 1
    else:
        direction = -1
    
    
    # initialization of lists
    B = [Bmag,Bmag,Bmag]
    r_list = [r_in[0],r_in[0],r_in[0]]
    t_list = [jno_t[0],jno_t[0],jno_t[0]]
    p_list = [jno_p[0],jno_p[0],jno_p[0]]


    # Breaking conditions
    n_max = 20000
    r_max = 30
    zm_max = 5.8
    
    for i in range(2,n_max):

        # direction of the magnetic field
        unitaire = np.array([Bx[0], By[0], Bz[0]] / np.sqrt(Bx[0]**2 + By[0]**2 + Bz[0]**2))

        # make the fictive particle follows the magnetic line, step by step
        fictive_particle_position_xyz = fictive_particle_position_xyz + \
direction*np.array([ds*unitaire])
        fictive_particle_position_rtp = Car2SphP(
            fictive_particle_position_xyz[:, 0], fictive_particle_position_xyz[:, 1], fictive_particle_position_xyz[:, 2])

        # Get Magnetic field
        Bxyz = GetJupiterMag(fictive_particle_position_rtp[:, 0], fictive_particle_position_rtp[:, 1], fictive_particle_position_rtp[:, 2], xyz_in=False, Bxyz_Out=True, CAN=True)
        Bx = Bxyz[:,0]
        By = Bxyz[:,1]
        Bz = Bxyz[:,2]
        Bmag = Bxyz[:,3]

        # Get fictive particle distance to Jupiter
        r_in = fictive_particle_position_rtp[0][0]
        mlat = magnetic_lat_CON2020(fictive_particle_position_rtp[0][2], 90-fictive_particle_position_rtp[0][1])
        z_mag = r_in*np.cos(np.radians(90-mlat))

        # Saving iteration data
        B.append(Bmag)
        r_list.append(r_in)
        t_list.append(fictive_particle_position_rtp[0][1])
        p_list.append(fictive_particle_position_rtp[0][2])

        
        # Breaking conditions 
        if B[i]>B[i-1] and B[i-1]>B[i-2]:
            M = r_list[i-2]
            t_end = t_list[i-2]
            p_end = p_list[i-2]
            S = (i-4)*ds  # Loop start at i=2 and end at i=imax-2
            break 
        
        if r_in > r_max :
            M = np.nan 
            t_end = np.nan
            p_end = np.nan
            S = np.nan
            break 
        
        if np.abs(z_mag) > zm_max : 
            M = np.nan 
            t_end = np.nan
            p_end = np.nan
            S = np.nan
            break 
        
        
        if i == n_max - 1 : 
            M = np.nan
            t_end = np.nan
            p_end = np.nan
            S = np.nan
            break
        
    return M, t_end, p_end, S


