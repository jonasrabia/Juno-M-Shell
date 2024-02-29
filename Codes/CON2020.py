import numpy as np 

# ------------------------------------------------------------------
#   Compute Jupiter's external magnetic field with the CON2020 model 
#   See Connerney et al.,2020 (DOI:10.1029/2020JA028138) for further information on the model
#
#   Adapted by Jonas Rabia (IRAP-CNRS) from the CAN81 model,contact : jonas.rabia@irap.omp.eu
#   Azimuthal magnetic field calculation is based on the CON2020 community code 
#   See https://github.com/gabbyprovan/con2020 and Wilson et al., 2023 (DOI:10.1007/s11214-023-00961-3) for 
#   further information on the numerical models
#
#   Required input : 
#       - r | Radial distance in Rj (1 Rj = 71,492 km)
#       - t | Colatitude S3RH in deg 
#       - p | Longitude S3RH in deg 
#
#   Output : 
#       - Br, Bt, Bp [nT] in S3RH
#
#-------------------------------------------------------------------


# ---------------- Coordinates conversion ---------------------------
def ConvInputPol(r,theta,phi):
    '''
    Converts input coordinates from spherical polar right-handed 
    System III to Cartesian current sheet coordinates.
        
    Inputs
    ======
    r : float
        System III radial distance (Rj).
    theta : float
        System III colatitude (°).
    phi : float
        System III east longitude (°).
            
    Returns
    =======
    x1 : float
        x current sheet coordinate
    y1 : float
        y current sheet coordinate
    z1 : float
        z current sheet coordinate
    rho1 : float
        distance from z-axis (Rj).
    abs_z1 : float
        abs(z1) (Rj).
    cost : float
        cos(theta) - where theta is the colatitude
    sint : float
        sin(theta)
    cosp : float
        cos(phi) - where phi is east longitude
    sinp : float    
        sin(phi)
    '''        
    sint = np.sin(np.radians(theta))
    cost = np.cos(np.radians(theta))
    sinp = np.sin(np.radians(phi))
    cosp = np.cos(np.radians(phi))
    
    cosxp = np.cos(np.radians(155.8-180))
    sinxp = np.sin(np.radians(155.8-180))
    cosxt = np.cos(np.radians(9.3))
    sinxt = np.sin(np.radians(9.3))


    x = r*sint*(cosp*cosxp + sinp*sinxp)
    y1 = r*sint*(sinp*cosxp - cosp*sinxp)
    z = r*cost

    x1 = x*cosxt + z*sinxt
    z1 = z*cosxt - x*sinxt    
    
    #some other bits we need for the model
    rho1 = np.sqrt(x1*x1 + y1*y1)
    abs_z1 = np.abs(z1)
          
    # Conversion to list 
    """
    x1 = list(x1)
    y1 = list(y1)
    z1 = list(z1)
    rho1 = list(rho1)
    abs_z1 = list(abs_z1)
    cost = list(cost)
    sint = list(sint)
    cosp = list(cosp)
    sinp = list(sinp) 
    """
    
    return x1,y1,z1,rho1,abs_z1,cost,sint,cosp,sinp    


def ConvOutputPol(cost,sint,cosp,sinp,x1,y1,rho1,Brho1,Bphi1,Bz1):
        '''
        Convert the output magnetic field from cylindrical current 
        sheet coordinates to spherical polar right-handed System III
        coordinates.
        
        Inputs
        ======
        cost : float
            cos(theta) - where theta is the colatitude
        sint : float
            sin(theta)
        cosp : float
            cos(phi) - where phi is east longitude
        sinp : float    
            sin(phi)
        x1 : float
            x-position in current sheet coords (Rj).
        y1 : float
            y-position in current sheet coords (Rj).
        rho1 : float
            distance from z-axis (Rj).
        Brho1 : float    
            Rho component of magnetic field (nT).
        Bphi1 : float
            Phi (azimuthal) component of the magnetic field (nT).
        Bz1 : float
            z component of the magnetic field (nT).
            
        Returns
        =======
        Br : float
            Radial component of magnetic field in right-handed System 
            III coordinates (nT).
        Bt : float
            Meridional component of magnetic field in right-handed 
            System III coordinates (nT).
        Bp : float
            Azimuthal component of magnetic field in right-handed System 
            III coordinates (nT).
            
        
        '''        
        cosxp = np.cos(np.radians(155.8-180))
        sinxp = np.sin(np.radians(155.8-180))
        cosxt = np.cos(np.radians(9.3))
        sinxt = np.sin(np.radians(9.3))
      
        #this now runs in about 60% of the time it used to
        cosphi1 = x1/rho1
        sinphi1 = y1/rho1
        
        
        
        Bx1 = Brho1*cosphi1 - Bphi1*sinphi1
        By1 = Brho1*sinphi1 + Bphi1*cosphi1         

        Bx = Bx1*cosxt - Bz1*sinxt
        Bz = Bx1*sinxt + Bz1*cosxt        

        Bx2 = Bx*cosxp - By1*sinxp
        By2 = By1*cosxp + Bx*sinxp    

        Br =  Bx2*sint*cosp+By2*sint*sinp+Bz*cost
        Bt =  Bx2*cost*cosp+By2*cost*sinp-Bz*sint
        Bp = -Bx2*     sinp+By2*     cosp
    
        return Br,Bt,Bp


# ---------------- External magnetic field calculation  ---------------------------

    
def CON2020(r_in,theta_in, phi_in):
    # --------------------------------------------
    # Description:
    #   Jupiter's external magnetic field model
    #   DOI : 10.1029/2020JA028138
    #
    # Input:
    #   r, theta, phi ['Rj','deg','deg'] (theta = colat)
    #   Format = float or list
    #
    # Output:
    #   Br,Bt,Bp [nT] in S3RH coordinate system 
    #   Format = list 
    #
    # History:
    #   Adapted by Jonas Rabia from the CAN81 function 
    #   Contact : jonas.rabia@irap.omp.eu
    #
    # --------------------------------------------

    # Conversion from SIIIRH to current sheet reference frame 
    x,y,z,rho,abs_z,cost,sint,cosp,sinp = ConvInputPol(r_in, theta_in, phi_in)
    

    # VALUES CON2020
    D = 3.6             # Rj
    a = 7.8             # Rj 
    b = 51.4            # Rj
    mu0_I0 = 139.6*2    # nT
    i_rho=16.7

      # VALUES MAX VOGT ET AL. 2017
    # mu0_I0 = 494    # nT
    # i_rho=0.0   
    # D = 2.5          # Rj
    # a = 5.0           # Rj 
    # b = 50.0            # Rj

    
    # [Approximate formulas given in Connerney+1981_The magnetic field in Jupiter.]
    # [Cylindrical coordinates: in nT]
    if (rho <= a):
        F1 = np.sqrt((z-D)**2 + a**2)
        F2 = np.sqrt((z+D)**2 + a**2)
        F3 = np.sqrt(z**2 + a**2)

        Brho = 0.5*rho*(1.0/F1 - 1.0/F2)
        tmp = (z - D)/F1**3 - (z + D)/F2**3
        Bz = 2.0*D/F3 - 0.25*rho**2*tmp
    else:
        F1 = np.sqrt((z-D)**2 + rho**2)
        F2 = np.sqrt((z+D)**2 + rho**2)
        F3 = np.sqrt(z**2 + rho**2)

        Brho = (F1-F2+2*D)/rho
        if (abs(z) >= D) and (z < 0):
            Brho = (F1-F2-2*D)/rho
        if (abs(z) < D):
            Brho = (F1-F2+2*z)/rho
        Brho -= 0.25*a**2*rho*(1.0/F1**3 - 1.0/F2**3)
        tmp = (z-D)/F1**3 - (z+D)/F2**3
        Bz = 2.0*D/F3 - 0.25*a**2*tmp

    F1 = np.sqrt((z-D)**2 + b**2)
    F2 = np.sqrt((z+D)**2 + b**2)
    F3 = np.sqrt(z**2 + b**2)
    Brho2 = 0.5*rho*(1/F1 - 1/F2)
    Bz2 = 2*D/F3 - 0.25*rho**2*((z-D)/F1**3 - (z+D)/F2**3)

    Brho -= Brho2
    Bz -= Bz2

    Brho *= 0.5*mu0_I0
    Bz *= 0.5*mu0_I0


    # Bphi calculation 
    # Function taken from https://github.com/gabbyprovan/con2020/blob/master/con2020/Model.py 
    abs_z = np.abs(z)
    Bphi = 2.7975*i_rho/rho
    if np.size(rho) == 1:
        if abs_z < D:
            Bphi *= (abs_z/D)
        if z > 0:
            Bphi = -Bphi
    else:
        ind = np.where(abs_z < D)[0]
        if ind.size > 0:
            Bphi[ind] *= (abs_z[ind]/D)
        ind = np.where(z > 0)[0]
        if ind.size > 0:
            Bphi[ind] = -Bphi[ind]
            
    
    # Conversion from cylindrical to spherical field 
    Br,Bt,Bp = ConvOutputPol(cost, sint, cosp, sinp, x, y, rho, Brho, Bphi, Bz)


    return [Br,Bt,Bp]

