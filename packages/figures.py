import meep as mp
import math

def create_ideal_prism(alpha_deg, n, prism_length, sxy, offs=0):
    # Implementation for creating an ideal prism (with a triangular cross-section).
    
    alpha = math.radians(alpha_deg)         # triangle angle
    pad = (sxy-prism_length)/2

    base_prism = [mp.Vector3(-prism_length/2 + offs, 0, -sxy/2 + pad), 
                  mp.Vector3(sxy/2, 0, -sxy/2 + pad)]


    # en funció d'alpha crític (quan cau just a la cantonada)
    alpha_c = math.atan((sxy/2 - pad)/(prism_length - offs))

    if alpha < alpha_c:
            base_prism.append(mp.Vector3(sxy/2, (prism_length - offs + pad)*math.tan(alpha), -sxy/2 + pad))

    else:
            base_prism.append(mp.Vector3(sxy/2, sxy/2, -sxy/2 + pad))

            if alpha > alpha_c:
                    base_prism.append(mp.Vector3(sxy/(2*math.tan(alpha))-(prism_length - offs)/2, sxy/2, -sxy/2 + pad))


    return mp.Prism(base_prism, height=sxy-2*pad, 
                                       axis=mp.Vector3(0, 0, 1), 
                                       material=mp.Medium(index=n))


def create_prism(alpha_deg, n, base_length, offsx=0, offsy=0, height=0):
        alpha = math.radians(alpha_deg)
        h = base_length

        if height > 0:
                h = height

        base_prism = [mp.Vector3(base_length/2 + offsx, 0 + offsy, -h/2),
                      mp.Vector3(offsx, base_length*math.tan(alpha)/2 + offsy, -h/2),
                      mp.Vector3(-base_length/2 + offsx, offsy, -h/2)]

        return mp.Prism(base_prism, height=h, axis=mp.Vector3(z=1), material=mp.Medium(index=n))


def create_h_waveguide(center_y=0, width=1, n=1):
        return mp.Block(center=mp.Vector3(0, center_y), 
                        size=mp.Vector3(mp.inf, width), 
                        material=mp.Medium(index=n))

def create_v_waveguide(center_x, width, n):
        return mp.Block(center=mp.Vector3(center_x), 
                        size=mp.Vector3(width, mp.inf), 
                        material=mp.Medium(index=n))