import meep as mp
import math

def create_ideal_prism(alpha_deg, n, left_vertex:mp.Vector3, sx, sy):
        # Implementation for creating an ideal prism (with a triangular cross-section).

        alpha = math.radians(alpha_deg)         # triangle angle

        height = (sx/2 - left_vertex.x) * math.tan(alpha) + left_vertex.y
        base_prism = [mp.Vector3(sx/2, left_vertex.y, -sx/2), 
                      mp.Vector3(sx/2, height if height < sy/2 else sy/2, -sx/2)]

        if height > sy/2:
                base_prism.append(mp.Vector3(left_vertex.x + (sy/2 - left_vertex.y)/math.tan(alpha), sy/2, -sx/2))


        if left_vertex.x < -sx/2:
                base_prism.append(mp.Vector3(-sx/2, (-sx/2 - left_vertex.x )*math.tan(alpha) + left_vertex.y, -sx/2))
                base_prism.append(mp.Vector3(-sx/2, left_vertex.y, -sx/2))

        else:   base_prism.append(mp.Vector3(left_vertex.x, left_vertex.y, -sx/2))

        return mp.Prism(base_prism, height=sx, 
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