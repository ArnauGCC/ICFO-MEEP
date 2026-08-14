import meep as mp
import math
import dataclasses


@dataclasses.dataclass
class Layer:
    width: float
    index:  float


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


def create_h_grating(gr_period, gr_height, gr_duty_cycle, n_cells, wg_width, n,
                     center=mp.Vector3(), gr_up=True, gr_down=False, n_ext=1, security_factor=1.1):
    """
    Creates an (infinite) horitzontal waveguide with a grating.
    gr_period:      Period of the grating
    gr_height:      How deep are the holes
    gr_duty_cycle:  
    ...
    security_factor:  [>= 1] This parameter doesn't affect to the hole deep- Factor to make the grating holes bigger than the grating height. This is to avoid that the holes are closed due to the resolution of the simulation.
    """

    wg_y = center.y
    geometry = [create_h_waveguide(wg_y, wg_width, n)]


    length = gr_period * n_cells
    left = -length/2 + gr_period*(1-gr_duty_cycle)/2 + center.x
    gdc = 1-gr_duty_cycle

    if gr_up:
        gy = wg_y + wg_width/2  + (security_factor/2 - 1)*gr_height

        for x in range(n_cells):
            geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*(1-gdc)), gy),
                                    size=mp.Vector3(gr_period*gdc, gr_height*security_factor),
                                    material=mp.Medium(index=n_ext)
                                    )
                            )

    if gr_down: 
        gy = wg_y - wg_width/2 + (1 - security_factor/2)*gr_height

        for x in range(n_cells):
            geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*(1-gdc)), gy),
                                    size=mp.Vector3(gr_period*gdc, gr_height*security_factor),
                                    material=mp.Medium(index=n_ext)
                                    )
                            )
    return geometry


def create_v_grating(gr_period, gr_height, gr_duty_cycle, n_cells, wg_width, n,
                     center=mp.Vector3(), gr_right=True, gr_left=False, n_ext=1):

    wg_x = center.x
    geometry = [create_v_waveguide(wg_x, wg_width, n)]


    length = gr_period * n_cells
    up = -length/2 + gr_period*(1-gr_duty_cycle)/2 + center.y
    gdc = 1-gr_duty_cycle

    if gr_right:
        gx = wg_x + (wg_width - gr_height)/2

        for y in range(n_cells):
            geometry.append(mp.Block(center=mp.Vector3(gx, up + gr_period*(y + 0.5*(1-gdc))),
                                    size=mp.Vector3(gr_height, gr_period*gdc),
                                    material=mp.Medium(index=n_ext)
                                    )
                            )

    if gr_left:
        gx = wg_x - (wg_width - gr_height)/2
    
        for y in range(n_cells):
            geometry.append(mp.Block(center=mp.Vector3(gx, up + gr_period*(y + 0.5*(1-gdc))),
                                    size=mp.Vector3(gr_height, gr_period*gdc),
                                    material=mp.Medium(index=n_ext)
                                    )
                            )
    return geometry