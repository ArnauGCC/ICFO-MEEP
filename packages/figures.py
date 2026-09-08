import meep as mp
import math
import dataclasses


@dataclasses.dataclass
class Layer:
    """
    A Layer can be either horizontal or vertical, so the width and height have different implications:

        In a horizontal wg, the fields propagate along the width of the Layer, while the height determines the supported modes.
        In a vertical wg, the fields propagate along the height of the Layer, while the width determines the supported modes.
    """
    index:  float=None
    width:  float=mp.inf
    height: float=mp.inf


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


def create_h_waveguide(center_y: float=None, height: float=None, n: float=None, layer: Layer=None, center_x: float=None, width: float=mp.inf):
    """
    Creates an horitzontal waveguide centered at center_y, parameters can be specified using a Layer.
    If width is not infinite, then center_x must be specified
    """
    if center_y is None:
        raise Exception("To create an horitzontal wg center_y must be specified.")

    if layer is not None:
        width = layer.width
        height = layer.height
        n = layer.index

    if height is None or n is None:
        raise Exception("To create an horitzontal waveguide height and n must be specified.")

    if width == mp.inf or center_x is not None:
        return mp.Block(center=mp.Vector3(0 if center_x is None else center_x, center_y), 
                        size=mp.Vector3(width, height), 
                        material=mp.Medium(index=n))
    else:
        raise Exception("If the width is not mp.inf (infinte waveguide), then center_x must be specified.")


    

def create_v_waveguide(center_x: float=None, width: float=None, n: float=None, layer: Layer=None, center_y: float=None, height: float=mp.inf):
    """
    Creates a vertical waveguide centered at center_x, parameters can be specified using a Layer.
    If height is not infinite, then center_y must be specified
    """
    if center_x is None:
        raise Exception("To create an horitzontal wg center_y must be specified.")

    if layer is not None:
        width = layer.width
        height = layer.height
        n = layer.index

    if width is None or n is None:
        raise Exception("To create an horitzontal waveguide width and n must be specified.")

    if height == mp.inf or center_y is not None:
        return mp.Block(center=mp.Vector3(center_x, 0 if center_y is None else center_y), 
                        size=mp.Vector3(width, height), 
                        material=mp.Medium(index=n))
    else:
        raise Exception("If the height is not mp.inf (infinte waveguide), then center_y must be specified.")

    

def create_h_grating(gr_period, gr_height, gr_duty_cycle, n_cells, layer: Layer,
                     center=mp.Vector3(), gr_up=True, excaved=True, n_ext=1, security_factor=1.1):
    """
    Creates an horitzontal waveguide with a grating.
    gr_period:      Period of the grating
    gr_height:      How deep are the holes (or how tall is the crest if not excaved)
    gr_duty_cycle:  --
    n_cells:        --
    layer:          Specifies the refraction index, the width (if necessary) and if excaved the total height of the layer
                    or if not excaved the height of the layer that will remain when there is no grating
    center:         The center of the layer, this is, if excaved: center of the total (grating + bottom) else: center of bottom
    gr_up:          If true the grating is on the top, else the grating is in the bottom
    excaved:        If true the creation porcess is: There is a layer and holes are excaved. 
                    If false the creation porcess is: There is a layer and cells are put on top
    n_ext:          Refraction index of the "ambient" where will be the holes of the grating
    security_factor:  [>= 1] This parameter doesn't affect to the hole depth- Factor to make the grating holes bigger than the grating height. This is to avoid that the holes are closed due to the resolution of the simulation.
    """

    wg_y = center.y
    length = gr_period * n_cells


    if excaved:
        geometry = [create_h_waveguide(wg_y, layer=layer, center_x=center.x)]
        left = -length/2 + gr_period*gdc/2 + center.x
        gdc = 1-gr_duty_cycle

        if gr_up:
            gy = wg_y + layer.height/2  + (security_factor/2 - 1)*gr_height


            for x in range(n_cells):
                geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*(1-gdc)), gy),
                                        size=mp.Vector3(gr_period*gdc, gr_height*security_factor),
                                        material=mp.Medium(index=n_ext)
                                        )
                                )

        else: 
            gy = wg_y - layer.height/2 + (1 - security_factor/2)*gr_height

            for x in range(n_cells):
                geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*(1-gdc)), gy),
                                        size=mp.Vector3(gr_period*gdc, gr_height*security_factor),
                                        material=mp.Medium(index=n_ext)
                                        )
                                )

    else:
        geometry = [create_h_waveguide(wg_y, layer=layer, center_x=center.x)]
        left = -length/2 + gr_period*(1-gr_duty_cycle)/2 + center.x

        if gr_up:
            gy = wg_y + layer.height/2 + gr_height/2
            for x in range(n_cells):
                geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*gr_duty_cycle), gy),
                                        size=mp.Vector3(gr_period*gr_duty_cycle, gr_height),
                                        material=mp.Medium(index=layer.index)
                                        )
                                )

        else:
            gy = wg_y - layer.height/2 - gr_height/2
            for x in range(n_cells):
                geometry.append(mp.Block(center=mp.Vector3(left + gr_period*(x + 0.5*gr_duty_cycle), gy),
                                        size=mp.Vector3(gr_period*gr_duty_cycle, gr_height),
                                        material=mp.Medium(index=layer.index)
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