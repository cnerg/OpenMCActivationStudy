import openmc

def make_spherical_shell(material, thickness, inner_radius):  
    '''
    Creates a spherical shell with its own material and inner/outer radius.
    
    inputs:
        material: OpenMC Material object/iterable of OpenMC Material
        thickness: radial thickness (float) of material
        inner_radius : inner radius of material
    
    outputs:
        geometry: OpenMC Geometry object
    '''
    inner_sphere = openmc.Sphere(r = inner_radius)
    cells = [openmc.Cell(fill = None, region = -inner_sphere)]
    outer_radius = inner_radius + thickness
    outer_sphere = openmc.Sphere(r = outer_radius, boundary_type = 'vacuum')
    cells.append(openmc.Cell(fill = material, region = +inner_sphere & -outer_sphere))   
    cells.append(openmc.Cell(fill = None, region = +outer_sphere)) 
    geometry = openmc.Geometry(cells)    
    return geometry