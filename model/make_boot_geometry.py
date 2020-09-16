from __future__ import division
import numpy as np

def make_boot_geometry(length, thickness, boot_length, boot_depth, hires, lores):
    "Exports counterclockwise ordered x and z coordinates of a boot geometry to use with meshing."

    #Make x and z coordiantes for the upstream part of the ice tongue, set to buoyancy.
    x1 = np.linspace(0, length-boot_length-hires, (length-boot_length-hires)//lores)
    z_bottom1 = (-0.9*thickness)*np.ones(len(x1))
    z_surface1 = z_bottom1 + thickness

    #Make x and z coordiantes for the downstream part of the ice tongue, set to buoyancy.
    x2 = np.linspace(length-boot_length, length, boot_length//hires)
    z_bottom2 = (-0.9*thickness)*np.ones(len(x2))
    z_surface2 = z_bottom2 + thickness - boot_depth

    #Merge two sets of coordinates
    x = np.append(x1, x2)
    z_bottom = np.append(z_bottom1, z_bottom2)
    z_surface = np.append(z_surface1, z_surface2)

    #Merge bottom and surface coordinates
    x_all = np.append(x, x[::-1])
    z_all = np.append(z_bottom, z_surface[::-1])

    return [x_all, z_all]