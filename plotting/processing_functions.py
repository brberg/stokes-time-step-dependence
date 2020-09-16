from __future__ import division
import numpy as np
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def get_cell_values(data_name, folders):
    "For a set of folders, gets the cell value of the given data name from a saved .vtu file."
    main_directory = os.getcwd()
    values_array = []
    for folder in folders:
        values = []
        os.chdir(folder)

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(data_name + '{0:06d}.vtu'.format(int(0)))
        reader.Update()
        data = reader.GetOutput()

        f = vtk_to_numpy(data.GetCellData().GetArray(0))

        values.append(f)
        values_array.append(np.asarray(values))

        os.chdir(main_directory)

    return values_array

def get_point_values(data_name, folders):
    "For a set of folders, gets the point value of the given data name from a saved .vtu file."
    main_directory = os.getcwd()
    values_array = []
    for folder in folders:
        values = []
        os.chdir(folder)

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(data_name + '{0:06d}.vtu'.format(int(0)))
        reader.Update()
        data = reader.GetOutput()

        f = vtk_to_numpy(data.GetPointData().GetArray(0))

        values.append(f)
        values_array.append(np.asarray(values))

        os.chdir(main_directory)

    return values_array

def l2_norm(field):
    "Calculates the L2 norm of the input data."
    return np.sqrt(np.sum(np.square(field)))/np.sqrt((len(field)-1))

def field_max(field):
    "Calculates the maximum of the input data."
    return np.max(field)