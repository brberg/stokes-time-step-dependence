from __future__ import division
import os
import numpy as np
from dolfin import  *

class FenicsMesh:
    def __init__(self, xz_boundary, hires_distance, cell_size, hires_cell_size, inflow_velocity):
        #Distance from the upstream boundary that high resolution mesh should begin.
        self.hires_distance = hires_distance
        #Cell size in low resolution portion of domain.
        self.cell_size = cell_size
        #Cell size in high resolution portion of domain.
        self.hires_cell_size = hires_cell_size
        #Inflow velocity at upstream boundary.
        self.inflow_velocity = inflow_velocity
        #Create initial mesh
        self.make_gmsh_unstructured(xz_boundary)

    def make_gmsh_unstructured(self, xz_boundary):
        "Makes mesh using gmsh by manually writing out a .geo file line by line and then calling gmsh using os.system and converting the resulting file to .xml. Does not remove any points."
        x = xz_boundary[0]
        z = xz_boundary[1]
        num_points = np.shape(x)[0]

        geo_file = open('new-mesh.geo', 'w')

        point_counter = 0
        for j in range(0, num_points):
            geo_file.write("point_list[%g] = newp; \n"%point_counter)
            if x[j] >= self.hires_distance:
                cell_size = self.hires_cell_size
            else:
                cell_size = self.cell_size
            geo_file.write("Point(point_list[%g]) = {%1.13f, %1.13f, 0, %g}; \n"%(point_counter, x[j], z[j], cell_size))
            point_counter += 1

        line_list = []
        line_counter = 1
        line_list.append(line_counter)
        for j in range(point_counter-1):
            geo_file.write("Line(%g) = {point_list[%g], point_list[%g]}; \n"%(line_counter, line_counter-1, line_counter))
            line_counter += 1
            line_list.append(line_counter)
        geo_file.write("Line(%g) = {point_list[%g], point_list[0]}; \n"%(line_counter, line_counter-1))
        line_counter += 1

        line_list_string = ""
        line_list_string += str(line_list[0])
        for j in range (1, len(line_list)):
            line_list_string += ", " + str(line_list[j])

        geo_file.write("Line Loop(%g) = {%s}; \n"%(line_counter, line_list_string))
        line_loop_index = line_counter
        line_counter += 1
        plane_surface_index = line_counter
        geo_file.write("Plane Surface(%g) = {%g}; \n"%(plane_surface_index, line_loop_index))

        geo_file.close()

        os.system("/Applications/Gmsh.app/Contents/MacOS/Gmsh new-mesh.geo -2")
        os.system("dolfin-convert new-mesh.msh new-mesh.xml")
        self.mesh = Mesh("new-mesh.xml")

        [os.remove("new-mesh" + extension) for extension in [".geo", ".msh", ".xml"]]

        self.define_function_spaces()
        self.define_boundaries()

    def define_function_spaces(self):
        "Defines function spaces."
        self.scalar_element = FiniteElement("CG", self.mesh.ufl_cell(), 1)
        self.vector_element = VectorElement("CG", self.mesh.ufl_cell(), 2)
        self.mixed_space = FunctionSpace(self.mesh, MixedElement([self.scalar_element, self.vector_element]))
        self.scalar_space_dg = FunctionSpace(self.mesh, "DG", 0)
        self.vector_space = VectorFunctionSpace(self.mesh, "CG", 2)
        self.tensor_space_dg = TensorFunctionSpace(self.mesh, "DG", 0)

    def define_boundaries(self):
        "Defines subclasses to use with variational form and boundary conditions."
        boundary_meshfunction = MeshFunction('size_t', self.mesh, self.mesh.topology().dim() - 1, 0)
        boundary_meshfunction.set_all(0)

        #Define class for underwater portion of boundary.
        class Underwater(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[1] <= 0.0
        
        Underwater().mark(boundary_meshfunction, 2)

        #Define class for upstream boundary.
        class Upstream(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 0.0, DOLFIN_EPS)

        Upstream().mark(boundary_meshfunction, 1)

        #Set upstream Dirichlet boundary condition.
        bc_velocity = DirichletBC(self.mixed_space.sub(1).sub(0), Constant(self.inflow_velocity), boundary_meshfunction, 1)

        self.boundary_conditions = [bc_velocity]
        self.ds = Measure("ds")(subdomain_data=boundary_meshfunction)

        #Save boundary to file.
        File("boundaries.pvd") << boundary_meshfunction