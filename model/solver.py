from __future__ import division
from meshing import *
from viscosity import *
from ufl import nabla_div

#Only display fenics warnings. Supress simple information.
set_log_level(20)

class SimulationSolver:

    def __init__(self, xz_boundary, solver_tolerance, creep_parameter, hires_distance, cell_size, hires_cell_size, inflow_velocity, include_inertial_term, rho_i=900, rho_w=1000, g=9.81, seconds_per_year=3.154E7):

        #Initial time at start of simulation.
        self.time = 0.0

        #Include inertial term or not.
        self.include_inertial_term = include_inertial_term
        
        #Sub-Classes
        self.FenicsMesh = FenicsMesh(xz_boundary, hires_distance, cell_size, hires_cell_size, inflow_velocity)
        self.ViscosityCalculator = ViscosityCalculator(creep_parameter)

        #Numerical Constants
        self.solver_tolerance = solver_tolerance
        self.seconds_per_year = Constant(seconds_per_year)

        #Physical Constants
        self.rho_i = rho_i
        self.rho_w = rho_w
        self.g = g
        self.body_force = Constant((0, -rho_i*g))

        #Initial velocity guesses used in solver loops
        self.iteration_velocity = interpolate(Constant((self.FenicsMesh.inflow_velocity, 0.0)), self.FenicsMesh.vector_space)
        self.step_velocity = interpolate(Constant((self.FenicsMesh.inflow_velocity, 0.0)), self.FenicsMesh.vector_space)

        #Variable to store solver solution
        self.solver_solution = Function(self.FenicsMesh.mixed_space)

        #Water pressure expression to use in variational form.
        class WaterPressure(UserExpression):
            def eval(self, value, x):
                if x[1] <= 0.0:
                    value[0] = rho_w*g*(0.0-x[1])
                else:
                    value[0] = 0.0
        self.water_pressure_expression = WaterPressure(degree=1)

        #Setup files for saving velocity and pressure solutions
        self.ux_file = File("ux.pvd", "compressed")
        self.uz_file = File("uz.pvd", "compressed")
        self.p_file = File("p.pvd", "compressed")

        #Setup files for saving diagnostic quantities
        self.strainrateE_file = File("strainrateE.pvd", "compressed")
        self.viscosity_file = File("viscosity.pvd", "compressed")
        self.sigmaP_file = File("sigmaP.pvd", "compressed")
        self.tauMax_file = File("tauMax.pvd", "compressed")

    def b_function(self, v, q):
        "Function used in variational form."
        return inner(nabla_div(v), q)*dx

    def norm_calculator(self, velocity):
        "Inner product for use in determining if the solver has converged."
        return assemble(inner(velocity, velocity)*dx)

    def save_solution(self):
        "Calculates diagnostic fields from velocity and pressure solutions and writes results to file."
        pressure_solution, velocity_solution = self.solver_solution.split()
        velocity_x_solution, velocity_z_solution = velocity_solution.split()

        #Strain Rate Tensor
        strainrate = strain_rate(velocity_solution)

        #Effective Strain Rate
        epsilonE_ufl = (0.5*strainrate[0,0]**2+0.5*strainrate[1,1]**2+strainrate[0,1]**2)**(1/2)
        epsilonE = project(epsilonE_ufl, self.FenicsMesh.scalar_space_dg)

        #Viscosity
        viscosity = project(self.ViscosityCalculator.compute_viscosity(velocity_solution), self.FenicsMesh.scalar_space_dg)

        #Cauchy Stress
        sigma = 2*self.ViscosityCalculator.compute_viscosity(velocity_solution)*strain_rate(velocity_solution)-pressure_solution*\
                Identity(self.FenicsMesh.tensor_space_dg.ufl_cell().topological_dimension())
        ddelta = ((sigma[0,0] - sigma[1,1])**(2) + 4*(sigma[0,1])**(2))**(1/2)
        sigmaP_ufl = 0.5*(sigma[0,0] + sigma[1,1] + ddelta)
        sigmaP = project(sigmaP_ufl, self.FenicsMesh.scalar_space_dg)

        #Deviatoric Stress
        tau = 2*self.ViscosityCalculator.compute_viscosity(velocity_solution)*strain_rate(velocity_solution)
        tau_max_ufl = (0.5*tau[0,0]**2+0.5*tau[1,1]**2+tau[0,1]**2)**(1/2)
        tau_max = project(tau_max_ufl, self.FenicsMesh.scalar_space_dg)

        #Rename variables for readability. 
        velocity_x_solution.rename('ux', 'ux')
        velocity_z_solution.rename('uz', 'uz')
        pressure_solution.rename('p', 'p')
        epsilonE.rename('e', 'e')
        viscosity.rename('viscosity', 'viscosity')
        sigmaP.rename('sP', 'sP')
        tau_max.rename('tM', 'tM')

        #Save solutions to file.
        self.ux_file << (velocity_x_solution, self.time)
        self.uz_file << (velocity_z_solution, self.time)
        self.p_file << (pressure_solution, self.time)
        self.strainrateE_file << (epsilonE, self.time)
        self.viscosity_file << (viscosity, self.time)
        self.sigmaP_file << (sigmaP, self.time)
        self.tauMax_file << (tau_max, self.time)

    def run(self, duration, dt):
        "Run simulation for a given duration with a given time step."
        while self.time <= duration:
            #Solve model and save solutions to file.
            self.solve_model(dt)
            self.save_solution()

            #Update mesh coordinates for next time step based on velocity solution
            pressure_solution, velocity_solution = self.solver_solution.split()
            velocity_x_solution, velocity_z_solution = velocity_solution.split()
            ux = velocity_x_solution.compute_vertex_values()
            uz = velocity_z_solution.compute_vertex_values()
            self.FenicsMesh.mesh.coordinates()[:, 0] = self.FenicsMesh.mesh.coordinates()[:, 0] + ux*dt
            self.FenicsMesh.mesh.coordinates()[:, 1] = self.FenicsMesh.mesh.coordinates()[:, 1] + uz*dt

            #Assign velocity solution to step_velocity for use calculation of inertial term at the next time step.
            assign(self.step_velocity, velocity_solution)
            self.time += dt

            print("Current time in years is " + str(self.time) + "/" + str(duration) + ".")

    def solve_model(self, dt):
        "Solve model with a given time step."
        #Define trial and test functions.
        trial_functions = TrialFunction(self.FenicsMesh.mixed_space)
        test_functions = TestFunction(self.FenicsMesh.mixed_space)
        (p, u) = split(trial_functions)
        (q, v) = split(test_functions)

        #Define necessary unit vectors.
        z_hat = interpolate(Constant((0.0, 1.0)), self.FenicsMesh.vector_space)
        facet_normal = FacetNormal(self.FenicsMesh.mesh)
        n_hat = as_vector([facet_normal[0], facet_normal[1]])

        norm = 1.0
        while float(norm) >= self.solver_tolerance:
            #Compute viscosity for current velocity guess.
            iteration_viscosity = self.ViscosityCalculator.compute_viscosity(self.iteration_velocity)

            #Define proper variational form based on whether the inertial term is included.
            if self.include_inertial_term:
                variational_form = 2*inner(iteration_viscosity*strain_rate(u), strain_rate(v))*dx - self.b_function(v, p) + self.b_function(u, q) \
                    + inner(self.rho_i*1/(self.seconds_per_year*dt)*(u-self.step_velocity), v)*dx \
                    - inner(self.body_force, v)*dx \
                    - inner(self.rho_w*self.g*dt*dot(u, z_hat)*n_hat, v)*self.FenicsMesh.ds(2) \
                    + inner(self.water_pressure_expression*n_hat, v)*self.FenicsMesh.ds(2)
            else:
                variational_form = 2*inner(iteration_viscosity*strain_rate(u), strain_rate(v))*dx - self.b_function(v, p) + self.b_function(u, q) \
                    - inner(self.body_force, v)*dx \
                    - inner(self.rho_w*self.g*dt*dot(u, z_hat)*n_hat, v)*self.FenicsMesh.ds(2) \
                    + inner(self.water_pressure_expression*n_hat, v)*self.FenicsMesh.ds(2)

            solve(lhs(variational_form) == rhs(variational_form), self.solver_solution, self.FenicsMesh.boundary_conditions)

            #Calculate velocity field change from previous time step and calculate norm to determine solver convergence.
            pressure_solution, velocity_solution = self.solver_solution.split()
            velocity_x_solution, velocity_z_solution = velocity_solution.split()
            velocity_solution_change = velocity_solution - self.iteration_velocity
            norm = self.norm_calculator(velocity_solution_change)/self.norm_calculator(velocity_solution)
            print("Current solver norm is " + str(float(norm)) + ".")

            #Update velocity for next iteration
            assign(self.iteration_velocity, velocity_solution)