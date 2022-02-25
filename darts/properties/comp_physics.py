import numpy as np
from darts.engines import *
from darts.physics import *
from properties.comp_operator_python import *
from properties.comp_properties import *


# Define our own operator evaluator class
class CO2Brine:
    def __init__(self, timer, n_points, min_p, max_p, min_z, max_z, property_container):

        # Obtain properties from user input during initialization:
        self.timer = timer.node["simulation"]
        self.components = property_container.component_name
        self.nr_components = property_container.nc
        self.phases = property_container.phase_name
        self.nr_phases = property_container.n_phases
        self.n_vars = self.nr_components
        self.vars = ['pressure', 'zCO2']

        """ Name of interpolation method and engine used for this physics: """
        # engine including gravity term
        engine_name = eval("engine_nc_cg_dif_cpu%d_%d" % (self.nr_components, self.nr_phases))
        self.nr_ops = self.nr_components + self.nr_components*self.nr_phases + self.nr_phases + self.nr_phases + self.nr_components

        acc_flux_itor_name = eval("operator_set_interpolator_i_d_%d_%d" % (self.nr_components, self.nr_ops))
        rate_interpolator_name = eval("operator_set_interpolator_i_d_%d_%d" % (self.nr_components, self.nr_phases))

        acc_flux_itor_name_long = eval("operator_set_interpolator_l_d_%d_%d" % (self.nr_components, self.nr_ops))
        rate_interpolator_name_long = eval("operator_set_interpolator_l_d_%d_%d" % (self.nr_components, self.nr_phases))

        """ reservoir evaluator """
        # Initialize main evaluator - Reservoir
        self.acc_flux_etor = AccFluxGravityEvaluator(property_container=property_container)
        self.acc_flux_etor_well = AccFluxGravityEvaluatorWell(property_container=property_container)
        # Initialize table entries (nr of points, axis min, and axis max):
        # nr_of_points for [pres, comp1, ..., compN-1]:
        self.acc_flux_etor.axis_points = index_vector([n_points, n_points])
        # axis_min for [pres, comp1, ..., compN-1]:
        self.acc_flux_etor.axis_min = value_vector([min_p, min_z])
        # axis_max for [pres, comp1, ..., compN-1]
        self.acc_flux_etor.axis_max = value_vector([max_p, max_z])

        # Create actual accumulation and flux interpolator:
        try:
            self.acc_flux_itor = acc_flux_itor_name(self.acc_flux_etor, self.acc_flux_etor.axis_points,
                                                    self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)
            self.acc_flux_itor_well = acc_flux_itor_name(self.acc_flux_etor_well, self.acc_flux_etor.axis_points,
                                                    self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)
        except RuntimeError:
            self.acc_flux_itor = acc_flux_itor_name_long(self.acc_flux_etor, self.acc_flux_etor.axis_points,
                                                         self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)
            self.acc_flux_itor_well = acc_flux_itor_name_long(self.acc_flux_etor_well, self.acc_flux_etor.axis_points,
                                                         self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)

        # set up timers
        self.timer.node["jacobian assembly"] = timer_node()
        self.timer.node["jacobian assembly"].node["interpolation"] = timer_node()
        # reservoir blocks
        self.timer.node["jacobian assembly"].node["interpolation"].node["acc flux interpolation"] = timer_node()
        self.acc_flux_itor.init_timer_node(self.timer.node["jacobian assembly"].node["interpolation"].node["acc flux interpolation"])
        # well
        self.timer.node["jacobian assembly"].node["interpolation"].node["acc flux w interpolation"] = timer_node()
        self.acc_flux_itor_well.init_timer_node(
            self.timer.node["jacobian assembly"].node["interpolation"].node["acc flux w interpolation"])

        # create engine according to physics selected
        self.engine = engine_name()

        # Create rate evaluator and interpolator: - for reservoir
        self.rate_etor = RateEvaluator(property_container=property_container)
        try:
            self.rate_itor = rate_interpolator_name(self.rate_etor, self.acc_flux_etor.axis_points,
                                                    self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)

        except RuntimeError:
            self.rate_itor = rate_interpolator_name_long(self.rate_etor, self.acc_flux_etor.axis_points,
                                                         self.acc_flux_etor.axis_min, self.acc_flux_etor.axis_max)

        # set up timers
        self.timer.node["jacobian assembly"].node["interpolation"].node["rate interpolation"] = timer_node()
        self.rate_itor.init_timer_node(
            self.timer.node["jacobian assembly"].node["interpolation"].node["rate interpolation"])

        # create engine according to physics selected
        self.engine = engine_name()

        # define well control factories
        # Injection wells (upwind method requires both bhp and inj_stream for bhp controlled injection wells):
        self.new_bhp_inj = lambda bhp, inj_stream: bhp_inj_well_control(bhp, value_vector(inj_stream))
        self.new_rate_gas_inj = lambda rate, inj_stream: rate_inj_well_control(self.phases, 0, self.nr_components,
                                                                               self.nr_components, rate,
                                                                               value_vector(inj_stream), self.rate_itor)
        self.new_rate_water_inj = lambda rate, inj_stream: rate_inj_well_control(self.phases, 1, self.nr_components,
                                                                               self.nr_components, rate,
                                                                               value_vector(inj_stream), self.rate_itor)

        # Production wells:
        self.new_bhp_prod = lambda bhp: bhp_prod_well_control(bhp)
        self.new_rate_gas_prod = lambda rate: rate_prod_well_control(self.phases, 0, self.nr_components,
                                                                     self.nr_components,
                                                                     rate, self.rate_itor)
        self.new_rate_water_prod = lambda rate: rate_prod_well_control(self.phases, 1, self.nr_components,
                                                                       self.nr_components,
                                                                       rate, self.rate_itor)

    # Define some class methods:
    def init_wells(self, wells):
        for w in wells:
            assert isinstance(w, ms_well)
            w.init_rate_parameters(self.nr_components, self.phases, self.rate_itor)

    def set_uniform_initial_conditions(self, mesh, uniform_pressure, uniform_composition: list):
        assert isinstance(mesh, conn_mesh)

        nb = mesh.n_blocks
        """ Uniform Initial conditions """
        # set initial pressure
        pressure = np.array(mesh.pressure, copy=False)
        pressure.fill(uniform_pressure)

        # set initial composition
        mesh.composition.resize(nb * (self.nr_components - 1))
        composition = np.array(mesh.composition, copy=False)
        # composition[:] = np.array(uniform_composition)
        if self.nr_components == 2:
            for c in range(self.nr_components - 1):                               # Denis
                composition[c::(self.nr_components - 1)] = uniform_composition[:] # Denis
        else:
            for c in range(self.nr_components - 1):  # Denis
                composition[c::(self.nr_components - 1)] = uniform_composition[c]  # Denis

    def set_boundary_conditions(self, mesh, uniform_pressure, uniform_composition):
        assert isinstance(mesh, conn_mesh)

        # Class methods which can create constant pressure and composition boundary condition:
        pressure = np.array(mesh.pressure, copy=False)
        pressure.fill(uniform_pressure)

        mesh.composition.resize(mesh.n_blocks * (self.nr_components - 1))
        composition = np.array(mesh.composition, copy=False)
        for c in range(self.nr_components - 1):
            composition[c::(self.nr_components - 1)] = uniform_composition[c]
