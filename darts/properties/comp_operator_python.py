import numpy as np
from darts.engines import *
from darts.physics import *
from properties.comp_params import *

import os.path as osp

physics_name = osp.splitext(osp.basename(__file__))[0]

# Define our own operator evaluator class
class AccFluxGravityEvaluator(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()  # Initialize base-class
        # Store your input parameters in self here, and initialize other parameters here in self
        self.num_comp = property_container.nc
        self.min_z = property_container.min_z
        self.phases = property_container.phase_name
        self.components = property_container.component_name

        self.property = property_container

        self.c_r = 1e-5
        self.p_ref = 1
        self.P_entry = 0.2

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))
        zc = self.comp_out_of_bounds(zc)

        # Perform Flash procedure here:
        (x, y, V) = self.property.flash_ev.flash(state)
        Mt_aq = 0
        Mt_g = 0
        # molar weight of mixture
        for i in range(self.num_comp):
            Mt_aq = Mt_aq + props(self.components[i], 'Mw') * x[i]
            Mt_g = Mt_g + props(self.components[i], 'Mw') * y[i]

        # print('Mt_aq:', Mt_aq, ', Mt_gas:', Mt_g)
        rho_g = 0
        mu_g = 1

        rho_aq = 0
        mu_aq = 1

        rho_g_m = 0
        rho_aq_m = 0

        """" PROPERTIES Evaluation """
        if V <= 0:
            sg = 0
            # single aqueous phase
            x = zc
            rho_aq = self.property.density_ev['Aq'].evaluate(x)      # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()      # output in [cp]
            rho_aq_m = rho_aq / Mt_aq                                # [kg/m3] / [kg/Kmol] = [Kmol/m3]
            kr_aq = self.property.rel_perm_ev['Aq'].evaluate(1-sg)
            kr_g = 0

        elif V >= 1:
            sg = 1
            # single vapor phase
            y = zc
            rho_g = self.property.density_ev['gas'].evaluate(pressure)        # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()               # in [cp]
            rho_g_m = rho_g / Mt_g                                            # [kg/m3] * [kg/Kmol]^-1 = [Kmol/m3]
            kr_g = self.property.rel_perm_ev['gas'].evaluate(sg)
            kr_aq = 0

        else:
            # two phases
            rho_aq = self.property.density_ev['Aq'].evaluate(x)      # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()      # output in [cp]
            rho_g = self.property.density_ev['gas'].evaluate(pressure)           # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()       # in [cp]
            rho_aq_m = rho_aq / Mt_aq
            rho_g_m = rho_g / Mt_g
            sg = rho_aq_m / (rho_g_m / V - rho_g_m + rho_aq_m)
            kr_aq = self.property.rel_perm_ev['Aq'].evaluate(1-sg)
            kr_g = self.property.rel_perm_ev['gas'].evaluate(sg)

        """ CONSTRUCT OPERATORS HERE """
        # Alpha operator represents accumulation term:
        num_alpha_op = self.num_comp
        for i in range(num_alpha_op):
            values[i] = (1+self.c_r*(pressure-self.p_ref))*((sg * rho_g_m * y[i]) + (rho_aq_m * (1-sg) * x[i]))

        # Beta operator represents flux term:
        # Gas phase
        num_beta_op_gas = self.num_comp + 1 + 1     # nr_comp + grav(gas) + pc(gas)
        values[num_alpha_op+0] = rho_g              # gravity operator of gas phase
        values[num_alpha_op+1] = 0                  # capillary operator of gas phase
        for i in range(self.num_comp):
            values[i+num_alpha_op+2] = y[i] * rho_g_m * kr_g/mu_g

        # Aqueous phase
        num_beta_op_aq = self.num_comp + 1 + 1
        values[num_alpha_op + num_beta_op_gas + 0] = rho_aq     # gravity operator of aq phase
        values[num_alpha_op + num_beta_op_gas + 1] = 0          # capillary operator of aq phase

        for i in range(self.num_comp):
            values[i+num_alpha_op+num_beta_op_gas+2] = x[i] * rho_aq_m * kr_aq/mu_aq

        """ add diffusion term """
        D = 0.0001728 # m2/day

        num_diff_op = self.num_comp
        for i in range(self.num_comp):
            compr = (1+self.c_r*(pressure-self.p_ref))          # compressible rock
            # values[i + num_alpha_op + num_beta_op_gas + num_beta_op_aq] = compr * (rho_g_m * sg * D * y[i] + rho_aq_m * (1-sg) * D * x[i])
            values[i + num_alpha_op + num_beta_op_gas + num_beta_op_aq] = 0

        return 0

class AccFluxGravityEvaluatorWell(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()  # Initialize base-class
        # Store your input parameters in self here, and initialize other parameters here in self
        self.num_comp = property_container.nc
        self.min_z = property_container.min_z
        self.phases = property_container.phase_name
        self.components = property_container.component_name

        self.property = property_container

        self.c_r = 1e-5
        self.p_ref = 1
        self.P_entry = 0.2

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))
        zc = self.comp_out_of_bounds(zc)

        # Perform Flash procedure here:
        (x, y, V) = self.property.flash_ev.flash(state)
        Mt_aq = 0
        Mt_g = 0
        # molar weight of mixture
        for i in range(self.num_comp):
            Mt_aq = Mt_aq + props(self.components[i], 'Mw') * x[i]
            Mt_g = Mt_g + props(self.components[i], 'Mw') * y[i]

        # print('Mt_aq:', Mt_aq, ', Mt_gas:', Mt_g)
        rho_g = 0
        mu_g = 1

        rho_aq = 0
        mu_aq = 1

        rho_g_m = 0
        rho_aq_m = 0

        """" PROPERTIES Evaluation """
        if V <= 0:
            sg = 0
            # single aqueous phase
            x = zc
            rho_aq = self.property.density_ev['Aq'].evaluate(x)      # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()      # output in [cp]
            rho_aq_m = rho_aq / Mt_aq                                # [kg/m3] / [kg/Kmol] = [Kmol/m3]
            kr_aq = self.property.rel_perm_ev['Aq'].evaluate(1-sg)
            kr_g = 0

        elif V >= 1:
            sg = 1
            # single vapor phase
            y = zc
            rho_g = self.property.density_ev['gas'].evaluate(pressure)        # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()               # in [cp]
            rho_g_m = rho_g / Mt_g                                            # [kg/m3] * [kg/Kmol]^-1 = [Kmol/m3]
            kr_g = self.property.rel_perm_ev['gas'].evaluate(sg)
            kr_aq = 0

        else:
            # two phases
            rho_aq = self.property.density_ev['Aq'].evaluate(x)      # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()      # output in [cp]
            rho_g = self.property.density_ev['gas'].evaluate(pressure)           # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()       # in [cp]
            rho_aq_m = rho_aq / Mt_aq
            rho_g_m = rho_g / Mt_g
            sg = rho_aq_m / (rho_g_m / V - rho_g_m + rho_aq_m)
            kr_aq = self.property.rel_perm_ev['Aq'].evaluate(1-sg)
            kr_g = self.property.rel_perm_ev['gas'].evaluate(sg)

        """ CONSTRUCT OPERATORS HERE """
        # Alpha operator represents accumulation term:
        num_alpha_op = self.num_comp
        for i in range(num_alpha_op):
            values[i] = (1+self.c_r*(pressure-self.p_ref))*((sg * rho_g_m * y[i]) + (rho_aq_m * (1-sg) * x[i]))

        # Beta operator represents flux term:
        # Gas phase
        num_beta_op_gas = self.num_comp + 1 + 1     # nr_comp + grav(gas) + pc(gas)
        values[num_alpha_op+0] = 0                  # gravity operator of gas phase
        values[num_alpha_op+1] = 0                  # capillary operator of gas phase
        for i in range(self.num_comp):
            values[i+num_alpha_op+2] = y[i] * rho_g_m * kr_g/mu_g

        # Aqueous phase
        num_beta_op_aq = self.num_comp + 1 + 1
        values[num_alpha_op + num_beta_op_gas + 0] = 0          # gravity operator of aq phase
        values[num_alpha_op + num_beta_op_gas + 1] = 0          # capillary operator of aq phase

        for i in range(self.num_comp):
            values[i+num_alpha_op+num_beta_op_gas+2] = x[i] * rho_aq_m * kr_aq/mu_aq

        """ add diffusion term """
        D = 0.0001728 # m2/day

        num_diff_op = self.num_comp
        for i in range(self.num_comp):
            compr = (1+self.c_r*(pressure-self.p_ref))          # compressible rock
            # values[i + num_alpha_op + num_beta_op_gas + num_beta_op_aq] = compr * (rho_g_m * sg * D * y[i] + rho_aq_m * (1-sg) * D * x[i])
            values[i + num_alpha_op + num_beta_op_gas + num_beta_op_aq] = 0

        return 0


class RateEvaluator(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()  # Initialize base-class
        # Store your input parameters in self here, and initialize other parameters here in self
        self.num_comp = property_container.nc
        self.min_z = property_container.min_z
        self.phases = property_container.phase_name
        self.components = property_container.component_name

        self.property = property_container

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))
        zc = self.comp_out_of_bounds(zc)

        # Perform Flash procedure here:
        (x, y, V) = self.property.flash_ev.flash(state)
        Mt_aq = 0
        Mt_g = 0
        # molar weight of mixture
        for i in range(self.num_comp):
            Mt_aq = Mt_aq + props(self.components[i], 'Mw') * x[i]
            Mt_g = Mt_g + props(self.components[i], 'Mw') * y[i]

        # print('Mt_aq:', Mt_aq, ', Mt_gas:', Mt_g)
        rho_g = 0
        mu_g = 1

        rho_aq = 0
        mu_aq = 1

        rho_g_m = 0
        rho_aq_m = 0

        """" PROPERTIES Evaluation """
        if V <= 0:
            sg = 0
            # single aqueous phase
            x = zc
            rho_aq = self.property.density_ev['Aq'].evaluate(x)      # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()      # output in [cp]
            rho_aq_m = rho_aq / Mt_aq                                # [kg/m3] / [kg/Kmol]^-1 = [Kmol/m3]
            kr_aq = self.property.rel_well_perm_ev['Aq'].evaluate(1-sg)
            kr_g = 0

        elif V >= 1:
            sg = 1
            # single vapor phase
            y = zc
            rho_g = self.property.density_ev['gas'].evaluate(pressure)       # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()              # in [cp]
            rho_g_m = rho_g / Mt_g                                           # [kg/m3] / [kg/Kmol]^-1 = [Kmol/m3]
            kr_g = self.property.rel_well_perm_ev['gas'].evaluate(sg)
            kr_aq = 0

        else:
            # two phases
            rho_aq = self.property.density_ev['Aq'].evaluate(x)                  # output in [kg/m3]
            mu_aq = self.property.viscosity_ev['Aq'].evaluate()                  # output in [cp]
            rho_g = self.property.density_ev['gas'].evaluate(pressure)           # in [kg/m3]
            mu_g = self.property.viscosity_ev['gas'].evaluate()                  # in [cp]
            rho_aq_m = rho_aq / Mt_aq
            rho_g_m = rho_g / Mt_g
            sg = rho_aq_m / (rho_g_m / V - rho_g_m + rho_aq_m)
            kr_aq = self.property.rel_well_perm_ev['Aq'].evaluate(1-sg)
            kr_g = self.property.rel_well_perm_ev['gas'].evaluate(sg)

        """ CONSTRUCT RATE OPERATORS HERE """
        num_rate_op = self.num_comp  # two p two c just for this case, otherwise need to import phases number

        # step-1
        flux = np.zeros(num_rate_op)
        for i in range(num_rate_op):
            flux[i] = (rho_g_m * kr_g * y[i] / mu_g) + (rho_aq_m * x[i] * kr_aq / mu_aq)
        # step-2
        flux_sum = np.sum(flux)
        # step-3
        total_density = sg * rho_g_m + (1 - sg) * rho_aq_m
        # step-4
        values[0] = sg * flux_sum / total_density
        values[1] = (1 - sg) * flux_sum / total_density

        return 0

""" only used for plotting in main.py """
class properties_evaluator(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()  # Initialize base-class
        # Store your input parameters in self here, and initialize other parameters here in self
        self.num_comp = property_container.nc
        self.min_z = property_container.min_z
        self.phases = property_container.phase_name
        self.components = property_container.component_name

        self.property = property_container

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))
        zc = self.comp_out_of_bounds(zc)

        # Perform Flash procedure here:
        (x, y, V) = self.property.flash_ev.flash(state)
        Mt_aq = 0
        Mt_g = 0
        # molar weight of mixture
        for i in range(self.num_comp):
            Mt_aq = Mt_aq + props(self.components[i], 'Mw') * x[i]
            Mt_g = Mt_g + props(self.components[i], 'Mw') * y[i]

        # print('Mt_aq:', Mt_aq, ', Mt_gas:', Mt_g)
        rho_g = 0
        mu_g = 1

        rho_aq = 0
        mu_aq = 1

        rho_g_m = 0
        rho_aq_m = 0

        """" PROPERTIES Evaluation """
        if V <= 0:
            sg = 0
            # single aqueous phase
            x = zc
            rho_aq = self.property.density_ev['Aq'].evaluate(x)  # output in [kg/m3]

        elif V >= 1:
            sg = 1
            # single vapor phase
            y = zc
            rho_g = self.property.density_ev['gas'].evaluate(pressure)  # in [kg/m3]

        else:
            # two phases
            state = value_vector(np.append(np.array([pressure]), x[:-1]))
            rho_aq = self.property.density_ev['Aq'].evaluate(x)  # output in [kg/m3]
            state = [pressure, [y[:-1]]]
            rho_g = self.property.density_ev['gas'].evaluate(pressure)  # in [kg/m3]
            rho_aq_m = rho_aq / Mt_aq
            rho_g_m = rho_g / Mt_g
            sg = rho_aq_m / (rho_g_m / V - rho_g_m + rho_aq_m)  # saturation using [Kmol/m3]

        return sg, x[0]
