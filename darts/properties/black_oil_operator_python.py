from darts.engines import *
from properties.black_oil_properties import *

# Define operator evaluator class
class black_oil_acc_flux_evaluator_python(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()
        self.property = property_container

    def evaluate(self, state, values):
        pbub = self.property.pb_ev.evaluate(state)
        Rs = self.property.rs_ev.evaluate(state, pbub)
        xgo = self.property.xgo_ev.evaluate(Rs)
        oil_dens = self.property.density_ev['oil'].evaluate(state, pbub, xgo)
        wat_dens = self.property.density_ev['water'].evaluate(state)
        gas_dens = self.property.density_ev['gas'].evaluate(state)
        oil_visco = self.property.viscosity_ev['oil'].evaluate(state, pbub)
        wat_visco = self.property.viscosity_ev['water'].evaluate(state)
        gas_visco = self.property.viscosity_ev['gas'].evaluate(state)
        wat_sat = self.property.sat_ev['water'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        oil_sat = self.property.sat_ev['oil'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        gas_sat = 1 - wat_sat - oil_sat
        wat_relp = self.property.rel_perm_ev['water'].evaluate(wat_sat)
        krow = self.property.rel_perm_ev['OW'].evaluate(wat_sat)
        gas_relp = self.property.rel_perm_ev['gas'].evaluate(gas_sat)
        krog = self.property.rel_perm_ev['GO'].evaluate(gas_sat)
        oil_relp = self.property.rel_perm_ev['oil'].evaluate(wat_sat, gas_sat, krow, krog)
        rock_cp = self.property.rock_compress_ev.evaluate(state)

        # accumulation operator
        values[0] = rock_cp * xgo * oil_dens * oil_sat + rock_cp * gas_sat * gas_dens
        values[1] = rock_cp * (1 - xgo) * oil_dens * oil_sat
        values[2] = rock_cp * wat_sat * wat_dens
        # flux operator
        values[3] = xgo * oil_dens * (oil_relp / oil_visco) + gas_dens * (gas_relp / gas_visco)
        values[4] = (1 - xgo) * oil_dens * (oil_relp / oil_visco)
        values[5] = wat_dens * (wat_relp / wat_visco)

        # print('State:', state, 'Ops', values)

        return 0

class black_oil_acc_flux_capillary_evaluator_python(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()
        self.property = property_container

    def evaluate(self, state, values):
        pbub = self.property.pb_ev.evaluate(state)
        Rs = self.property.rs_ev.evaluate(state, pbub)
        xgo = self.property.xgo_ev.evaluate(Rs)
        oil_dens = self.property.density_ev['oil'].evaluate(state, pbub, xgo)
        wat_dens = self.property.density_ev['water'].evaluate(state)
        gas_dens = self.property.density_ev['gas'].evaluate(state)
        oil_visco = self.property.viscosity_ev['oil'].evaluate(state, pbub)
        wat_visco = self.property.viscosity_ev['water'].evaluate(state)
        gas_visco = self.property.viscosity_ev['gas'].evaluate(state)
        wat_sat = self.property.sat_ev['water'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        oil_sat = self.property.sat_ev['oil'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        gas_sat = 1 - wat_sat - oil_sat
        wat_relp = self.property.rel_perm_ev['water'].evaluate(wat_sat)
        krow = self.property.rel_perm_ev['OW'].evaluate(wat_sat)
        gas_relp = self.property.rel_perm_ev['gas'].evaluate(gas_sat)
        krog = self.property.rel_perm_ev['GO'].evaluate(gas_sat)
        oil_relp = self.property.rel_perm_ev['oil'].evaluate(wat_sat, gas_sat, krow, krog)
        rock_cp = self.property.rock_compress_ev.evaluate(state)

        # accumulation operator
        values[0] = rock_cp * xgo * oil_dens * oil_sat + rock_cp * gas_sat * gas_dens
        values[1] = rock_cp * (1 - xgo) * oil_dens * oil_sat
        values[2] = rock_cp * wat_sat * wat_dens

        # flux operator
        # (1) gas phase
        values[3] = gas_dens             # gas density operator
        values[4] = 0                    # gas and oil, pg=po+pcgo=po-(-pcgo)
        values[5] = gas_dens * (gas_relp / gas_visco)   # gas component in gas phase
        values[6] = 0                                   # oil component in gas phase
        values[7] = 0                                   # water component in gas phase

        # (2) oil phase
        values[8] = oil_dens             # oil density operator
        values[9] = 0                    # reference phase, pc = 0
        values[10] = xgo * oil_dens * (oil_relp / oil_visco)         # gas component in oil phase
        values[11] = (1 - xgo) * oil_dens * (oil_relp / oil_visco)   # oil component in oil phase
        values[12] = 0                                               # wat component in oil phase

        # (2) oil phase
        values[13] = wat_dens              # oil density operator
        values[14] = 0                     # water and oil, pw = po - pcow
        values[15] = 0                                   # gas component in water phase
        values[16] = 0                                   # oil component in water phase
        values[17] = wat_dens * (wat_relp / wat_visco)   # wat component in water phase
        return 0

class black_oil_acc_flux_capillary_evaluator_python_well(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()
        self.property = property_container

    def evaluate(self, state, values):
        pbub = self.property.pb_ev.evaluate(state)
        Rs = self.property.rs_ev.evaluate(state, pbub)
        xgo = self.property.xgo_ev.evaluate(Rs)
        oil_dens = self.property.density_ev['oil'].evaluate(state, pbub, xgo)
        wat_dens = self.property.density_ev['water'].evaluate(state)
        gas_dens = self.property.density_ev['gas'].evaluate(state)
        oil_visco = self.property.viscosity_ev['oil'].evaluate(state, pbub)
        wat_visco = self.property.viscosity_ev['water'].evaluate(state)
        gas_visco = self.property.viscosity_ev['gas'].evaluate(state)
        wat_sat = self.property.sat_ev['water'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        oil_sat = self.property.sat_ev['oil'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        gas_sat = 1 - wat_sat - oil_sat
        wat_relp = self.property.rel_perm_ev['water'].evaluate(wat_sat)
        krow = self.property.rel_perm_ev['OW'].evaluate(wat_sat)
        gas_relp = self.property.rel_perm_ev['gas'].evaluate(gas_sat)
        krog = self.property.rel_perm_ev['GO'].evaluate(gas_sat)
        oil_relp = self.property.rel_perm_ev['oil'].evaluate(wat_sat, gas_sat, krow, krog)
        rock_cp = self.property.rock_compress_ev.evaluate(state)

        # accumulation operator
        values[0] = rock_cp * xgo * oil_dens * oil_sat + rock_cp * gas_sat * gas_dens
        values[1] = rock_cp * (1 - xgo) * oil_dens * oil_sat
        values[2] = rock_cp * wat_sat * wat_dens

        # flux operator
        # (1) gas phase
        values[3] = 0                    # gas density operator
        values[4] = 0                    # gas and oil, pg=po+pcgo=po-(-pcgo)
        values[5] = gas_dens * (gas_relp / gas_visco)   # gas component in gas phase
        values[6] = 0                                   # oil component in gas phase
        values[7] = 0                                   # water component in gas phase

        # (2) oil phase
        values[8] = 0                    # oil density operator
        values[9] = 0                    # reference phase, pc = 0
        values[10] = xgo * oil_dens * (oil_relp / oil_visco)         # gas component in oil phase
        values[11] = (1 - xgo) * oil_dens * (oil_relp / oil_visco)   # oil component in oil phase
        values[12] = 0                                               # wat component in oil phase

        # (2) oil phase
        values[13] = 0                     # oil density operator
        values[14] = 0                     # water and oil, pw = po - pcow
        values[15] = 0                                   # gas component in water phase
        values[16] = 0                                   # oil component in water phase
        values[17] = wat_dens * (wat_relp / wat_visco)   # wat component in water phase
        return 0

class black_oil_rate_evaluator_python(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()
        self.property = property_container

    def evaluate(self, state, values):
        pbub = self.property.pb_ev.evaluate(state)
        Rs = self.property.rs_ev.evaluate(state, pbub)
        xgo = self.property.xgo_ev.evaluate(Rs)
        oil_dens = self.property.density_ev['oil'].evaluate(state, pbub, xgo)
        wat_dens = self.property.density_ev['water'].evaluate(state)
        gas_dens = self.property.density_ev['gas'].evaluate(state)
        oil_visco = self.property.viscosity_ev['oil'].evaluate(state, pbub)
        wat_visco = self.property.viscosity_ev['water'].evaluate(state)
        gas_visco = self.property.viscosity_ev['gas'].evaluate(state)
        wat_sat = self.property.sat_ev['water'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        oil_sat = self.property.sat_ev['oil'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        gas_sat = 1 - wat_sat - oil_sat
        wat_relp = self.property.rel_perm_ev['water'].evaluate(wat_sat)
        krow = self.property.rel_perm_ev['OW'].evaluate(wat_sat)
        gas_relp = self.property.rel_perm_ev['gas'].evaluate(gas_sat)
        krog = self.property.rel_perm_ev['GO'].evaluate(gas_sat)
        oil_relp = self.property.rel_perm_ev['oil'].evaluate(wat_sat, gas_sat, krow, krog)

        # flux in reservoir condition
        gas_flux = xgo * oil_dens * (oil_relp / oil_visco) + gas_dens * (gas_relp / gas_visco)
        oil_flux = (1 - xgo) * oil_dens * (oil_relp / oil_visco)
        wat_flux = wat_dens * (wat_relp / wat_visco)

        # surface density
        self.surface_wat_dens = self.property.surf_wat_dens
        self.surface_oil_dens = self.property.surf_oil_dens
        self.surface_gas_dens = self.property.surf_gas_dens

        # convert to surface condition
        values[0] = gas_flux / self.surface_gas_dens
        values[1] = oil_flux/ self.surface_oil_dens
        values[2] = wat_flux/ self.surface_wat_dens

        return 0

class Properties(operator_set_evaluator_iface):
    def __init__(self, property_container):
        super().__init__()
        self.property = property_container

    def evaluate(self, state):
        pbub = self.property.pb_ev.evaluate(state)
        Rs = self.property.rs_ev.evaluate(state, pbub)
        xgo = self.property.xgo_ev.evaluate(Rs)
        oil_dens = self.property.density_ev['oil'].evaluate(state, pbub, xgo)
        wat_dens = self.property.density_ev['water'].evaluate(state)
        gas_dens = self.property.density_ev['gas'].evaluate(state)
        oil_visco = self.property.viscosity_ev['oil'].evaluate(state, pbub)
        wat_visco = self.property.viscosity_ev['water'].evaluate(state)
        gas_visco = self.property.viscosity_ev['gas'].evaluate(state)
        wat_sat = self.property.sat_ev['water'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        oil_sat = self.property.sat_ev['oil'].evaluate(state, pbub, gas_dens, oil_dens, wat_dens, xgo)
        gas_sat = 1 - wat_sat - oil_sat
        wat_relp = self.property.rel_perm_ev['water'].evaluate(wat_sat)
        krow = self.property.rel_perm_ev['OW'].evaluate(wat_sat)
        gas_relp = self.property.rel_perm_ev['gas'].evaluate(gas_sat)
        krog = self.property.rel_perm_ev['GO'].evaluate(gas_sat)
        oil_relp = self.property.rel_perm_ev['oil'].evaluate(wat_sat, gas_sat, krow, krog)
        rock_cp = self.property.rock_compress_ev.evaluate(state)


        return oil_sat, wat_sat, gas_sat, pbub