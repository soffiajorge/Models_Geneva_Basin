import numpy as np
from darts.physics import *
from interpolation import TableInterpolation
from darts.tools.keyword_file_tools import *

class PropertyContainer(property_evaluator_iface):
    def __init__(self, phase_name, component_name, pvt):
        super().__init__()
        # This class contains all the property evaluators required for simulation
        self.n_phases = len(phase_name)
        self.nc = len(component_name)
        self.component_name = component_name
        self.phase_name = phase_name
        self.pvt = pvt
        self.surf_dens = get_table_keyword(self.pvt, 'DENSITY')[0]
        self.surf_oil_dens = self.surf_dens[0]
        self.surf_wat_dens = self.surf_dens[1]
        self.surf_gas_dens = self.surf_dens[2]

        # Allocate (empty) evaluators
        self.pb_ev = []
        self.rs_ev = []
        self.xgo_ev = []
        self.density_ev = []
        self.viscosity_ev = []
        self.saturation_ev = []
        self.rel_perm_ev = []
        self.rock_compress_ev = []


class DensityGas(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.gas_dens = get_table_keyword(self.pvt, 'DENSITY')[0][2]
        self.table = get_table_keyword(self.pvt, 'PVDG')

    def evaluate(self, state):
        pres = state[0]
        pres_index = 0    # first column in the table - pressure
        Bg_index = 1      # second column in the table - volume factor
        Table = TableInterpolation()

        if (pres < self.table[0][0] or pres > self.table[len(self.table) - 1][0]):
            Bg = Table.LinearExtraP(self.table, pres, pres_index, Bg_index)
        else:
            Bg = Table.LinearInterP(self.table, pres, pres_index, Bg_index)

        return self.gas_dens / Bg


class DensityWater(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.wat_dens = get_table_keyword(self.pvt, 'DENSITY')[0][1]
        self.table = get_table_keyword(self.pvt, 'PVTW')[0]

    def evaluate(self, state):
        pres = state[0]
        X = self.table[2] * (pres - self.table[0])
        Bw = self.table[1] / (1 + X + X * X / 2)

        return self.wat_dens / Bw


class VsicoGas(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.table = get_table_keyword(self.pvt, 'PVDG')

    def evaluate(self, state):
        pres = state[0]
        pres_index = 0  # first column in the table - pressure
        vgas_index = 2   # third column in the table - viscosity
        Table = TableInterpolation()

        if (pres < self.table[0][0] or pres > self.table[len(self.table) - 1][0]):
            visco_gas = Table.LinearExtraP(self.table, pres, pres_index, vgas_index)
        else:
            visco_gas = Table.LinearInterP(self.table, pres, pres_index, vgas_index)

        return visco_gas


class ViscoWat(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.table = get_table_keyword(self.pvt, 'PVTW')[0]

    def evaluate(self, state):
        pres = state[0]
        Y = -self.table[4] * (pres - self.table[0])
        visco_wat = self.table[3] / (1 + Y + Y * Y / 2)

        return visco_wat


class BubblePointPres(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.dens = get_table_keyword(self.pvt, 'DENSITY')[0]
        self.oil_dens = self.dens[0]
        self.gas_dens = self.dens[2]
        self.table = get_table_keyword(self.pvt, 'PVTO')
        self.len_table = len(self.table)

    def evaluate(self, state):
        pres = state[0]
        zg = state[1]
        zo = state[2]
        sat_rs = self.oil_dens / (self.gas_dens * (zo+zg) / zg - self.gas_dens)
        rs_index = 0        # first column in the table - rs
        pres_index = 1      # second column in the table - pressure
        Table = TableInterpolation()

        # find the index of saturated Rs
        for i in range(self.len_table - 1):
            if (self.table[i][0] == self.table[i+1][0]):
                num = i
                break

        if (sat_rs < self.table[0][0]):
            pbub = Table.LinearExtraP(self.table, sat_rs, rs_index, pres_index)
        elif (sat_rs > self.table[self.len_table-1][0]):
            pbub = Table.SatExtrapolation(self.table, sat_rs, rs_index, pres_index, num)
        else:
            pbub = Table.LinearInterP(self.table, sat_rs, rs_index, pres_index)

        if pres>pbub:
            pbub = pbub
        else:
            pbub = pres

        return pbub


class Rs(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.dens = get_table_keyword(self.pvt, 'DENSITY')[0]
        self.oil_dens = self.dens[0]
        self.gas_dens = self.dens[2]
        self.table = get_table_keyword(self.pvt, 'PVTO')
        self.len_table = len(self.table)

    def evaluate(self, state, pbub):
        pres = state[0]
        zg = state[1]
        zo = state[2]
        rs_index = 0  # first column in the table - rs
        pres_index = 1    # second column in the table - pressure
        Table = TableInterpolation()

        # find the index of saturated Rs
        for i in range(self.len_table - 1):
            if (self.table[i][0] == self.table[i + 1][0]):
                num = i
                break

        # undersaturated condition
        if (pres>pbub):
            rs = self.oil_dens / (self.gas_dens * (zo+zg) / zg - self.gas_dens)
        # saturated condition
        else:
            if pres<self.table[0][1]:
                rs = Table.LinearExtraP(self.table, pres, pres_index, rs_index)
            elif pres>self.table[num][1]:
                rs = Table.SatExtrapolation(self.table, pres, pres_index, rs_index, num)
            else:
                rs = Table.LinearInterP(self.table, pres, pres_index, rs_index)

        return rs


class Xgo(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.dens = get_table_keyword(self.pvt, 'DENSITY')[0]
        self.oil_dens = self.dens[0]
        self.gas_dens = self.dens[2]
        self.table = get_table_keyword(self.pvt, 'PVTO')

    def evaluate(self, rs):

        xgo = self.gas_dens * rs / (self.oil_dens + self.gas_dens * rs)

        return xgo


class DensityOil(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.oil_dens = get_table_keyword(self.pvt, 'DENSITY')[0][0]
        self.table = get_table_keyword(self.pvt, 'PVTO')
        self.len_table = len(self.table)

    def evaluate(self, state, pbub, xgo):
        pres = state[0]
        pres_index = 1
        Bo_index = 2

        Table = TableInterpolation()

        # find the index of max Bo
        for i in range(self.len_table - 1):
            if (self.table[i][2] >= self.table[i + 1][2]):
                num = i
                break

        # calculate the saturated Bo
        if (pbub < self.table[0][1]):
            Bo_bub = Table.LinearExtraP(self.table, pbub, pres_index, Bo_index)
        elif (pbub > self.table[num][1]):
            Bo_bub = Table.SatExtrapolation(self.table, pbub, pres_index, Bo_index, num)
        else:
            Bo_bub = Table.LinearInterP(self.table, pbub, pres_index, Bo_index)

        # calculate Bo in current pressure
        # (1) saturated condition
        if pres < pbub:
            if (pres < self.table[0][1]):
                Bo = Table.LinearExtraP(self.table, pres, pres_index, Bo_index)
            elif (pres > self.table[num][1]):
                Bo = Table.SatExtrapolation(self.table, pres, pres_index, Bo_index, num)
            else:
                Bo = Table.LinearInterP(self.table, pres, pres_index, Bo_index)
        # (2) undersaturated condition
        else:
            pres_undersat = pres + self.table[num][1] - pbub
            Bo_undersat = Table.SatExtrapolation(self.table, pres_undersat, pres_index, Bo_index, num+1)
            if (pbub < self.table[num][1]):
                Bo = Bo_undersat * Bo_bub / self.table[num][2]
            else:
                Bo = Bo_undersat - (self.table[num][2] - Bo_bub)

        return self.oil_dens / Bo / (1 - xgo)


class ViscoOil(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.table = get_table_keyword(self.pvt, 'PVTO')
        self.len_table = len(self.table)

    def evaluate(self, state, pbub):
        pres = state[0]
        pres_index = 1
        visco_index = 3

        Table = TableInterpolation()

        # find the index of min visco
        for i in range(self.len_table - 1):
            if (self.table[i][3] <= self.table[i + 1][3]):
                num = i
                break

        # calculate the saturated viscosity
        if (pbub < self.table[0][1]):
            visco_bub = Table.LinearExtraP(self.table, pbub, pres_index, visco_index)
        elif (pbub > self.table[num][1]):
            visco_bub = Table.SatExtrapolation(self.table, pbub, pres_index, visco_index, num)
        else:
            visco_bub = Table.LinearInterP(self.table, pbub, pres_index, visco_index)

        # calculate viscosity in current pressure
        # (1) saturated condition
        if pres < pbub:
            if (pres < self.table[0][1]):
                visco = Table.LinearExtraP(self.table, pres, pres_index, visco_index)
            elif (pres > self.table[num][1]):
                visco = Table.SatExtrapolation(self.table, pres, pres_index, visco_index, num)
            else:
                visco = Table.LinearInterP(self.table, pres, pres_index, visco_index)
        # (2) undersaturated condition
        else:
            pres_undersat = pres + self.table[num][1] - pbub
            visco_undersat = Table.SatExtrapolation(self.table, pres_undersat, pres_index, visco_index, num + 1)
            if (pbub < self.table[num][1]):
                visco = visco_undersat * visco_bub / self.table[num][3]
            else:
                visco = visco_undersat - (self.table[num][3] - visco_bub)

        return visco


class Watsat(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt

    def evaluate(self, state, pbub, gas_dens, oil_dens, wat_dens, xgo):
        pres = state[0]
        zg = state[1]
        zo = state[2]

        Nu_w = 1 - zg - zo
        # (1) two phase undersaturated condition
        if (pres > pbub):
            Nu_g = 0
            Nu_o = zo
            wat_sat = oil_dens * Nu_w / (oil_dens * Nu_w + wat_dens * Nu_o)

        # (2) three phases saturated condition
        else:
            Nu_o = zo / (1 - xgo)
            Nu_g = 1 - Nu_o - Nu_w
            wat_sat = gas_dens * oil_dens * Nu_w / (gas_dens * oil_dens * Nu_w + gas_dens * wat_dens * Nu_o +
                                                    oil_dens * wat_dens * Nu_g)
        wat_sat = max(wat_sat, 0.0)
        wat_sat = min(wat_sat, 1.0)

        return wat_sat


class Oilsat(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt

    def evaluate(self, state, pbub, gas_dens, oil_dens, wat_dens, xgo):
        pres = state[0]
        zg = state[1]
        zo = state[2]

        Nu_w = 1 - zg - zo
        # (1) two phase undersaturated condition
        if (pres > pbub):
            Nu_g = 0
            Nu_o = zo
            oil_sat = wat_dens * Nu_o / (oil_dens * Nu_w + wat_dens * Nu_o)

        # (2) three phases saturated condition
        else:
            Nu_o = zo / (1 - xgo)
            Nu_g = 1 - Nu_o - Nu_w
            oil_sat = gas_dens * wat_dens * Nu_o / (gas_dens * oil_dens * Nu_w + gas_dens * wat_dens * Nu_o +
                                                    oil_dens * wat_dens * Nu_g)
        oil_sat = max(oil_sat, 0.0)
        oil_sat = min(oil_sat, 1.0)

        return oil_sat


class Gassat(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt

    def evaluate(self, oil_sat, wat_sat):
        gas_sat = 1 - oil_sat - wat_sat
        gas_sat = max(gas_sat, 0.0)
        gas_sat = min(gas_sat, 1.0)

        return gas_sat


class WatRelPerm(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.SWOF = get_table_keyword(self.pvt, 'SWOF')

    def evaluate(self, wat_sat):
        wat_index = 0
        krw_index = 1

        Table = TableInterpolation()
        if (wat_sat < self.SWOF[0][0] or wat_sat > self.SWOF[len(self.SWOF) - 1][0]):
            krw = Table.SCALExtraP(self.SWOF, wat_sat, wat_index, krw_index)
        else:
            krw = Table.LinearInterP(self.SWOF, wat_sat, wat_index, krw_index)

        return krw


class Krow(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.SWOF = get_table_keyword(self.pvt, 'SWOF')

    def evaluate(self, wat_sat):
        wat_index = 0
        krow_index = 2

        Table = TableInterpolation()
        if (wat_sat < self.SWOF[0][0] or wat_sat > self.SWOF[len(self.SWOF) - 1][0]):
            krow = Table.SCALExtraP(self.SWOF, wat_sat, wat_index, krow_index)
        else:
            krow = Table.LinearInterP(self.SWOF, wat_sat, wat_index, krow_index)

        return krow


class GasRelPerm(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.SGOF = get_table_keyword(self.pvt, 'SGOF')

    def evaluate(self, gas_sat):
        gas_index = 0
        krg_index = 1

        Table = TableInterpolation()
        if (gas_sat < self.SGOF[0][0] or gas_sat > self.SGOF[len(self.SGOF) - 1][0]):
            krg = Table.SCALExtraP(self.SGOF, gas_sat, gas_index, krg_index)
        else:
            krg = Table.LinearInterP(self.SGOF, gas_sat, gas_index, krg_index)

        return krg


class Krog(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.SGOF = get_table_keyword(self.pvt, 'SGOF')

    def evaluate(self, gas_sat):
        gas_index = 0
        krog_index = 2

        Table = TableInterpolation()
        if (gas_sat < self.SGOF[0][0] or gas_sat > self.SGOF[len(self.SGOF) - 1][0]):
            krog = Table.SCALExtraP(self.SGOF, gas_sat, gas_index, krog_index)
        else:
            krog = Table.LinearInterP(self.SGOF, gas_sat, gas_index, krog_index)

        return krog


# here we use Stone I model to calculate oil relperm
class OilRelPerm(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.SGOF = get_table_keyword(self.pvt, 'SGOF')
        self.SWOF = get_table_keyword(self.pvt, 'SWOF')

    def evaluate(self, wat_sat, gas_sat, krow, krog):
        len_SWOF = len(self.SWOF)
        len_SGOF = len(self.SGOF)
        Sorw = 1 - self.SWOF[len_SWOF - 1][0]
        Sorg = 1 - self.SGOF[len_SGOF - 1][0]
        Swc = self.SWOF[0][0]
        MINIMAL_FOR_COMPARE = 1e-12
        Krocw = self.SWOF[0][2]
        SatLimit = 1 - 2 * MINIMAL_FOR_COMPARE
        oil_sat = 1 - wat_sat - gas_sat

        if (gas_sat < MINIMAL_FOR_COMPARE):
            kro = krow        # water-oil two phases
        elif (wat_sat < Swc):
            kro = krog        # water phase not mobile
        else:
            # Stone I model -> alpha
            alpha = 1 - gas_sat / (1 - min(Swc + Sorg, SatLimit))
            # -> Som
            Som = alpha * Sorw + (1.0 - alpha) * Sorg
            # -> denom
            denom = 1.0 / (1.0 - min(Swc + Som, SatLimit))
            # Normalized saturations
            if (oil_sat - Som) > 0:
                Ma = oil_sat - Som
            else:
                Ma = 0
            SoStar = Ma * denom
            SwStar = min(wat_sat - Swc, SatLimit) * denom
            SgStar = min(gas_sat, SatLimit) * denom

            kro = (SoStar * krow * krog) / (Krocw * (1.0 - SwStar) * (1.0 - SgStar))

        return kro


class RockCompactionEvaluator(property_evaluator_iface):
    def __init__(self, pvt):
        super().__init__()
        self.pvt = pvt
        self.rock_table = get_table_keyword(self.pvt, 'ROCK')

    def evaluate(self, state):
        pressure = state[0]
        pressure_ref = self.rock_table[0][0]
        compressibility = self.rock_table[0][1]

        return (1.0 + compressibility * (pressure - pressure_ref))