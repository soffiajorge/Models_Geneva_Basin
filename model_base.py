from darts.models.reservoirs.struct_reservoir import StructReservoir
from darts.models.physics.dead_oil import DeadOil
from darts.models.darts_model import DartsModel
from darts.engines import value_vector, sim_params
import numpy as np
from darts.tools.keyword_file_tools import load_single_keyword, get_table_keyword
import os

class BaseModel(DartsModel):
    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()
        self.n_points = n_points
        # measure time spend on reading/initialization
        self.timer.node["initialization"].start()

        # create reservoir 
        filename = 'grid.grdecl'
        self.permx = load_single_keyword(filename, 'PERMX')
        self.permy = load_single_keyword(filename, 'PERMY')
        self.permz = load_single_keyword(filename, 'PERMZ')
        self.poro = load_single_keyword(filename, 'PORO')
        self.depth = load_single_keyword(filename, 'DEPTH')

        self.dx = load_single_keyword(filename, 'DX')
        self.dy = load_single_keyword(filename, 'DY')
        self.dz = load_single_keyword(filename, 'DZMTRXV')

        # Import other properties from files
        self.actnum = load_single_keyword(filename, 'ACTNUM')
        self.coord = load_single_keyword(filename, 'COORD')
        self.zcorn = load_single_keyword(filename, 'ZCORN')

        is_CPG = False  # True for re-calculation of dx, dy and dz from CPG grid

        self.reservoir = StructReservoir(self.timer, nx=68, ny=70, nz=55, dx=self.dx, dy=self.dy, dz=self.dz,
                                         permx=self.permx, permy=self.permy, permz=self.permz, poro=self.poro,
                                         depth=self.depth, actnum=self.actnum, coord=self.coord, zcorn=self.zcorn)

        if is_CPG:
            dx, dy, dz = self.reservoir.get_cell_cpg_widths()
            get_table_keyword('width.in', ['DX', 'DY', 'DZ'], [dx, dy, dz])

        well_dia = 0.152
        well_rad = well_dia / 2
        # """producers"""
        self.reservoir.add_well("PROD1", wellbore_diameter=well_dia)
        for i in range(53):
            #if (i + 1) #not in [4, 5, 14, 20]:
            self.reservoir.add_perforation(self.reservoir.wells[-1], 34, 15, i + 1, well_radius=well_rad, multi_segment=False)

        self.reservoir.add_well("PROD2", wellbore_diameter=well_dia)
        for i in range(53):
            #if (i + 1) #not in [4, 5, 9, 13, 14, 15, 16, 17, 18, 19, 20]:
            self.reservoir.add_perforation(self.reservoir.wells[-1], 34, 45, i + 1, well_radius=well_rad, multi_segment=False)

        self.reservoir.add_well("INJ", wellbore_diameter=well_dia)
        for i in range(53):
            #if (i + 1) #not in [4, 5, 9, 14, 20]:
            self.reservoir.add_perforation(self.reservoir.wells[-1], 34, 35, i + 1, well_radius=well_rad, multi_segment=False)

        self.timer.node["initialization"].stop()


    def set_op_list(self):
        self.op_num = np.array(self.reservoir.mesh.op_num, copy=False)
        n_res = self.reservoir.mesh.n_res_blocks
        self.op_num[n_res:] = 1
        self.op_list = [self.physics.acc_flux_itor, self.physics.acc_flux_itor_well]
