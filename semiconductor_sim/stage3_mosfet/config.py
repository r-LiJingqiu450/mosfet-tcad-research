"""Single production configuration shared by study, freeze and regression tests."""
from dataclasses import asdict, dataclass
from common import constants as C

@dataclass(frozen=True)
class TableConfig:
    nps: int = 400
    npn: int = 220
    interpolation: str = 'log'
    phi_s_min: float = -.3
    phi_s_max: float = 1.6
    phi_n_max: float = 2.
    x_max: float = 1.2e-6
    d0: float = .125e-9
    hmax: float = 4e-9

    def to_dict(self):
        return asdict(self)

PRODUCTION = TableConfig()


def production_metadata():
    return dict(table=PRODUCTION.to_dict(), physics=dict(
        na_m3=C.NA, ni_m3=C.NI, temperature_K=C.TEMP,
        mu_m2_Vs=C.MU_N, tox_m=C.EPS_OX/C.COX, length_m=C.L_CH,
        width_m=C.W_CH, vfb_V=C.VFB, empirical_clm=False,
        density_clipping=False), transfer=dict(vg_min=0.,vg_max=2.,vg_nodes=81,
        vd_linear=.05,vd_high=1.5,ith_A=1e-7))
