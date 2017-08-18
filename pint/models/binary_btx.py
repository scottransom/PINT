
"""This model provides the BTX model (D. Nice, unpublished).
    """
from pint import ls, GMsun, Tsun
from .stand_alone_psr_binaries.BT_model import BTmodel
from .pulsar_binary import PulsarBinary
from . import parameter as p
from .timing_model import TimingModel, MissingParameter
import astropy.units as u

class BinaryBTX(PulsarBinary):
    """This is a PINT pulsar binary BTX model class, a subclass of PulsarBinary.
    It is a wrapper for stand alone BTXmodel class defined in
    ./stand_alone_psr_binary/BTX_model.py
    All the detailed calculations are in the stand alone BTXmodel.
    The aim for this class is to connect the stand alone binary
    model with PINT platform.  Currently only the binary frequency
    derivatives, and not the jumps, are implemented.
    BTXmodel special parameters:
      FB, the orbital frequency
      FBx, where x is x'th derivative of the orbital frequency
      XDOT_x, where x is the x+1'th derivative of X (i.e. asini/c)
    """
    register = True
    def __init__(self):
        super(BinaryBTX, self).__init__()
        self.binary_model_name = 'BTX'
        self.binary_model_class = BTXmodel
        self.add_param(p.floatParameter(name="FB", value=0.0,
             units="Hz",
             description="Orbital Frequency"))
        self.add_param(p.prefixParameter(name="FB1",
                       value=0.0, units='Hz/s^1',
                       description="Orbital Frequency Derivative",
                       unit_template=self.F_unit,
                       description_template=self.F_description,
                       type_match='float'))

    def setup(self):
        super(BinaryBTX, self).setup()

        # Check if FB terms are in order.
        FB_mapping = self.get_prefix_mapping_component('FB')
        FB_terms = list(FB_mapping.keys())
        FB_terms.sort()
        FB_in_order = list(range(1, max(FB_terms)+1))
        if not FB_terms == FB_in_order:
            diff = list(set(FB_in_order) - set(FB_terms))
            raise MissingParameter("FB", "FB%d" % diff[0])
        self.num_FB_terms = len(FB_terms)

        # set up derivative functions
        for ii, val in FB_mapping.items():
            self.register_deriv_funcs(self.d_delay_FB_d_FBX, val)

        # If any necessary parameter is missing, raise MissingParameter.
        for p in ("FB", "T0", "A1"):
            if getattr(self, p).value is None:
                raise MissingParameter("BTX", p,
                                       "%s is required for BTX" % p)

        # If any *DOT is set, we need T0
        for p in ("FB1", "OMDOT", "EDOT", "XDOT"):
            if getattr(self, p).value is None:
                getattr(self, p).set("0")
                getattr(self, p).frozen = True

            if getattr(self, p).value is not None:
                if self.T0.value is None:
                    raise MissingParameter("BTX", "T0",
                        "T0 is required if FB1 or *DOT is set")

        if self.GAMMA.value is None:
            self.GAMMA.set("0")
            self.GAMMA.frozen = True
