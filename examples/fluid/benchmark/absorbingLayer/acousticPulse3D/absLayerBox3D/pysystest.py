
__pysys_title__   = r""" Absorbing layer 3D pulse, box """ 
#                        ==========================
__pysys_purpose__ = r""" Testing the box absorbing layer for a 3D pulse """ 
    
__pysys_created__ = "2025-05-06"
#__pysys_skipped_reason__   = "Skipped until Bug-1234 is fixed"

#__pysys_traceability_ids__ = "Bug-1234, UserStory-456" 
#__pysys_groups__           = "myGroup, disableCoverage, performance"
#__pysys_modes__            = lambda helper: helper.inheritedModes + [ {'mode':'MyMode', 'myModeParam':123}, ]
#__pysys_parameterized_test_modes__ = {'MyParameterizedSubtestModeA':{'myModeParam':123}, 'MyParameterizedSubtestModeB':{'myModeParam':456}, }

import pysys.basetest, pysys.mappers
from pysys.constants import *

from apes.apeshelper import ApesHelper
class PySysTest(ApesHelper, pysys.basetest.BaseTest):
    def setup(self):
        self.apes.setupMusubi()

    def execute(self):
        musrun = self.apes.runMusubi(np = 8)

    def validate(self):
        self.apes.checkMusLog()
        self.assertPathExists('tracking/pulse3D_probe_p00000.res',
                              abortOnError = True)
        self.apes.assertIsClose('pulse3D_probe_p00000.res',
                                dir = 'tracking')

