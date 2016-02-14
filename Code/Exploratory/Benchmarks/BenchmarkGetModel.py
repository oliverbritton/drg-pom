# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 16:31:06 2016

@author: Oliver Britton
"""

import time
import PopulationOfModels as pom

start = time.time()

model = pom.GetModel('DeterminedOcelot')

end = time.time()

print end-start