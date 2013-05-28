import pcssTools
import pcssCluster
import configobj

import sys
import os


configFileName = sys.argv[1]

pcssConfig = configobj.ConfigObj(configFileName)

runner = pcssTools.TrainingAnnotationRunner(pcssConfig)
runner.execute()
