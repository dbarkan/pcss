import pcssTools
import pcssCluster
import configobj

import sys
import os


configFileName = sys.argv[1]
tempPcssConfig = configobj.ConfigObj(configFileName)
configSpec = tempPcssConfig["user_config_spec_file"]

pcssConfig = configobj.ConfigObj(configFileName, configspec=configSpec)

runner = pcssTools.TrainingAnnotationRunner(pcssConfig)
runner.execute()
