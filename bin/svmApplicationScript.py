import pcssTools
import pcssCluster
import configobj

import sys
import os


configFileName = sys.argv[1]
tempPcssConfig = configobj.ConfigObj(configFileName)
configSpec = tempPcssConfig["config_spec_file"]

pcssConfig = configobj.ConfigObj(configFileName, configspec=configSpec)
pcssConfig["attribute_file_name"] = os.path.join(pcssConfig["pcss_directory"], "data", "context", "svmApplicationFileAttributes.txt")
runner = pcssTools.SvmApplicationFeatureRunner(pcssConfig)
runner.execute()
