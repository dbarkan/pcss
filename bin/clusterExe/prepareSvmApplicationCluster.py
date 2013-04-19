import sys
import os
import pcssTools
import pcssCluster
import unittest
import configobj


configFileName = sys.argv[1]
tempPcssConfig = configobj.ConfigObj(configFileName)
configSpecFile = tempPcssConfig["user_config_spec_file"]

pcssConfig = configobj.ConfigObj(configFileName, configspec=configSpecFile)
runner = pcssTools.SvmApplicationClusterRunner(pcssConfig)

runner.execute()
