import sys
import os
import configobj
import pcssTools
import logging
#config file

configFileName = sys.argv[1]
tempPcssConfig = configobj.ConfigObj(configFileName)
configSpec = tempPcssConfig["disopred_standalone_config_spec_file"]

pcssConfig = configobj.ConfigObj(configFileName, configspec=configSpec)

runner = pcssTools.DisopredStandaloneRunner(pcssConfig)
logging.basicConfig(filename=runner.pdh.getFullOutputFile("%s_clusters.log" % pcssConfig["run_name"]), level=logging.DEBUG,
                    filemode="w", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

runner.execute()


