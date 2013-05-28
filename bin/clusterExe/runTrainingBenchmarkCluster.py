import pcssTools
import pcssCluster
import configobj

import sys
import os


configFileName = sys.argv[1]

pcssConfig = configobj.ConfigObj(configFileName)
#pcssConfig["attribute_file_name"] = os.path.join(pcssConfig["pcss_directory"], "data", "context", "svmTrainingAttributes.txt")
runner = pcssTools.TrainingBenchmarkRunner(pcssConfig)
runner.execute()

runner = pcssTools.LeaveOneOutBenchmarkRunner(pcssConfig)
runner.execute()

runner = pcssTools.CompleteSvmRunner(pcssConfig)
runner.execute()
