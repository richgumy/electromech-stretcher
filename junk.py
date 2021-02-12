"""
Just absolute garbage.py
"""
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

step_profile = [1,2,3]
filename = "test this juan"

log_file = open("log.txt","a")
log_file.write(str(datetime.date(datetime.now()))+'\n')
log_file.write(str(datetime.time(datetime.now()))+'\n')
log_file.write('filename='+filename+'\n')
log_file.write('step profile='+str(step_profile)+'\n')
log_file.write('velocity profile='+str(velocity_profile)+'\n')
log_file.write('diff_min(convergence checker)='+str(diff_min)+'\n')
log_file.write('iter_max(convergence timeout)='+str(iter_max)+'\n')
log_file.close()


# def main(input_filename):
#     if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input
#     print(input_filename)
#     Rout = np.array([])
#     Rin = np.array([])
#     tR = np.array([])
#     P = np.array([])
#     tP = np.array([])
#     F = np.array([])
#     tF = np.array([])
#
#     with open(input_filename, 'r', newline='') as csvfile:
#         data = csv.reader(csvfile, delimiter=',')
#         line_count = 0
#         for row in data:
#             if line_count == 0:
#                 # print(f'Column names are {", ".join(row)}')
#                 line_count = 1
#             else:
#                 R = row[0].strip('][').split(', ')
#                 Rout = np.append(Rout,float(R[0]))
#                 Rin = np.append(Rin,float(R[1]))
#                 tR = np.append(tR,float(row[1]))
#                 P = np.append(P,float(row[2]))
#                 tP = np.append(tP,float(row[3]))
#                 F = np.append(F,float(row[4]))
#                 tF = np.append(tF,float(row[5]))
#
#     plt.plot(tR, Rout)
#     # plt.plot(Rin, tR)
#     plt.show()
#
# # The real main driver
# if __name__ == "__main__":
# 	# Default input parameters
# 	input_filename = ""
# 	if len(sys.argv)>1: input_filename   = (sys.argv[1])
# 	# if len(sys.argv)>2: input2   =int(sys.argv[2])
#
# 	main(input_filename)
