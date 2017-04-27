import math
import numpy as np 

def tokenval(fname, token):
	file_in = open(fname, 'r')
	lines = file_in.readlines()
	file_in.close()
	found = 0
	for line in lines:
		if (token in line):
			strnext = line.split(token)[1]
			val = strnext.split()[0]
			valf = float(val)
			found = 1
			break

	if (found == 0):
		print "could not find token : ", token , " in file : ", fname
	return valf



def findconstants(filelist):
	hlist = []
	errlist = []
	constants = []
	
	for fname in  filelist:
		token = 'hmin : '
		val_hmin = tokenval(fname, token)
		token = 'hmax : '
		val_hmax = tokenval(fname, token)
		if val_hmin <> val_hmax :
			print "non uniform mesh used for checking convergence for log : ", fname

		token = '0.49) : '
		val_err = tokenval(fname, token)
		hlist.append(val_hmin)
		errlist.append(val_err)
	
	
 
	for i in range(len(hlist)-1):
		constants.append(math.fabs(math.log( errlist[i+1] / errlist[i]  ) / math.log(hlist[i+1] / hlist[i] )) )
		

	return hlist, errlist, constants


if __name__ == '__main__': 
	#2D 1 process tests
	file = open("convergence_test_details.txt","w")
	
	hlist, errlist, constants_2D_1 = findconstants(['log_convergence_2D_numprocs_1_lvl_1', 'log_convergence_2D_numprocs_1_lvl_2', 'log_convergence_2D_numprocs_1_lvl_3', 'log_convergence_2D_numprocs_1_lvl_4', 'log_convergence_2D_numprocs_1_lvl_5'])
	if ( np.mean(np.asarray(constants_2D_1))  > 1.99):
		print "2D 1 prcess test : PASSED" 
	else: 
		print "2D 1 prcess test : FAILED" 
		
	file.write('2D mesh 1 process case :\n')
	file.write('mesh sizes :\n')
	file.write('----------------------\n')
	file.write(str(hlist))
	file.write('\n')
	file.write('errors     :\n')
	file.write('----------------------\n')
	file.write(str(errlist))
	file.write('\n')
	file.write('ln(e-1 / e) / ln(h-1 / h)     :\n')
	file.write('----------------------\n')
	file.write(str(constants_2D_1))
	file.write('\n')
	


	#2D 5 process tests

	hlist, errlist, constants_2D_5 = findconstants(['log_convergence_2D_numprocs_5_lvl_1', 'log_convergence_2D_numprocs_5_lvl_2', 'log_convergence_2D_numprocs_5_lvl_3', 'log_convergence_2D_numprocs_5_lvl_4', 'log_convergence_2D_numprocs_5_lvl_5'])
	if ( np.mean(np.asarray(constants_2D_5)) > 1.99):
		print "2D 5 prcess test : PASSED" 
	else: 
		print "2D 5 prcess test : FAILED" 


	file.write('2D mesh 5 process case :\n')
	file.write('mesh sizes :\n')
	file.write('----------------------\n')
	file.write(str(hlist))
	file.write('\n')
	file.write('errors     :\n')
	file.write('----------------------\n')
	file.write(str(errlist))
	file.write('\n')
	file.write('ln(e-1 / e) / ln(h-1 / h)     :\n')
	file.write('----------------------\n')
	file.write(str(constants_2D_5))
	file.write('\n')
	

	#print constants_2D_5


	#3D different num process tests
	hlist, errlist, constants_3D_5 = findconstants(['log_convergence_3D_numprocs_15_lvl_1', 'log_convergence_3D_numprocs_20_lvl_2', 'log_convergence_3D_numprocs_30_lvl_3', 'log_convergence_3D_numprocs_32_lvl_4', 'log_convergence_3D_numprocs_32_lvl_5'])
	if ( np.mean(np.asarray(constants_3D_5)) > 1.99):
		print "2D 5 prcess test : PASSED" 
	else: 
		print "2D 5 prcess test : FAILED" 

	#print constants_3D_5
	file.write('3D mesh case :\n')
	file.write('mesh sizes :\n')
	file.write('----------------------\n')
	file.write(str(hlist))
	file.write('\n')
	file.write('errors     :\n')
	file.write('----------------------\n')
	file.write(str(errlist))
	file.write('\n')
	file.write('ln(e-1 / e) / ln(h-1 / h)     :\n')
	file.write('----------------------\n')
	file.write(str(constants_3D_5))
	file.write('\n')
	

	file.close();
