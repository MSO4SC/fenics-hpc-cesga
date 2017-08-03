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

                token = 'xxxx time step :0.998000, l2 error :'
                val_err = tokenval(fname, token)
                hlist.append(val_hmin)
                errlist.append(val_err)
    
    
 
        for i in range(len(hlist)-1):
                constants.append(math.fabs(math.log( errlist[i+1] / errlist[i]  ) / math.log(hlist[i+1] / hlist[i] )) )
    

        return hlist, errlist, constants


if __name__ == '__main__': 
        #2D tests
        file = open("convergence_test_details.txt","w")
    
        hlist, errlist, constants_2D_1 = findconstants(['logx1_r0', 'logx1_r1', 'logx1_r2'])
        if ( np.mean(np.asarray(constants_2D_1))  > 1.2):
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
    


        file.close();

