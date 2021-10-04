import pyemma 
import numpy as np 
import matplotlib.pyplot as plt 
'''

We can use this method when the fraction of states/counts used is not equlal to 1 

'''
M = pyemma.load('M.h5')

dtraj = M.dtrajs_active
dtraj_concat = np.concatenate(dtraj)
print(np.shape(dtraj))
print(np.shape(dtraj_concat))
np.savetxt('dtraj.txt',dtraj_concat,'%6.0d')


j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] > 0.0):   
        j1=j1+1   
 
    if(M.eigenvectors_right()[i,1] < 0.0):   
        j2=j2+1   
print(j1,j2)
Mp = np.zeros((j1,4))   
Mn = np.zeros((j2,4))  
j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] > 0.0):   
        Mp[j1][0] = i+1   
        Mp[j1][1] = i   
        Mp[j1][2] = M.eigenvectors_right()[i,1]   
        Mp[j1][3] = M.eigenvectors_left()[0,i]   
        j1=j1+1  
    if(M.eigenvectors_right()[i,1] < 0.0):    
        Mn[j2][0] = i+1    
        Mn[j2][1] = i    
        Mn[j2][2] = M.eigenvectors_right()[i,1]    
        Mn[j2][3] = M.eigenvectors_left()[0,i]    
        j2=j2+1    
print("p:  ",np.sum(Mp[:,3])*100 )   
print("n:  ",np.sum(Mn[:,3])*100 )  
##############################################################################################  
j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0):   
        j1=j1+1   
      
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0):   
        j2=j2+1   
print(j1,j2)
Mpp = np.zeros((j1,4))   
Mpn = np.zeros((j2,4))  
j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0):    
        Mpp[j1][0] = i+1   
        Mpp[j1][1] = i   
        Mpp[j1][2] = M.eigenvectors_right()[i,1]   
        Mpp[j1][3] = M.eigenvectors_left()[0,i]   
        j1=j1+1  
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0):     
        Mpn[j2][0] = i+1    
        Mpn[j2][1] = i    
        Mpn[j2][2] = M.eigenvectors_right()[i,1]    
        Mpn[j2][3] = M.eigenvectors_left()[0,i]    
        j2=j2+1    
print("pp:  ",np.sum(Mpp[:,3])*100 )   
print("pn:  ",np.sum(Mpn[:,3])*100 )   
 
##############################################################################################  

j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0):   
        j1=j1+1   
      
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0):   
        j2=j2+1   
 
print(j1,j2)         
Mnp = np.zeros((j1,4))   
Mnn = np.zeros((j2,4))  
  
j1=0;j2=0   
for i in range(np.shape(M.eigenvectors_right())[0]):   
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0):    
        Mnp[j1][0] = i+1   
        Mnp[j1][1] = i   
        Mnp[j1][2] = M.eigenvectors_right()[i,1]   
        Mnp[j1][3] = M.eigenvectors_left()[0,i]   
        j1=j1+1  
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0):     
        Mnn[j2][0] = i+1    
        Mnn[j2][1] = i    
        Mnn[j2][2] = M.eigenvectors_right()[i,1]    
        Mnn[j2][3] = M.eigenvectors_left()[0,i]    
        j2=j2+1    
  
  
   
print("np:  ",np.sum(Mnp[:,3])*100 )   
print("nn:  ",np.sum(Mnn[:,3])*100 )   

##############################################################################################
np.savetxt('pp.txt',Mpp[:,1],'%6.0d')
np.savetxt('pn.txt',Mpn[:,1],'%6.0d')
np.savetxt('np.txt',Mnp[:,1],'%6.0d')
np.savetxt('nn.txt',Mnn[:,1],'%6.0d')
##############################################################################################
j1=0;j2=0    
for i in range(np.shape(M.eigenvectors_right())[0]):    
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] > 0.0):    
        j1=j1+1    
       
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] < 0.0):    
        j2=j2+1    
print(j1,j2)          
Mpnp = np.zeros((j1,4))    
Mpnn = np.zeros((j2,4))   
   
j1=0;j2=0    
for i in range(np.shape(M.eigenvectors_right())[0]):    
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] > 0.0):      
        Mpnp[j1][0] = i+1    
        Mpnp[j1][1] = i    
        Mpnp[j1][2] = M.eigenvectors_right()[i,1]    
        Mpnp[j1][3] = M.eigenvectors_left()[0,i]    
        j1=j1+1   
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] < 0.0):   
        Mpnn[j2][0] = i+1     
        Mpnn[j2][1] = i     
        Mpnn[j2][2] = M.eigenvectors_right()[i,1]     
        Mpnn[j2][3] = M.eigenvectors_left()[0,i]     
        j2=j2+1     
   
   
    
print("pnp:  ",np.sum(Mpnp[:,3])*100 )    
print("pnn:  ",np.sum(Mpnn[:,3])*100 )    

##############################################################################################
j1=0;j2=0    
for i in range(np.shape(M.eigenvectors_right())[0]):    
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0):    
        j1=j1+1    
       
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0):    
        j2=j2+1    
print(j1,j2) 
Mppp = np.zeros((j1,4))    
Mppn = np.zeros((j2,4))   
   
j1=0;j2=0    
for i in range(np.shape(M.eigenvectors_right())[0]):    
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0):      
        Mppp[j1][0] = i+1    
        Mppp[j1][1] = i    
        Mppp[j1][2] = M.eigenvectors_right()[i,1]    
        Mppp[j1][3] = M.eigenvectors_left()[0,i]    
        j1=j1+1   
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0):   
        Mppn[j2][0] = i+1     
        Mppn[j2][1] = i     
        Mppn[j2][2] = M.eigenvectors_right()[i,1]     
        Mppn[j2][3] = M.eigenvectors_left()[0,i]     
        j2=j2+1     
   
   
    
print("ppp:  ",np.sum(Mppp[:,3])*100 )    
print("ppn:  ",np.sum(Mppn[:,3])*100 )    
  
##############################################################################################
np.savetxt('pnp.txt',Mpnp[:,1],'%6.0d')
np.savetxt('pnn.txt',Mpnn[:,1],'%6.0d')
np.savetxt('ppp.txt',Mppp[:,1],'%6.0d')
np.savetxt('ppn.txt',Mppn[:,1],'%6.0d')
##############################################################################################
j1=0;j2=0 
for i in range(np.shape(M.eigenvectors_right())[0]): 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0): 
        j1=j1+1 
 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0): 
        j2=j2+1 
print(j1,j2) 
Mnpp = np.zeros((j1,4)) 
Mnpn = np.zeros((j2,4)) 
 
j1=0;j2=0 
for i in range(np.shape(M.eigenvectors_right())[0]): 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0): 
        Mnpp[j1][0] = i+1 
        Mnpp[j1][1] = i 
        Mnpp[j1][2] = M.eigenvectors_right()[i,1] 
        Mnpp[j1][3] = M.eigenvectors_left()[0,i] 
        j1=j1+1 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0): 
        Mnpn[j2][0] = i+1 
        Mnpn[j2][1] = i 
        Mnpn[j2][2] = M.eigenvectors_right()[i,1] 
        Mnpn[j2][3] = M.eigenvectors_left()[0,i] 
        j2=j2+1 
 
 
 
print("npp:  ",np.sum(Mnpp[:,3])*100 ) 
print("npn:  ",np.sum(Mnpn[:,3])*100 )                                                                                      

##############################################################################################
j1=0;j2=0 
for i in range(np.shape(M.eigenvectors_right())[0]): 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] > 0.0): 
        j1=j1+1 
 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] < 0.0): 
        j2=j2+1 
print(j1,j2) 
Mnnp = np.zeros((j1,4)) 
Mnnn = np.zeros((j2,4)) 
 
j1=0;j2=0 
for i in range(np.shape(M.eigenvectors_right())[0]): 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] > 0.0): 
        Mnnp[j1][0] = i+1 
        Mnnp[j1][1] = i 
        Mnnp[j1][2] = M.eigenvectors_right()[i,1] 
        Mnnp[j1][3] = M.eigenvectors_left()[0,i] 
        j1=j1+1 
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0) and (M.eigenvectors_right()[i,3] < 0.0): 
        Mnnn[j2][0] = i+1 
        Mnnn[j2][1] = i 
        Mnnn[j2][2] = M.eigenvectors_right()[i,1] 
        Mnnn[j2][3] = M.eigenvectors_left()[0,i] 
        j2=j2+1 
 
 
 
print("nnp:  ",np.sum(Mnnp[:,3])*100 ) 
print("nnn:  ",np.sum(Mnnn[:,3])*100 )                                                                                         

##############################################################################################
np.savetxt('nnp.txt',Mnnp[:,1],'%6.0d')
np.savetxt('nnn.txt',Mnnn[:,1],'%6.0d')
np.savetxt('npp.txt',Mnpp[:,1],'%6.0d')
np.savetxt('npn.txt',Mnpn[:,1],'%6.0d')
##############################################################################################

