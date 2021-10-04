import numpy as np

dtrajs_1d=np.loadtxt('./dtrajs_1d.txt')

dtrajs_2d=np.load('../dtrajs.npy',allow_pickle=True)
dtrajs_2d=list(dtrajs_2d)

start = 0
dtrajs_list = []
for i in range(len(dtrajs_2d)): #for number of trajactories
    end = start+len(dtrajs_2d[i])+0 #getting end frame number
    temp = dtrajs_1d[start:end:1].reshape(len(dtrajs_2d[i]))+0
    dtrajs_list.append(np.array(temp,dtype=np.int32))
    start = end+0


np.save('./dtrajs_2d',dtrajs_list)
