import numpy as np
import pandas as pd

class ECC:

    def __init__(self):

        # print('Please input g(0) (like 101、111): ')
        # self.g0 = np.array([int(s) for s in input() ])
        # print('Please input g(1) (like 101、111): ')
        # self.g1 = np.array([int(s) for s in input() ])
        # print('Please input series (Binary): ')
        # self.u = np.array([int(s) for s in input() ])
        #self.StructBlock(self.g0,self.g1,self.u)
        g0 = np.array([1,0,1])
        g1 = np.array([1,1,0])
        u = np.array([1,1,0,1,0,0])
        self.Encode(g0,g1,u)
        
        pass



    def Encode(self,g0,g1,u):
        v = []
        G = np.zeros([len(u),( len(u)+len(g0)-1) * 2 ])
        g = []

        for i in range(len(g0)):
            g.append(g0[i])
            g.append(g1[i]) 

        print('g = ', g)

        print('u=\n',u)
        
        for i in range(len(G)):
            G[i][2*i:2*i+2*len(g1)] = g    
        print('G =\n',G)



        v = np.matmul(u,G)
        v = v.astype(np.int8)
        v = np.mod(v,2)    
        print('v=\n', v)
        return v

    def AWGNPass(self):
        R = []      



        return R

if __name__ == '__main__':
    ECC = ECC()

    

