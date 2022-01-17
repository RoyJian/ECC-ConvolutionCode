import numpy as np
from numpy.random.mtrand import rand
import networkx as nx
import pylab
import matplotlib.pyplot as plt
import scipy as sp

class ECC:  

    def __init__(self):

        # print('Please input g(0) (like 101、111): ')
        # self.g0 = np.array([int(s) for s in input() ])
        # print('Please input g(1) (like 101、111): ')
        # self.g1 = np.array([int(s) for s in input() ])
        # print('Please input series (Binary): ')
        # self.u = np.array([int(s) for s in input() ])
        #self.StructBlock(self.g0,self.g1,self.u)
        self.g0 = np.array([1,1,1])
        self.g1 = np.array([1,0,1])
        self.u = np.array([1,1,0,1,1])
        v = self.Encode(self.g0,self.g1,self.u)
        r = self.AWGNPass(v,1)
        print('r = \n',r)
        
        pass



    def Encode(self,g0,g1,u):
        v = []
        G = np.zeros([len(u),(len(u)+len(g0)-1) * 2 ])
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
        # v = v[:][:2*len(u)]
        print('v=\n', v)
        return v

    def AWGNPass(self,v,SNR): # SNR：雜訊強度
        N = len(v)    # 訊息長度
        s = 2*v-np.ones((1,N))      # 將1、0變成 1、 -1 (0不能做運算要換成-1)
        noise_var = 1/(10**(SNR/10))
        r = s + np.random.randn(1,N)*np.sqrt(noise_var)
        r = r[0].astype(np.int8)
        
        for i in range(N):
            if r[i] > 0:
                r[i] = 1
            else:
                r[i] = 0 
        return r
    
    def StructMap(self,g0,g1):
        allSignal  = np.array([(0,0),(1,0),(0,1),(1,1)]).astype(np.int8)
        temp = np.zeros([1,3]).astype(np.int8) # [input temp0 temp1]
        outputs = np.zeros([8,2]).astype(int) # [V0 V1] x8  
        nextState = np.zeros([8,2])  # [next0 next1] x8
        for  i in range(8):
            temp[0][1:] = allSignal[int(i/2)] # 放入暫存器
            if np.mod(i,2) == 0:  
                temp[0][0] = 0  # 偶數列是input = 0
                for j in range(3):
                    # 根據g0計算V0 output
                    if g0[j] == 1: 
                        outputs[i][0] += temp[0][j] 
                    # 根據g1計算V1 output
                    if g1[j] == 1:
                        outputs[i][1] += temp[0][j]
                outputs[i][0] = np.mod(outputs[i][0],2) # V1轉為0、1
                outputs[i][1] = np.mod(outputs[i][1],2) # V2轉為0、1
                nextState[i][1] = temp[0][1] # next1 是 NowState0(temp0)
                nextState[i][0] = 0 # next0是 input
                
            else:
                temp[0][0] = 1 # 奇數列是input = 1
                for j in range(3):
                    if g0[j] == 1: 
                        outputs[i][0] +=  temp[0][j]
                    if g1[j] == 1:
                        outputs[i][1] +=  temp[0][j]
                outputs[i][0] = np.mod(outputs[i][0],2)
                outputs[i][1] = np.mod(outputs[i][1],2)
                nextState[i][1] = temp[0][1]
                nextState[i][0] = 1
        print(outputs)
        print(nextState)
        outputsInAllSignal = [(allSignal == tuple(outputs[i])).all(axis=1).nonzero()[0][0] for i in range(8)] # outputs 是哪一個singal
        nextStateInAllSignal =  [(allSignal == tuple(nextState[i])).all(axis=1).nonzero()[0][0] for i in range(8)]
        print(nextStateInAllSignal)

        Base = [0,1,2,3]
        allEdges = []
        # for i in range(8):
        #     allEdges.append((int(i/2),nextStateInAllSignal[i]))
        allEdges = [ (int(i/2),nextStateInAllSignal[i]) for i in range(8)]

        G = nx.DiGraph()
        for i in range(8):
            G.add_edges_from([allEdges[i]], weight=outputsInAllSignal[i])  
        
        pos = nx.circular_layout(G)
        
        options = {
        'node_color': 'blue',
        'node_size': 800,
        'width': 1,
        'arrowstyle': '-|>',
        'arrowsize': 8,
        }
        print(dict([((u,v,),d['weight'])
                 for u,v,d in G.edges(data=True)]))
                
        edge_labels = dict([((u,v,),d['weight'])
                 for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)
        nx.draw_networkx(G,pos, connectionstyle='arc3, rad = 0.1', arrows=True, **options)
        pylab.show()
       
        
        

        # G = nx.DiGraph()
        # G.add_nodes_from([(0,0),(1,0),(0,1),(1,1)])
        
        

        return [outputs,nextState]

    def Decode(self,r):
        Dm = []
        return Dm 


if __name__ == '__main__':
   # ECC = ECC()
   g0 = np.array([1,0,0])
   g1 = np.array([0,1,1])
   ECC.StructMap(None,g0,g1)


    
