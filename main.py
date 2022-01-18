from cProfile import label
import numpy as np
from numpy.lib.function_base import append
from numpy.random.mtrand import rand
import networkx as nx
import matplotlib.pyplot as plt
import scipy as sp
from decimal import Decimal


class ECC:

    def __init__(self,g0,g1,N):
        self.N = N
        self.g0 = np.array(g0)
        self.g1 = np.array(g1)
        [self.G,outputsInAllSignal, nextStateInAllSignal] = self.StructMap(
            self.g0, self.g1)
        self.base = np.array([(0, 0), (1, 0), (0, 1), (1, 1)]).astype(np.int8)
        self.allPath = []
        self.allDistance = []
        self.poptimes = 0
        self.sTimes = int(self.N + len(self.g0) -1)
        self.DFS(self.G,0,[],0)
        self.allEdges = [(int(i/2), nextStateInAllSignal[i]) for i in range(8)]
        pass

    def Encode(self, u):
        g0 = self.g0
        g1 = self.g1
        v = []
        G = np.zeros([len(u), (len(u)+len(g0)-1) * 2])
        g = []

        for i in range(len(g0)):
            g.append(g0[i])
            g.append(g1[i])


        for i in range(len(G)):
            G[i][2*i:2*i+2*len(g1)] = g

        v = np.matmul(u, G).astype(np.int8)
        v = np.mod(v, 2)
        return v

    def AWGNPass(self, v, SNR):  # SNR：雜訊強度
        N = len(v)    # 訊息長度
        s = 2*v-np.ones((1, N))      # 將1、0變成 1、 -1 (0不能做運算要換成-1)
        noise_var = 1/(10**(SNR/10))
        r = s + np.random.randn(1, N) * np.sqrt(noise_var)
        r = r[0]

        for i in range(N):
            if r[i] > 0:
                r[i] = 1
            else:
                r[i] = 0
        return r

    def StructMap(self, g0, g1):
        allSignal = np.array([(0, 0), (1, 0), (0, 1), (1, 1)]).astype(np.int8)
        temp = np.zeros([1, 3]).astype(np.int8)  # [input temp0 temp1]
        outputs = np.zeros([8, 2]).astype(int)  # [V0 V1] x8
        nextState = np.zeros([8, 2])  # [next0 next1] x8
        for i in range(8):
            temp[0][1:] = allSignal[int(i/2)]  # 放入暫存器
            if np.mod(i, 2) == 0:
                temp[0][0] = 0  # 偶數列是input = 0
                for j in range(3):
                    # 根據g0計算V0 output
                    if g0[j] == 1:
                        outputs[i][0] += temp[0][j]
                    # 根據g1計算V1 output
                    if g1[j] == 1:
                        outputs[i][1] += temp[0][j]
                outputs[i][0] = np.mod(outputs[i][0], 2)  # V1轉為0、1
                outputs[i][1] = np.mod(outputs[i][1], 2)  # V2轉為0、1
                nextState[i][1] = temp[0][1]  # next1 是 NowState0(temp0)
                nextState[i][0] = 0  # next0是 input

            else:
                temp[0][0] = 1  # 奇數列是input = 1
                for j in range(3):
                    if g0[j] == 1:
                        outputs[i][0] += temp[0][j]
                    if g1[j] == 1:
                        outputs[i][1] += temp[0][j]
                outputs[i][0] = np.mod(outputs[i][0], 2)
                outputs[i][1] = np.mod(outputs[i][1], 2)
                nextState[i][1] = temp[0][1]
                nextState[i][0] = 1

        outputsInAllSignal = [(allSignal == tuple(outputs[i])).all(axis=1).nonzero()[
            0][0] for i in range(8)]  # outputs 是哪一個singal
        nextStateInAllSignal = [(allSignal == tuple(nextState[i])).all(
            axis=1).nonzero()[0][0] for i in range(8)]

        allEdges = [(int(i/2), nextStateInAllSignal[i]) for i in range(8)]

        G = nx.DiGraph()
        for i in range(8):
            G.add_edges_from([allEdges[i]], weight=outputsInAllSignal[i])

        pos = nx.circular_layout(G)

        options = {
            'node_color': 'tab:blue',
            'node_size': 800,
            'width': 1,
            'arrowstyle': '-|>',
            'arrowsize': 8,
        }
        print(dict([((u,v,),d['weight'])
                 for u,v,d in G.edges(data=True)]))

        edge_labels = dict([((u, v,), d['weight'])
                            for u, v, d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        nx.draw_networkx(
            G, pos, connectionstyle='arc3, rad = 0.1', arrows=True, **options)

        
        return [G,outputsInAllSignal, nextStateInAllSignal]  # input從index判斷所以省略

    
    def BitDistance(self,bit1,bit2): # (0,1) 
        return int(bit1[0]^bit2[0]) + int(bit1[1]^bit2[1]) # ^是XOR運算子
    
  
    def WhichBase(self,l):
        return [(self.base == tuple(l[i])).all(
            axis=1).nonzero()[0][0] for i in range(len(l))]
    
    def DFS(self,G,vertex,queue,stack):
        queue.append(vertex)
        stack += 1
        for i in list(G.neighbors(vertex)):
            if stack < self.sTimes+1:
                self.DFS(G,i,queue.copy(),stack)
            if stack == self.sTimes+1:
                self.allPath.append(queue.copy())
                queue.pop()
                self.poptimes += 1
                if self.poptimes >= 2:
                    self.poptimes = 0
                    queue.pop()
                    return 
                return

        return 
    
    def Decode(self, v ):
        G = self.G
        v = [(v[2*i], v[2*i+1]) for i in range(int(len(v)/2))]
        
        v2base = self.WhichBase(v)

        for i in self.allPath:
            theDistance = 0
            for j in range(len(i)-1):
                pathSelect = (i[j],i[j+1])
                pathSelect = G.get_edge_data(*pathSelect)['weight']
                bit1 = self.base[pathSelect]
                bit2 = self.base[v2base[j]]
                theDistance += self.BitDistance(bit1,bit2)    
              
            self.allDistance.append(theDistance)
            theDistance = 0  

        minIndex = self.allDistance.index(min(self.allDistance))
        path = self.allPath[minIndex]
        vDecode = []
        self.allDistance.clear()
        for i in range(len(path)-1):
            s = (path[i],path[i+1])
            vDecode.append(np.mod(self.allEdges.index(s),2))

        return vDecode




if __name__ == '__main__':
    N = 5
    ecc = ECC([1,1,1],[1,0,1],N)
    errTimes = 0
    errTimes2 = 0
    CodeNum =5000
    dB = [i for i in range(1,11)]
    ber = []
    ber2 = []
    for i in dB:
        for j in range(CodeNum):
            u = np.round(np.random.rand(1,N)).astype(np.int8)[0]
            v = ecc.Encode(u)
            r = ecc.AWGNPass(v,i)
            r2 = ecc.AWGNPass(u,i)
            du = ecc.Decode(r)
            du = du[:N]
            errTimes += np.sum(np.mod(du-u,2))
            errTimes2 += np.sum(np.mod(r2-u,2))
        ber.append(Decimal(errTimes/CodeNum/N))
        ber2.append(Decimal(errTimes2/CodeNum/N))
        errTimes = 0
        errTimes2 = 0

    plt.figure(2)   
    plt.title("Conv Code BER Plot ")
    plt.plot(dB,ber,label="Coded")
    plt.plot(dB,ber2,label="Uncoded")
    plt.legend(loc='upper right')
    plt.yscale('log')
    plt.ylim(bottom=10**(-4))
    plt.ylabel("BER")
    plt.xlabel("Eb/No (dB)")
    plt.show()
    


