import numpy as np

def AWGNPass( v, SNR):  # SNR：雜訊強度
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
        return r.astype(np.int8)

if __name__ == "__main__":
    print(AWGNPass(np.array([0,1,1,0,1]),1))