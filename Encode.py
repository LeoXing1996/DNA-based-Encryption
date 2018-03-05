from PIL import Image
from  pylab import *
#from pandas import DataFrame,Series


class DNA_Encode:

    def __init__(self,raw):
        self.RAW = raw
        self.M = raw.shape[0]
        self.N = raw.shape[1]
        self.R_raw = self.RAW[:, :, 0].tolist()
        self.G_raw = self.RAW[:, :, 1].tolist()
        self.B_raw = self.RAW[:, :, 2].tolist()
        self.B1 = self.Generate_cha(self.M * self.N)
        self.B21 = self.Generate_cha(self.M)
        self.B22 = self.Generate_cha(self.N)
        self.B31 = self.Generate_cha(self.M)
        self.B32 = self.Generate_cha(self.N)
        self.X = self.Generate_Seq(self.B1)
        self.Y1 = self.Generate_Seq(self.B21)
        self.Y2 = self.Generate_Seq(self.B22)
        self.Z1 = self.Generate_Seq(self.B31)
        self.Z2 = self.Generate_Seq(self.B32)
        self.DNA = self.Get_DNA()

    # 生成索引序列
    def Generate_Seq(self, B):
        return [B.index(i) for i in sorted(B)]

    # 生成混沌序列，暂时用随机数组代替
    def Generate_cha(self,length):
        return list(rand(length))

    # 暂时使用现成的DNA序列
    def Get_DNA(self):
        self.DNA = 'cacttctcctacaacttaaatgtcttctccactgttacagcttctctaaatggcttctcc'\
                    'aatgttacagcttgtaccgatgacactgttacaggttgtacctacaaggaccagacggat' \
                    'ctcatagttggtccatacggagtaaaatgtgaaagcaattaggttggaatgtgaaggaag' \
                    'aacactgaggctaccatgaaactctcacaagacgcgactgaatagcaagcgaagctaatg' \
                    'ttacattataccggtgtgaaaactctagcccaacgtgaattattgataaccttgatatac' \
                    'caattcacacataacagcctcatgtgccacaaaacaatgcctaccactgcggttgttttc'

    # 图像整体置乱
    def plot_replacement(self):
        R_step1 = [[0 for i in range(self.N)] for j in range(self.M)]
        G_step1 = [[0 for i in range(self.N)] for j in range(self.M)]
        B_step1 = [[0 for i in range(self.N)] for j in range(self.M)]
        T = reshape(self.X, (self.M, self.N)).astype(int)
        for i in range(self.M):
            for j in range(self.N):
                u = (T[i][j]-1) // self.M
                v = (T[i][j]-1) % self.M
                R_step1[u][v] = self.R_raw[i][j]
                G_step1[u][v] = self.G_raw[i][j]
                B_step1[u][v] = self.B_raw[i][j]

        return R_step1, G_step1, B_step1

    # 输入二维list 进行比特置换
    def bits_replacement(self, c, bits):
        for i in range(self.M):
            for j in range(self.N):
                bits[i][j] = self.Butterfly_Network(bits[i][j], c)
        return bits

    # 子图划分
    def div2subplot(self, plt):
        # 按照从 左->右 87654321 的顺序排列
        # p1 是 1357 p2 是 2468
        p1 = [[0 for i in range(self.N)] for j in range(self.M)]
        p2 = [[0 for i in range(self.N)] for j in range(self.M)]
        for i in range(self.M):
            for j in range(self.N):
                p1[i][j] = [plt[i][j][2 * k + 1] for k in range(4)]
                p2[i][j] = [plt[i][j][2 * k] for k in range(4)]
        return p1, p2

    # 子图上移置乱
    def upward(self, plt, seq):
        plt_up = [[0 for i in range(self.N)] for j in range(self.M)]
        up = [i%self.M for i in seq]
        for col in range(self.N):
            if up[col] == 0:
                new_col = [plt[row][col] for row in range(self.M)]
            else:
                new_col = [plt[row][col] for row in range(up[col],self.M)] + \
                          [plt[row][col] for row in range(up[col])]
            for row in range(self.M):
                plt_up[row][col] = new_col[row]
        return plt

    # 子图左移置乱
    def leftward(self, plt, seq):
        plt_left = [[0 for i in range(self.N)] for j in range(self.M)]
        left = [i % self.N for i in seq]
        for i in range(self.M):
            if left[i] == 0:
                plt_left[i] = plt[i]
            else:
                plt_left[i] = plt[i][left[i]:] + plt[i][:left[i]]
        return plt

    # 子图下移
    def downward(self, plt, seq):
        plt_down = [[0 for i in range(self.N)] for j in range(self.M)]
        down = [i % self.M for i in seq]
        for col in range(self.N):
            if down[col] == 0:
                new_col = [plt[row][col] for row in range(self.M)]
            else:
                new_col = [plt[row][col] for row in range(self.M - down[col],self.M)] + \
                          [plt[row][col] for row in range(4-down[col])]
            for row in range(self.M):
                plt_down[row][col] = new_col[row]
        return plt_down

    # 子图右移
    def rightward(self, plt, seq):
        plt_right = [[0 for i in range(self.N)] for j in range(self.M)]
        right = [i % self.N for i in seq]
        for i in range(self.M):
            if right[i] == 0:
                plt_right[i] = plt[i]
            else:
                plt_right[i] = plt[i][self.N - right[i] : ] + \
                               plt[i][ : self.N - right[i]]
        return plt_right

    # 子图置乱
    def Subplot_replacement(self, plt, seq1, seq2, flag):
        if flag == 1:
            plt_left = self.leftward(plt,seq1)
            plt_up = self.upward(plt_left, seq2)
            return plt_up
        else:
            plt_right = self.rightward(plt, seq1)
            plt_down = self.downward(plt_right, seq2)
            return plt_down

    # 子图合并
    def subplot_merge(self, sub1, sub2):
        merge = [[0 for i in range(self.N)] for j in range(self.M)]
        for i in range(self.M):
            for j in range(self.N):
                merge[i][j] = []
                for k in range(4):
                    merge[i][j] += [sub2[i][j][k], sub1[i][j][k]]
        return merge

    # 像素替代
    def pixel_replacement(self, DNA, plt):
        # 分离每个像素加密需要的 DNA 序列
        DNA_seq = [[[] for i in range(self.N)] for j in range(self.M)]
        plt_xor = [[[] for i in range(self.N)] for j in range(self.M)]

        # 处理 DNA 序列，使之成与 plt 相同形态的 list
        # DNA 数据量不够 暂时使用单一序列代替
        for i in range(self.M):
            # DNA_seq[i] = [DNA[4*(j+4*i):4*(j+4*i)+4] for j in range(4)]
            DNA_seq[i] = [DNA for j in range(4)]
        for i in range(self.M):
            for j in range(self.N):
                xor = self.DNA2Bin(DNA_seq[i][j])
                plt_xor[i][j] = [plt[i][j][k]^xor[k] for k in range(8)]
        return plt_xor

    # 密文扩散
    def cipher_diffusion(self, plt):
        s = [0] + [int(i) for i in bin(127)[2:]]
        new_plt = [[[] for i in range(self.N)] for j in range(self.M)]
        for i in range(self.M):
            for j in range(self.N):
                new_plt[i][j] = [plt[i][j][k]^s[k] for k in range(8)]
                s = new_plt[i][j]
        return new_plt

    # 像素合并
    def pixel_merge(self,R, G, B):
        raw = [[[R[i][j],G[i][j],B[i][j]] for j in range(self.N)] for i in range(self.M)]
        return raw

    # 将十进制数转化为二进制数组
    def Int2Bin(self, Int):
        Bin = [int(i) for i in bin(Int)[2:] ]
        if len(Bin) < 8:
            Bin = [0 for i in range(8-len(Bin))] + Bin
        return Bin

    # 将DNA转化为二进制
    def DNA2Bin(self, DNA):
        c = []
        for i in DNA:
            if i == 'a':
                c += [0, 0]
            elif i == 'c':
                c += [0, 1]
            elif i == 'g':
                c += [1, 0]
            elif i == 't':
                c += [1, 1]
        return c

    # 将十进制数组转化为二进制数组
    def Raw2Bin(self, R_raw, G_raw, B_raw):
        for i in range(self.M):
            for j in range(self.N):
                R_raw[i][j] = self.Int2Bin(R_raw[i][j])
                G_raw[i][j] = self.Int2Bin(G_raw[i][j])
                B_raw[i][j] = self.Int2Bin(B_raw[i][j])
        return R_raw, G_raw, B_raw

    # 将二进制数组转化回十进制
    def bin2raw(self, Bin):
        raw = [[0 for i in range(self.N)] for j in range(self.M)]
        for i in range(self.M):
            for j in range(self.N):
                Int = ''
                for k in range(8):
                    Int += str(Bin[i][j][k])
                raw[i][j] = int(Int, 2)
        return raw

    # 比特置换网络
    def Butterfly_Network(self,b,c):
        s = self.Butterfly_layer1(b, c)
        t = self.Butterfly_layer2(s, c)
        d = self.Butterfly_layer3(t, c)
        return d

    # Butterfly置换网络——第一层
    def Butterfly_layer1(self, b, c):
        s = [-1 for i in range(8)]
        for i in range(4):
            if c[-1 - i] == 0:
                s[-1 - i] = b[-1 - i]
                s[-5 - i] = b[-5 - i]
            else:
                s[-5 - i] = b[-1 - i]
                s[-1 - i] = b[-5 - i]
        return s

    # Butterfly置换网络——第二层
    def Butterfly_layer2(self, s, c):
        t = [-1 for i in range(8)]
        for i in range(2):
            if c[-5 - i] == 0:
                t[-1 - i] = s[-1 - i]
                t[-3 - i] = s[-3 - i]
            else:
                t[-3 - i] = s[-1 - i]
                t[-1 - i] = s[-3 - i]

            if c[-7 - i] == 0:
                t[-5 - i] = s[-5 - i]
                t[-7 - i] = s[-7 - i]
            else:
                t[-7 - i] = s[-5 - i]
                t[-5 - i] = s[-7 - i]
        return t

    # Butterfly置换网络——第三层
    def Butterfly_layer3(self, t, c):
        d = [-1 for i in range(8)]
        for i in range(4):
            if c[-1 - i] ^ c[-5 - i] == 0:
                d[-2 * i - 1] = t[-2 * i - 1]
                d[-2 * i - 2] = t[-2 * i - 2]
            else:
                d[-2 * i - 1] = t[-2 * i - 2]
                d[-2 * i - 2] = t[-2 * i - 1]
        return d

    # 加密函数，用于外部调用
    def Encoding(self):
        c = self.DNA2Bin(self.DNA[:8])
        R_step1, G_step1, B_step1 = self.plot_replacement() #像素置乱
        R_bin ,G_bin, B_bin = self.Raw2Bin(R_step1, G_step1, B_step1) # 转化为二进制

        # 比特置换
        R_step2 = self.bits_replacement(c, R_bin)
        G_step2 = self.bits_replacement(c, G_bin)
        B_step2 = self.bits_replacement(c, B_bin)

        # 子图划分
        R_p1, R_p2 = self.div2subplot(R_step2)
        G_p1, G_p2 = self.div2subplot(G_step2)
        B_p1, B_p2 = self.div2subplot(B_step2)

        # 子图置乱
        R_p1_re = self.Subplot_replacement(R_p1, self.Y1, self.Y2, 1)
        G_p1_re = self.Subplot_replacement(G_p1, self.Y1, self.Y2, 1)
        B_p1_re = self.Subplot_replacement(B_p1, self.Y1, self.Y2, 1)

        R_p2_re = self.Subplot_replacement(R_p1, self.Z1, self.Z2, 2)
        G_p2_re = self.Subplot_replacement(G_p1, self.Z1, self.Z2, 2)
        B_p2_re = self.Subplot_replacement(B_p1, self.Z1, self.Z2, 2)

        # 子图合并
        R_step3 = self.subplot_merge(R_p1_re, R_p2_re)
        G_step3 = self.subplot_merge(G_p1_re, G_p2_re)
        B_step3 = self.subplot_merge(B_p1_re, B_p2_re)

        # 像素替代
        DNA_pix = self.DNA[:4]    #临时使用前4位作为整个图的加密序列
        R_step4 = self.pixel_replacement(DNA_pix, R_step3)
        G_step4 = self.pixel_replacement(DNA_pix, G_step3)
        B_step4 = self.pixel_replacement(DNA_pix, B_step3)

        # 密文扩散
        R_step5 = self.cipher_diffusion(R_step4)
        G_step5 = self.cipher_diffusion(G_step4)
        B_step5 = self.cipher_diffusion(B_step4)

        # 二进制转化回十进制
        R_int = self.bin2raw(R_step5)
        G_int = self.bin2raw(G_step5)
        B_int = self.bin2raw(B_step5)

        # 像素合并成一个图
        encoded = self.pixel_merge(R_int, G_int, B_int)
        return encoded



raw = Image.open('lena_std.jpg')
raw = array(raw)
Encoder = DNA_Encode(raw)
Encoder.Encoding()