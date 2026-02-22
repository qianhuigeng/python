import math

# 全局变量定义（对应Fortran的common块）
mozu = 0.0
MSDT = 0  # 全局控制变量
FAX = [0.0 for _ in range(27)]  # 索引1-26使用（对应Fortran的1-based）
FAY = [0.0 for _ in range(27)]
x = [0.0 for _ in range(1362)]  # 索引1-1361使用
y = [0.0 for _ in range(1362)]


def moditu(QM, SS, DD):
    """对应原Fortran的MODITU子程序：计算mozu2（从x,y数组插值）"""
    global x, y
    RE1 = abs(QM) / SS * DD / (1e-6)  # 计算雷诺数相关值
    i = 0
    while True:
        i += 1
        j = i + 1
        if i >= 1361:
            j = 1361
            return y[1361]
        x1 = x[i]
        x2 = x[j]
        y1 = y[i]
        y2 = y[j]
        if RE1 >= x1 and RE1 < x2:
            # 线性插值计算
            return y1 + (RE1 - x1) * (y2 - y1) / (x2 - x1)
        if i >= 1350:
            break
    return y[1361]  # 超出范围时返回最后一个值


def nodeshangku(HU, S1, B2, H2, Q2, H1U, Q1U, H1D, Q1D, DR, DT):
    """对应原Fortran的NODESHANGKU子程序：上库节点计算"""
    global mozu, MSDT
    if MSDT == 1:
        Q1U_val = Q1D
        H1U_val = H1D
        # 调用MODITU更新mozu
        mozu = moditu(Q1D, S1, DR)

        R1 = mozu * 1000.0 * DT / (S1 ** 2) / DR / 2.0 / 9.81
        BM = B2 + R1 * abs(Q2)
        CM = H2 - B2 * Q2
        BMX = BM + 1.0 / 19.62 / (S1 ** 2) * Q1U_val
        Q1D_new = (HU - CM) / BMX
        H1D_new = CM + BM * Q1D_new
    else:
        Q1U_val = Q1D
        H1U_val = H1D
        # 调用MODITU更新mozu
        mozu = moditu(Q1D, S1, DR)

        R1 = mozu * 1000.0 * DT / (S1 ** 2) / DR / 2.0 / 9.81
        BM = B2 + R1 * abs(Q2)
        CM = H2 - B2 * Q2
        H1D_new = HU
        Q1D_new = (HU - CM) / BM
        # if Q1D_new < 0.0:
        #     Q1D_new = 0.0
        #     H1D_new = CM
    return H1U_val, Q1U_val, H1D_new, Q1D_new


def taotk(XLVK, TAO0, DT):
    """对应原Fortran的TAOTK子程序：阀门开启计算"""
    TAO = TAO0 + DT / XLVK
    if TAO >= 1.0:
        TAO = 1.0
    return TAO


def taotg(XLVG1, TAO0, DT):
    """对应原Fortran的TAOTG子程序：阀门关闭计算"""
    ZDIAN1 = 0.0
    if TAO0 >= ZDIAN1:
        TAO = TAO0 - DT / XLVG1
    else:
        TAO = TAO0 - DT / XLVG1
    if TAO <= 0.001:
        TAO = 0.001
    return TAO


def nodetyfa(S1, DM, T, DT, TAO, DDR, KG, XLVK, XLVG, B2, R2, H2, Q2, H2U, Q2U, H2D, Q2D, HD):
    """对应原Fortran的NODETYFA子程序：阀门节点计算"""
    global mozu, FAX, FAY, x, y
    H2U_val = H2D
    Q2U_val = Q2D

    # 调用MODITU更新mozu
    mozu = moditu(Q2D, S1, DM)

    R2_val = mozu * 1000.0 * DT / (S1 ** 2) / DM / 2.0 / 9.81
    CP2 = H2 + B2 * Q2
    BP2 = B2 + R2_val * abs(Q2)
    XX = TAO

    # 根据阀门动作参数更新TAO
    if KG == 2 and T >= 0.0:
        TAO = taotk(XLVK, XX, DT)
    if KG == 3 and T >= 0.0:
        TAO = taotg(XLVG, XX, DT)

    if TAO <= 0.001:
        return H2U_val, Q2U_val, HD, 0.0, TAO  # Q2D=0.0

    # 读取DIEFA.DAT文件获取FAX和FAY
    with open('DIEFA.DAT', 'r') as f:
        for i in range(1, 27):
            line = f.readline().strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    FAX[i] = float(parts[0])
                    FAY[i] = float(parts[1])

    ALFA = 0.0
    # 查找插值区间
    for i in range(1, 27):
        if TAO <= FAX[i]:
            if TAO <= FAX[1]:
                TAO = FAX[1]
                ALFA = FAY[1]
                break
            XYZ1 = FAX[i]
            XYZ2 = FAX[i - 1]
            FXYZ1 = FAY[i]
            FXYZ2 = FAY[i - 1]
            # 线性插值
            ALFA = (FXYZ1 - FXYZ2) / (XYZ1 - XYZ2) * (TAO - XYZ2) + FXYZ2
            break

    FA = (DDR ** 2) * 3.1415926 / 4.0
    XSU1 = ALFA * abs(Q2D) / (2.0 * 9.81) / (FA ** 2)
    QP = (CP2 - HD) / (BP2 + XSU1)
    Q2D_new = QP
    H2D_new = CP2 - BP2 * QP

    return H2U_val, Q2U_val, H2D_new, Q2D_new, TAO


def nbdm1(QG, HG, BG, RG):
    """对应原Fortran的NBDM1子程序：内部节点计算"""
    CP = HG[0] + BG * QG[0]
    BP = BG + RG * abs(QG[0])
    CM = HG[1] - BG * QG[1]
    BM = BG + RG * abs(QG[1])
    QJ = (CP - CM) / (BP + BM)
    HJ = CP - BP * QJ
    return QJ, HJ


def nbdm(MNS, M, S1, D, Q, H, B, R, NS, MSS, DT):
    """对应原Fortran的NBDM子程序：内部节点批量计算"""
    global mozu, x, y
    for i in range(1, M + 1):
        QG = [0.0, 0.0]
        HG = [0.0, 0.0]
        NB = int(NS[i - 1]) - 2
              # Python列表0-based

        if NB == 0:
            continue
        K1 = int(MSS[i - 1])  # 支路第一断面位置

        # 复制上一时刻数据
        for j in range(1, NB + 1):
            Q[0][int(K1 + j - 1)] = Q[1][int(K1 + j - 1)]  # 调整索引适配Python
            H[0][int(K1 + j - 1)] = H[1][int(K1 + j - 1)]

        # 计算内部节点
        for j in range(1, NB + 1):
            QG[0] = Q[1][K1 + j - 2] if (K1 + j - 2) >= 0 else 0.0
            HG[0] = H[1][K1 + j - 2] if (K1 + j - 2) >= 0 else 0.0
            QG[1] = Q[1][K1 + j] if (K1 + j) < len(Q[1]) else 0.0
            HG[1] = H[1][K1 + j] if (K1 + j) < len(H[1]) else 0.0

            # 调用MODITU更新mozu
            current_Q = Q[1][K1 + j - 1]
            mozu = moditu(current_Q, S1, D)
            R[i - 1] = mozu * 1000.0 * DT / (S1 ** 2) / D / 2.0 / 9.81

            BG = B[i - 1]
            RG = R[i - 1]
            QJ, HJ = nbdm1(QG, HG, BG, RG)

            Q[1][K1 + j - 1] = QJ
            H[1][K1 + j - 1] = HJ


def main():
    """主程序：对应原Fortran主程序逻辑"""
    global mozu, MSDT, FAX, FAY, x, y

    # 初始化变量
    NXX0 = 1
    # 分配数组（1-based索引适配，实际用0-based列表）
    NS = [0 for _ in range(NXX0)]
    MSS = [0 for _ in range(NXX0)]
    NFD = [0.0 for _ in range(NXX0)]
    S = [0.0 for _ in range(NXX0)]
    XSU1 = [0.0 for _ in range(NXX0)]
    BSU = [0.0 for _ in range(NXX0)]
    B = [0.0 for _ in range(NXX0)]
    R = [0.0 for _ in range(NXX0)]
    MJD = [0.0 for _ in range(NXX0)]

    # 读取管道数据文件moditu4.DATA（1-based索引）
    with open('moditu4.DATA', 'r') as f:
        for i in range(1, 1362):
            line = f.readline()
            if not line:
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                x[i] = float(parts[0])
                y[i] = float(parts[1])

    # 参数初始化
    MSDT = 1
    DT = 10 / 1000.0 / 16.0
    PL = 10
    BSU[0] = 1000.0  # BSU(1)


    NFD[0] = round(PL / BSU[0] / DT)

    DR = 0.1
    S[0] = 3.1415926 * DR * DR / 4.0  # S(1)

    DDR = 0.1  # 阀参数
    KG = 3.0  # 动作参数
    XLVK = 30.0
    XLVG = 0.009

    # 初始化NS数组
    for i in range(NXX0):
        NS[i] = NFD[i] + 1

    MSS[0] = 1  # MSS(1)

    # 初始化MJD数组

    MJD = [int(MSS[i] + NFD[i]) for i in range(NXX0)]

    MNS = MJD[NXX0 - 1]

    # 分配H和Q二维列表（2行，MNS列，1-based索引）
    H = [
        [0.0 for _ in range(int(MNS) + 1)],  # H(1,I)
        [0.0 for _ in range(int(MNS) + 1)]  # H(2,I)
    ]
    Q = [
        [0.0 for _ in range(int(MNS) + 1)],  # Q(1,I)
        [0.0 for _ in range(int(MNS) + 1)]  # Q(2,I)
    ]

    # 计算Q0并调用MODITU
    Q0 = 0.3 * 3.1415926 * (0.022 ** 2) / 4.0
    mozu = moditu(Q0, S[0], DR)

    # 计算XSU1(1)
    XSU1[0] = mozu * PL / (S[0] ** 2) / DR / 2 / 9.81

    # 其他参数设置
    HU = 100
    HD = 0.0
    V1 = Q0 / S[0]

    # 再次调用MODITU
    mozu = moditu(Q0, S[0], DR)
    print(f"mozu: {mozu}")
    input("Press Enter to continue...")

    # 初始化B和R数组
    for i in range(NXX0):
        B[i] = BSU[i] / 9.81 / S[i]
        R[i] = XSU1[i] / NFD[i]

    # 初始化Q数组
    for i in range(1, int(MJD[0]) + 1):
        Q[0][i] = Q0  # Q(1,I)

    # 第一段总水头计算
    for i in range(1, int(MJD[0]) + 1):
        H[0][i] = HU - XSU1[0] / NFD[0] * (i - 1) * Q0 * Q0

    # 不考虑流速水头时的测压管水头计算
    if MSDT == 1:
        for i in range(int(MSS[0]), int(MJD[0]) + 1):
            H[0][i] = H[0][i] - (V1 ** 2) / (2.0 * 9.81)

    # 复制数据到H(2,I)和Q(2,I)
    for i in range(1, int(MNS) + 1):
        Q[1][i] = Q[0][i]
        H[1][i] = H[0][i]

    # 计算过阀损失
    DH = H[1][int(MJD[0])]
    DDH = DH - HD
    ALFAM = DDH * 2.0 * 9.81 * (S[0] ** 2) / (Q0 ** 2)

    # 读取DIEFA.DAT文件
    with open('DIEFA.DAT', 'r') as f:
        for i in range(1, 27):
            line = f.readline().strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    FAX[i] = float(parts[0])
                    FAY[i] = float(parts[1])

    # 查找ALFAM对应的TAO0
    TAO0 = 0.0
    for i in range(26, 0, -1):
        if ALFAM <= FAY[i]:
            if ALFAM <= FAY[26]:
                ALFA = FAY[26]
                TAO0 = FAX[26]
                break
            XYZ1 = FAX[i]
            XYZ2 = FAX[i + 1] if (i + 1) <= 26 else FAX[i]
            FXYZ1 = FAY[i]
            FXYZ2 = FAY[i + 1] if (i + 1) <= 26 else FAY[i]
            TAO0 = (XYZ1 - XYZ2) / (FXYZ1 - FXYZ2) * (ALFAM - FXYZ2) + XYZ2
            break

    TAO = TAO0
    print(f"DDH: {DDH}, TAO: {TAO}, Q0: {Q0}, DH: {DH}, V1: {V1}, ALFAM: {ALFAM}")
    input("Press Enter to continue...")

    # 初始化时间和输出文件
    T = 0.0
    with open('UNSTEAY.DAT', 'w') as f:
        # 写入初始数据
        line = f"{T:10.4f}  {TAO:12.5f}  {mozu:12.5f}  {H[1][1]:12.5f}  {Q[1][1]:12.5f}  {H[1][9]:12.5f}  {Q[1][9]:12.5f}  {H[1][int(MJD[0])]:12.5f}  {Q[1][int(MJD[0])]:12.5f}\n"
        f.write(line)

    NN = 0
    # 时间循环
    while True:
        T += DT
        # 上库节点计算
        H1U, Q1U, H1D, Q1D = nodeshangku(
            HU, S[0], B[0], H[1][2], Q[1][2],
            H[0][1], Q[0][1], H[1][1], Q[1][1], DR, DT
        )
        H[0][1], Q[0][1], H[1][1], Q[1][1] = H1U, Q1U, H1D, Q1D

        # 阀门节点计算
        NX = 1
        NXY = int(MJD[NX - 1])
        H2U, Q2U, H2D, Q2D, TAO = nodetyfa(
            S[NX - 1], DR, T, DT, TAO, DDR, KG, XLVK, XLVG,
            B[NX - 1], R[NX - 1], H[1][NXY - 1], Q[1][NXY - 1],
            H[0][NXY], Q[0][NXY], H[1][NXY], Q[1][NXY], HD
        )
        H[0][NXY], Q[0][NXY], H[1][NXY], Q[1][NXY] = H2U, Q2U, H2D, Q2D

        # 内部节点计算
        nbdm(int(MNS), NXX0, S[0], DR, Q, H, B, R, NS, MSS, DT)

        # 输出数据
        NN += 1
        YYY = NN
        NNN = int(YYY)
        if YYY == NNN:
            with open('UNSTEAY.DAT', 'a') as f:
                line = f"{T:10.4f}  {TAO:12.5f}  {mozu:12.5f}  {H[1][1]:12.5f}  {Q[1][1]:12.5f}  {H[1][9]:12.5f}  {Q[1][9]:12.5f}  {H[1][int(MJD[0])]:12.5f}  {Q[1][int(MJD[0])]:12.5f}\n"
                f.write(line)

        # 时间判断
        if T > 3.0:
            break

    print("计算完成，结果已写入UNSTEAY.DAT")


if __name__ == "__main__":
    main()