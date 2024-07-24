import itertools

from __init__ import *


def graph(data=None, label=['JW mapping', 'BK mapping', 'TT mapping'], N=[2 * i for i in range(2, 40)]):
    import matplotlib.pyplot as plt
    import matplotlib
    def n2(n):
        return n * (n + 1)

    def n4(n):
        return n * (n + 1) * n * (n + 1)
    if data is None:
        # for i in range(K):
        JW = [1.8461538461538463, 2.2777777777777777, 2.72, 3.1671641791044776, 3.6129032258064515, 4.054878048780488,
              4.492569002123142, 4.926153846153846, 5.3560414269275025, 5.782685512367491, 6.206514206514207,
              6.627906976744186, 7.047191011235955, 7.464644970414201, 7.880505081613797, 8.294971487817522,
              8.708214049768774, 9.120377358490567, 9.531584948688712, 9.941942792410082, 10.351542177629135,
              10.760462074978204, 11.168771084337349, 11.576529038889841, 11.983788330656878, 12.390595009596929,
              12.79698969835297, 13.203008356545961, 13.608682921904185, 14.014041850220265, 14.419110571894741,
              14.823911879446186, 15.228466257668712, 15.632792165957994, 16.036906280594536, 16.440823703379156,
              16.844558141891035, 17.2481220657277]
        BK = [1.7307692307692308, 2.685185185185185, 2.683333333333333, 3.8656716417910446, 4.185867895545315,
              4.308797909407666, 4.290074309978769, 5.3364102564102565, 5.834522439585731, 6.1044812078380986,
              6.349387849387849, 6.60418263906636, 6.750786516853933, 6.6767258382643, 6.6326801663073605, 7.38764370444912,
              7.916146712667302, 8.18364448857994, 8.497165662159961, 8.678945664926973, 8.88728421337117,
              8.951290701641332, 9.16463453815261, 9.331893095768374, 9.531952753219098, 9.586753191359717,
              9.747382729539908, 9.729547593891077, 9.752745385094318, 9.644443358912415, 9.594618091871313,
              10.039672329506338, 10.456435531267722, 10.70467078984667, 11.04200083728075, 11.17850878426485,
              11.35944604933219, 11.431674491392801]
        TT1 = [1.8461538461538463, 2.2685185185185186, 3.0633333333333335, 3.497014925373134, 4.056067588325653,
               4.461672473867596, 4.6289808917197455, 4.6290598290598295, 5.123590333716916, 5.608095085126887,
               5.8093093093093096, 5.867620751341682, 6.110497592295345, 6.499753451676528, 6.810093932861102,
               7.049492269691703, 7.20746776285204, 7.460734856007944, 7.631975891839062, 7.777528286874081,
               7.87637241985068, 7.971470944998294, 8.056622489959839, 8.109313003255096, 8.138908006423492,
               8.11334277793011, 8.33485419952942, 8.500268946306791, 8.70057886146028, 8.913426353559755,
               9.080576967326277, 9.250969497622906, 9.382843738722483, 9.459082040237819, 9.494489275224732,
               9.550730998263722, 9.570108339921001, 9.587527988443481]
        TT2 = [2.3076923076923075, 2.462962962962963, 3.0733333333333333, 3.4776119402985075, 4.276497695852535,
               4.69904181184669, 4.804936305732484, 4.774188034188034, 5.26156501726122, 5.660857693543206,
               5.819415569415569, 5.882826475849732, 6.062054574638844, 6.539077909270217, 6.913285340314136,
               7.1194614704357635, 7.326347598424234, 7.5809831181727905, 7.712607916598794, 7.902118649782201,
               7.9779286735808475, 8.055759827148327, 8.137911646586346, 8.202131231797155, 8.213428156304962,
               8.18432196935076, 8.374080210086625, 8.526258764768034, 8.754343985057297, 8.948889204680025,
               9.104605070493244, 9.281160105607439, 9.40782182811775, 9.465179376737165, 9.50062960875255,
               9.556560777876294, 9.569634427752723, 9.596029854339713]
        TT = [1.3846153846153846, 2.2962962962962963, 2.7066666666666666, 3.1970149253731344, 3.5929339477726576,
              3.9050522648083623, 4.271231422505308, 4.605470085470086, 4.87433831990794, 5.158849983938323,
              5.451489951489951, 5.671948206831928, 5.912295345104334, 6.155325443786983, 6.34139205420388,
              6.538651541487512, 6.7287675255083315, 6.883296921549156, 7.034696204593582, 7.184340062844735,
              7.305023478936523, 7.467552405140063, 7.6491084337349395, 7.793661127291417, 7.936612137575661,
              8.07641999004763, 8.189641953876116, 8.301534914993757, 8.41006949702997, 8.49982828856994, 8.622338200225913,
              8.76622879912102, 8.883377841934319, 8.99726103042686, 9.106782962984939, 9.197049067250065,
              9.284652577722262, 9.369589502828939]
        TT = [1.4615384615384615, 2.1296296296296298, 2.6166666666666667, 3.0880597014925373, 3.4846390168970816,
              3.8192508710801394, 4.164808917197452, 4.52051282051282, 4.789067894131185, 5.1025538066174105,
              5.408235158235159, 5.596217735752619, 5.833322632423756, 6.079585798816568, 6.271327379119187,
              6.466258042875004, 6.662360224130759, 6.818550148957299, 6.972405929304447, 7.124273442030451,
              7.248020900194813, 7.407793487737386, 7.591421686746988, 7.737221175261264, 7.880753642625631,
              8.022260813047762, 8.137093860429866, 8.249405436557486, 8.359722011880121, 8.45074842499171,
              8.575769777852152, 8.722669207179266, 8.842220962004435, 8.95814236015241, 9.070117318196353,
              9.162391593630982, 9.25404758318125, 9.342771156855664]
        data = [JW,BK,TT]
    # BK_theory = [(log(n)/log(2) * 4 * n4(n) + log(n)/log(2) * 2 * n2(n)) / (n2(n) + n4(n)) for n in N]
    # TT_theory =  [(log(n)/log(3) * 4 * n4(n) + log(n)/log(3) * 2 * n2(n)) / (n2(n) + n4(n)) for n in N]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    f = 25
    k = 0
    font = {'family': "serif",
            'weight': "normal",
            'size': 10}

    matplotlib.rc('font', **font)

    # ax.plot(R[k:f], GE4[k:f], , label="TT in UCCSD, 4 кубита")
    ax.set_xlabel('Number of modes (qubits)', fontsize=12)
    ax.set_ylabel('Average Pauli weight', fontsize=12)
    for i, dat in enumerate(data):
        ax.plot(N, dat,  label=label[i])

    # ax.plot(N, BK_theory, '-', label="BK theory")
    # ax.plot(N, TT_theory, '-', label="TT theory")
    # ax.plot(N, JW, 'bs', label="JW mapping")
    # ax.set_xscale("log", base = 2 )
    # ax.set_yscale("log", base = 2)

    ax.yaxis.get_label()
    ax.xaxis.get_label()
    ax.grid()
    ax.legend()
    # ax.set_title(label="Mоделирование молекулы H_2")
    plt.savefig("compare.png", dpi=200)


def st_enumeration1(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []
    j = 0
    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, nmodes - 1):
            if j <= 1:
                j += 1
                num_list = num_list + [i, i + 2]
            else:
                if j == 2:
                    j = 3
                else:
                    j = 0
        if nmodes // 2 % 2 != 0:
            num_list = num_list + [nmodes - 1, nmodes]
        j = 0
        for i in range(1, nmodes - 1):
            if j <= 1:
                j += 1
                num_list = num_list + [i + nmodes, i + 2 + nmodes]
            else:
                if j == 2:
                    j = 3
                else:
                    j = 0
        if nmodes // 2 % 2 != 0:
            num_list = num_list + [nmodes * 2 - 1, 2 * nmodes]
    return num_list


def st_enumeration2(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        j = 1
        for i in range(1, 18):
            num_list = num_list + [i * j + (nmodes - i + 1) * (1 - j)]
            j = (j + 1) % 2
        j = 1
        for i in range(1, 19):
            num_list = num_list + [(i + 1) * j + (nmodes - i + 2) * (1 - j)]
            j = (j + 1) % 2
        j = 1
        for i in range(37, 55):
            num_list = num_list + [i * j + (3 * nmodes - i + 1) * (1 - j)]
            j = (j + 1) % 2
        j = 1
        for i in range(37, 55):
            num_list = num_list + [(i + 1) * j + (3 * nmodes - i + 2) * (1 - j)]
            j = (j + 1) % 2

    return num_list


def st_enumeration3(nmodes: int = 30):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    # nmodes = 36
    k = nmodes // 4 + 1
    j = nmodes // 2 % 2
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, k):
            num_list = num_list + [2 * i - 1] + [nmodes - 2 * i + 1]
        if j == 1:
            num_list = num_list + [nmodes // 2]

        for i in range(1, k):
            num_list = num_list + [2 * i] + [nmodes - 2 * i + 2]
        if j == 1:
            num_list = num_list + [nmodes // 2 + 1]

        for i in range(1, k):
            num_list = num_list + [nmodes + 2 * i - 1] + [2 * nmodes - 2 * i + 1]
        if j == 1:
            num_list = num_list + [nmodes // 2 + nmodes]

        for i in range(1, k):
            num_list = num_list + [nmodes + 2 * i] + [2 * nmodes - 2 * i + 2]
        if j == 1:
            num_list = num_list + [nmodes // 2 + nmodes + 1]

    return num_list


def st_enumeration4(nmodes: int = 30):
    # nmodes = 36
    k = nmodes // 4 + 1
    j = nmodes // 2 % 2
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, nmodes + 1):
            num_list = num_list + [i,i + nmodes]
    return num_list


def st_enumeration5(nmodes: int = 30):
    # nmodes = 36
    k = nmodes // 4 + 1
    j = nmodes // 2 % 2
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, nmodes + 1):
            num_list = num_list + [2*i - 1, 2*i]
    return num_list


def main(tt: TernaryTree, only2prod=False, only4prod=False):
    nmodes = len(tt.branches()) // 2

    def valid_2prod(i, j):
        par = (i + j) % 2 == 1 or ((i - j) % 2 == 1)
        return par

    def valid_4prod(i, j, k, l):
        # par = (i + j + k + l) % 2 == 0
        par = (i == (k - nmodes)) and (j == (l - nmodes))
        return par


    # draw(tt, k = 1)

    s = 0
    n = 0
    k = 2
    if only2prod:
        # nmodes = nmodes//2
        k = 2
    if not only4prod:
        for i in range(1, k * nmodes + 1):
            if i <= nmodes:
                m1 = nmodes + 1
            else:
                m1 = 2 * nmodes + 1
            for j in range(i, m1):
                if valid_2prod(i, j):
                    r = pauli_weight(tt, [i, j])
                    if r > 0:
                        s += r
                        n += 1
                    else:
                        print('пиздец')

    if not only2prod:
        for i in range(1, 2 * nmodes + 1):
            if i <= nmodes:
                m1 = nmodes + 1
            else:
                m1 = 2 * nmodes + 1
            for j in range(i, m1):
                for k in range(j, 2 * nmodes + 1):
                    if k <= nmodes:
                        m2 = nmodes + 1
                    else:
                        m2 = 2 * nmodes + 1
                    for l in range(k, m2):
                        if valid_4prod(i, j, k, l):
                            r = pauli_weight(tt, [i, j, k, l])
                            if r > 0:

                                s += r
                                n += 1
    print(n, ", ", s, ", ", s / n)
    return [s / n, n]

def perm_check():
    n = 0
    k = 0
    l = 1
    for i in range(1, 9):
        l = l * i
    TT = [0] * l
    for nmodes in range(4, 5):
        tt = TernaryTree(2 * nmodes)
        tt.build_alphabeta_tree(2 * nmodes)
        print("tt", end=' : ')
        for tr in itertools.permutations(range(1, 2 * nmodes + 1)):
            print("tr = ", tr, end=":  ")
            tt.update_branchnum(list(tr) + list(range(2 * nmodes + 1, 4 * nmodes + 1)))
            TT[n] = main(tt, only2prod=True)
            n += 1
            if n / l > 0.01 * k:
                print(k, end=' ')
                k += 1

# dr = True
if __name__ == "__main__":
    # tt = TernaryTree()
    # tt.build_bk_tree(12)
    # draw(tt)
    BK = []
    TT2 = []
    TT1 = []
    JW = []
    N = range(2, 20)
    for nmodes in N:
        tt = TernaryTree(2 * nmodes)
        print("jw",end=' : ')
        JW = JW + [main(tt)[0]]

        tt.build_alphabeta_tree(2 * nmodes)
        print("tt", end=' : ')
        tt.update_branchnum(st_enumeration5(2*nmodes))
        prod2 = main(tt, only2prod=True)
        tt.update_branchnum(st_enumeration4(2 * nmodes))
        prod4 = main(tt, only4prod=True)
        TT2 = TT2 + [(prod2[0]*prod2[1] + prod4[0]*prod4[1])/(prod2[1] + prod4[1])]

        tt.build_alphabeta_tree(2 * nmodes)
        print("tt4", end=' : ')
        tt.update_branchnum(st_enumeration3(2 * nmodes))
        TT1 = TT1 + [main(tt)[0]]

        print("bk", end=' : ')
        tt.build_bk_tree(2 * nmodes)
        BK = BK + [main(tt)[0]]
    print("JW = ", JW)
    print("BK = ", BK)
    # print("TT1 = ", TT1)
    print("TT2 = ", TT2)
    # print("TT = ", TT)
    graph(data=[JW,BK,TT2,TT1], label=['JW mapping', 'BK mapping', 'TT mapping', 'TT1 mapping'], N=N)
