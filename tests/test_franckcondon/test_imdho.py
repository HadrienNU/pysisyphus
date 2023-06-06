import numpy as np

from pysisyphus.franckcondon.imdho import imdho_abs_cross_section


def test_imdho():
    gamma = 160
    dE_exc = 36024

    displs = np.array(
        """1     0.000860
    2     0.000101
    3    -0.000042
    4    -0.000531
    5    -0.411010
    6    -0.568713
    7     0.000058
    8     0.000000
    9    -0.001358
    10     0.000009
    11    -0.000231
    12     0.000018
    13     0.046431
    14    -0.000829
    15     0.000006
    16    -0.000015
    17    -0.000147
    18     0.000041
    19     0.503103
    20    -0.000029
    21    -0.516329
    22    -0.000100
    23     0.066056
    24    -0.200802
    25     0.000079
    26    -0.181164
    27    -1.107052
    28    -0.000812
    29    -0.026333
    30     0.000026
    31    -0.000978
    32    -0.000009
    33     0.015524
    34     0.058773
    35    -0.015335
    36    -0.001115""".split(),
        dtype=float,
    ).reshape(-1, 2)[:, 1]

    nus = np.array(
        """1 100.985227
    2 144.001276
    3 229.749397
    4 251.485668
    5 341.194669
    6 436.280333
    7 529.543611
    8 597.376354
    9 685.937137
    10 865.855734
    11 880.735871
    12 887.407602
    13 922.814195
    14 926.641605
    15 957.477214
    16 985.65681
    17 1019.0173
    18 1136.567469
    19 1191.786311
    20 1245.073446
    21 1274.295136
    22 1287.088528
    23 1289.451551
    24 1392.306925
    25 1424.365418
    26 1578.612388
    27 1624.373846
    28 1627.149545
    29 3043.639615
    30 3046.091489
    31 3055.47366
    32 3060.29147
    33 3066.032777
    34 3066.058264
    35 3155.342259
    36 3155.353628""".split(),
        dtype=float,
    ).reshape(-1, 2)[:, 1]

    dEs_inc = np.linspace(dE_exc - 1000, dE_exc + 6000, num=50)
    cross_secs = imdho_abs_cross_section(dEs_inc, dE_exc, gamma, displs, nus)

    # import matplotlib.pyplot as plt
    # plt.plot(dEs_inc, cross_secs)
    # plt.show()
