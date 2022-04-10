function retVal = J1_RootDivRad(num, radius)
    J1Roots = [3.83171
            7.01559
            10.17347
            13.32369
            16.47063
            19.61586
            22.76008
            25.90367
            29.04683
            32.18968
            35.33231
            38.47477
            41.61709
            44.75932
            47.90146
            51.04354
            54.18555
            57.32753
            60.76946
            63.61136
            66.75322
            69.89507
            73.03689
            76.178699
            79.32048
            82.46225
            85.60401
            88.74576
            91.88750
            95.02923];
    retVal = J1Roots(1:num)' ./ radius;
