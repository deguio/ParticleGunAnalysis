#! /usr/bin/env python
import os
import subprocess


RADLOWRANGE  = ["1", "9", "17", "25", "33", "41"]
RADHIGHRANGE = ["8", "16", "24", "32", "40", "48"]

LAYVALUES    = ["15", "22"]


for lay in LAYVALUES:
    for low, high in zip(RADLOWRANGE, RADHIGHRANGE):
        os.system("root -l -q \'/Applications/root_init/setTDRStyle.C\' \'macros/plotCumulativeRDF.C(\"HGCDigiADC[HGCDigiIndex==2 && HGCDigiIEta < 0 && HGCDigiLayer=="+lay+" && abs(HGCDigiIEta)>="+low+" && abs(HGCDigiIEta)<="+high+"]\",64)\'")
