#time ./bin/PGunAnalysis.exe cfg/noiseScenario_0_algo_0_scaleByArea_False_scaleByDose_False.cfg >& /tmp/deguio/1.txt &
time ./bin/PGunAnalysis.exe cfg/20190330/noiseScenario_0_algo_1_scaleByArea_False.cfg >& /tmp/deguio/2.txt &
time ./bin/PGunAnalysis.exe cfg/20190330/noiseScenario_0_algo_2_scaleByArea_False.cfg >& /tmp/deguio/3.txt &
time ./bin/PGunAnalysis.exe cfg/20190330/noiseScenario_0_algo_2_scaleByArea_True.cfg  >& /tmp/deguio/4.txt &

time ./bin/PGunAnalysis.exe cfg/20190330/noiseScenario_3000_algo_2_scaleByArea_False.cfg >& /tmp/deguio/5.txt &
time ./bin/PGunAnalysis.exe cfg/20190330/noiseScenario_3000_algo_2_scaleByArea_True.cfg >& /tmp/deguio/6.txt &
