rm -rf print_testres.txt
./int2-test-helmrouts3d
./int2-test-hfmm3d
./int2-test-hfmm3d-zkbig
./int2-test-hfmm3d-scale
./int2-test-hfmm3d-vec
./int2-test-hfmm3d-mps
mv print_testres.txt ../../print_testreshelm.txt
rm -rf fort.13
