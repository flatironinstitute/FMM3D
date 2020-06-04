rm -rf print_testres.txt
./int2-test-laprouts3d
./int2-test-lfmm3d
./int2-test-lfmm3d-scale
./int2-test-lfmm3d-vec
mv print_testres.txt ../../print_testreslap.txt
rm -rf fort.13
