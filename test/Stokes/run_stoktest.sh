rm -rf print_testres.txt
./int2-test-stokkernels
./int2-test-stfmm3d
mv print_testres.txt ../../print_testresstok.txt
rm -rf fort.13
