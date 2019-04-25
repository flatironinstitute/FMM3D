rm -rf print_testres.txt
./test_helmrouts3d
./test_hfmm3dpart
./test_hfmm3dpart_vec
mv print_testres.txt ../../print_testreshelm.txt
rm -rf fort.13
