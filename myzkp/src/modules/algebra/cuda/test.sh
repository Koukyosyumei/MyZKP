mkdir build

nvcc -std=c++17 test_fr.cu -o build/test_fr
nvcc -std=c++17 test_sumcheck.cu -o build/test_sumcheck
./build/test_fr
./build/test_sumcheck

rm -rf build
