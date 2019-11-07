g++ -O3 test_mystructs.cpp -std=c++14
valgrind --tool=massif --time-unit=ms ./a.out 28000 20
massif-visualizer massif.out.*
