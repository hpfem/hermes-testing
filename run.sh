cmake .
echo "This test suite will fail if you do not have Hermes, or Valgrind, etc.."
echo "CMake generation done."
echo "Run memory-leaks tests? (Will take some time) [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing memory leaks tests..."
  cd memory-leaks
  cd 01-simple
  echo "Building the test..."
  make
  echo "Test built, running..."
  valgrind --leak-check=full --log-file=../01-simple-ValgrindLogfile ./01-simple 3 3
  echo "Valgrind output '01-simple-ValgrindLogfile' available in memory-leaks/"
  cd ../02-adapt
  echo "Building the test..."
  make
  echo "Test built, running..."
  valgrind --leak-check=full --log-file=../02-adapt-ValgrindLogfile ./02-adapt
  echo "Valgrind output '02-adapt-ValgrindLogfile' available in memory-leaks/"
  echo "Memory leaks tests - Done."
  cd ../..
fi
echo "Run visualization tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing visualization tests..."
  cd visualization/views
  make
  ./01-all
  echo "Outputs saved in visualization/views/*.bmp"
  echo "Visualization tests - Done."
  cd ../../..
fi
echo "Run inner-funcionality tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing internal tests..."
  cd inner-functionality
  make -j4
  ctest -j4
  echo "Internal-functionality tests - Done."
  cd ..
fi
echo "Run calculations tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing calculations tests..."
  cd calculations
  make -j4
  ctest -j4
  echo "Calculations tests - Done."

fi
echo "Quitting..."
