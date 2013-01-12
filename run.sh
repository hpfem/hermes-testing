cmake .
echo "This test suite will fail if you do not have Hermes, or Valgrind, etc.."
echo "CMake generation done."
echo "Run memory-leaks tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing memory leaks tests..."
  cd memory-leaks
  cd 01-simple
  make
  valgrind --leak-check=full --log-file=../01-simple-ValgrindLogfile ./01-simple 3 3
  echo "Valgrind output '01-simple-ValgrindLogfile' available in memory-leaks/"
  cd ../02-adapt
  make
  valgrind --leak-check=full --log-file=../02-adapt-ValgrindLogfile ./02-adapt
  echo "Valgrind output '02-adapt-ValgrindLogfile' available in memory-leaks/"
  echo "Memory leaks tests - Done."
fi
echo "Run visualization tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing visualization tests..."
  cd visualization/views
  cd 01-all
  make
  ./01-all
  echo "Outputs saved in visualization/views/*.bmp"
  echo "Visualization tests - Done."
fi
echo "Quitting"
