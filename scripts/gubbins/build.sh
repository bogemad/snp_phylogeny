export CFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"
sed -i 's~setup.py install~setup.py install --single-version-externally-managed --record=/tmp/record.txt~g' $SRC_DIR/python/Makefile.am

autoreconf -i
./configure --prefix=$PREFIX
make
make install
cd python
$PYTHON setup.py install --single-version-externally-managed --record=/tmp/record.txt
