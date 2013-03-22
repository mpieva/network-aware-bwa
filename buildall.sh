set -e
version=$( awk '/#define PACKAGE_VERSION/ { print $3 }' main.h | tr -d '"')
version="${version}-`git describe --always`"
echo "Building and installing bwa-${version}"

make clean
make
make install prefix=/home/public/usr/stow/bwa-${version}

cd /home/public/usr/stow
stow -v -D bwa*
stow -v    bwa-${version}
