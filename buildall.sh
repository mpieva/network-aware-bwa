set -e
version=$( awk '/#define PACKAGE_VERSION/ { print $3 }' main.h | tr -d '"')
echo "Building and installing bwa-${version}"

make clean
make
make install prefix=/home/public/usr64/stow/bwa-${version}

cd /home/public/usr64/stow
stow -v -D bwa*
stow -v    bwa-${version}
