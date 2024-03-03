# CppAD & IPOPT Install Tips



## 1. Eigen

```bash
sudo apt-get install libeigen3-dev
```



## 2. fmt

```bash
sudo apt-get install libfmt-dev
```



## 3. CppAD

```bash
sudo apt-get install cppad
```



## 4. IPOPT

```bash
sudo apt-get install -y gfortran libncurses5-dev libncursesw5-dev coinor-libipopt-dev libmetis-dev
```

To avoid strange compilation errors after installation, we need to make appropriate modifications in `/usr/include/coin/IpSmartPtr.hpp`, like this :

```c
#define HAVE_CSTDDEF	// add code
#ifdef  HAVE_CSTDDEF
#  include <cstddef>
#else
#  ifdef HAVE_STDDEF_H
#    include <stddef.h>
#  else
#    error "don't have header file for stddef"
#  endif
#endif
#undef HAVE_CSTDDEF 	// add code
```



## 5. osqp

```bash
git clone --recursive https://github.com/osqp/osqp
cd osqp
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON
make
sudo make install
```



## 6. osqp-eigen

```bash
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build 
cd build
cmake ..
make
sudo make install
```

