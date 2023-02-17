# cgl_rs
An implementation of the CGL hash function in rust (using modular polynomials). Warning: this is my first project in rust so it is not perfect!

See https://eprint.iacr.org/2006/021.pdf for the original paper on the CGL hash function.

See a pseudo code snippet below for the general idea of the algorithm.

![cgl](https://user-images.githubusercontent.com/17739301/219529234-55e23af6-97fd-4dfd-a646-c32eb8b3ace7.png)


Note that the sage version does not implement the same square root algorithm and so they output different results, this will be fixed in future version.

Requires sagemath and rust, simply run ./benchmark.sh to test for 256 bit input.
