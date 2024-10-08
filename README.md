# Transferable E-cash

A rust implementation of Transferable e-cash according to the paper [Anonymous Transferable E-Cash](https://www.iacr.org/archive/pkc2015/90200211/90200211.pdf). In this library, the code is mainly* written based on the paper [Transferable E-cash: A Cleaner Model and the First Practical Instantiation](https://eprint.iacr.org/2020/1400) which provides concret construction of the scheme.

****some minor parts have been modified/TBC***

> Transferable e-cash is the digital analog of physical cash, which allows users to transfer coins between them in isolation without interacting with a bank or a “ledger”. 
It protects user privacy and, provides means to trace fraudulent behavior (double-spending of coins).

The entire work consists of some building blocks. The current progress towards their completeness:

- [x] Double Spending Tag
- [ ] Encryption Scheme E' (Replayable-CCA encryption scheme)
    - [X] One-time linearly homomorphic structure-preserving signature
- [ ] Encryption Scheme E
- [ ] Transferable E-cash construction


This work is built on top of [Arkworks](https://github.com/arkworks-rs/), a rust ecosystem for cryptographic libraries such as elliptic curve arithmetic.

Note: This library has not been thoroughly audited for production use. Please take your own risk to use it in production.


## Reference:

- Baldimtsi, F., Chase, M., Fuchsbauer, G., Kohlweiss, M.: Anonymous transferable e-cash. (2015).
- Balthazar Bauer, Georg Fuchsbauer, and Chen Qian: Transferable E-cash: A Cleaner Model and the First Practical Instantiation. (2020).