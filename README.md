# hackthew3

A proposal for 'hack the W3' hackaton. The repository is intended to be used as the workplace for the initiative.

# Aim
The aim of the initiative is to replace the calls of the bottom layers (modular arithmetic) of the blst library by bolos (Ledger Nano cryptographic device).

# How we proceed ?
- A first level of integration is to brutally integrate the blst on Nano with a 'ZKP-boilerApp' demonstrating the PoC
- Test vectors should be extracted from here and tested for compliance : https://gitlab.inria.fr/zk-curves/snark-2-chains/-/tree/master/sage
- the montgomery functions are the primary target for bolos replacement (look for all mont_* functions and replace them by bolos equivalent)
- Tis' a bit challenging for a 2 day hackaton but it will be necessary at some point, let the shiny monkeys to others.

# License

The files in this project are provided under a dual BSD/GPLv2 license. When
using or redistributing this software, you may do so under either license.




