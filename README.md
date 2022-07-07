# Results (6/7/2022)

- blst has been integrated in the speculos framework without tests (not tested= it fails)
- no integration or test where performed over Nano
- final presentation is available in final_talk repo

We failed fast and had fun. Congrats to the 'Ledger Your crypto' winner team. The work now continues in the cylib repo.


# (Expired) 
What follows was the initial objectives of the hackaton
## Aim
The aim of the initiative is to replace the calls of the bottom layers (modular arithmetic) of the blst library by bolos (Ledger Nano cryptographic device).

## How we proceed ?
- A first level of integration is to brutally integrate the blst on Nano with a 'ZKP-boilerApp' demonstrating the PoC
- Test vectors should be extracted from here and tested for compliance : https://gitlab.inria.fr/zk-curves/snark-2-chains/-/tree/master/sage
- the montgomery functions are the primary target for bolos replacement (look for all mont_* functions and replace them by bolos equivalent)
- Tis' a bit challenging for a 2 day hackaton but it will be necessary at some point, let the shiny monkeys to others.

## License

The files in this project are provided under a dual BSD/GPLv2 license.
