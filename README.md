## Code for the slim simulations in the preprint:
Femerling, G., Van Oosterhout, C., Feng, S., Bristol, R., Zhang, G., Groombridge, J., ... & Morales, H. E. (2022). Genetic load and adaptive potential of a recovered avian species that narrowly avoided extinction. bioRxiv.

### Dependencies:
- SLiM3 (https://messerlab.org/slim/)
- Parallel for bash

### test run
Small run of 10 replicates for a population N=500 in parallel 2 computing cores
```./SPF_slim_parallel_run.sh 500 1_out out_test 0.2 0.25 2 10```
