### About
A MCMC-based sampler to sample directed graphs G such that
- the underlying graph pr(G) is retained
- the number of directed edges is retained
- the number of simplices in its directed flag complex stays within given borders.


See TODO:ADD_PAPER for details.

### Usage:
Install rust (e.g. via rustup.rs) an type

`cargo build --release --bin sample`

then run the code, e.g. via

`target/release/sample -i input.flag -l experiment_label` 

for other features, call

`target/release/sample --help`

### Roadmap:
- Fix Doku.



