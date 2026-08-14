# dc34sim

A browser-native Monte Carlo simulator of the DEF CON 34 badge genetics system.

Live at **[dc34sim.com](https://dc34sim.com)**.

> **Want the standalone, browser-native static HTML page so you can tinker with it and run it locally?** It lives on the `gh-pages` branch. Grab just that branch with:
>
> ```
> git clone -b gh-pages --single-branch https://github.com/charlesreid1/dc34sim.git
> ```

The DC34 badge runs a population-genetics simulator under the hood: diploid genome, meiosis with independent assortment, tunable mutation, an assortative-mating penalty, and eight badge types with distinct starting allele ranges. The physical badge is one node; the "mate choice" function is two humans pointing badges at each other. That works fine for a weekend, but it doesn't scale to the questions the model actually invites — the population is too small, mating is too sparse, and one weekend is not enough generations to see equilibria.

This repo pulls the mathematical core out of the badge firmware, ports it to JavaScript, and runs it at 10,000 individuals in your browser. No dependencies, no build step, no server. Open the page and watch the population evolve.

## What's in here

- **[`html/index.html`](html/index.html)** — the simulator. Single file, client-side only. Renders each of 10,000 badges as a 5×5 animated cell in a 100×100 grid. Seeded RNG so runs are reproducible. Includes the panmictic simulator, a "vim gene" specimen editor, and a **borg** tab that pins every mating to a fixed ROOT genotype.
- **[`synthetic-population-genetics.md`](synthetic-population-genetics.md)** — the writeup. Genome layout, meiosis and its five linkage groups, the phenotype expression map (including the asymmetric `nonlin` line), Gray-code mutation, the inbreeding pass, and what to measure. Start here if you want to understand what the sim is doing.
- **[`borg-genetics.md`](borg-genetics.md)** — supplement covering the borg regime. What happens when every colony member's slot-1 chromosome is a fresh gamete drawn from a frozen ROOT diploid? Predictions and experiments.

## Running it

There is no build step. Either open `html/index.html` in a browser directly, or serve the directory statically:

```
python3 -m http.server -d html 8000
```

Then visit `http://localhost:8000`.

## Provenance

The genetics core (`meiosis`, `mutate`, `phenotype`, badge-type ranges, Gray-code LUTs, mutation-rate tables) is a JavaScript port of the Rust in `dc34-api/src/lib.rs` and `dc34-vault/src/`. All ranges, expressions, and probabilities are lifted directly. Everything labeled *policy* in the docs (selection regimes, mating schedulers, the borg pairing rule) is a simulator extension — the physical badge imposes no such constraints on the mathematical core.

See [`synthetic-population-genetics.md` §12](synthetic-population-genetics.md) for exact file/line citations.

## License

See [`LICENSE`](LICENSE).
