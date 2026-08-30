# badge-genetics-sim

A browser simulation of the DEF CON 34 badge's population-genetics system:
10,000+ virtual badges evolving in your tab, plus a widget that can talk to
your physical DC34 badge over the QR exchange.

Live at **[dc34sim.com](https://dc34sim.com)**.

Open [`index.html`](index.html) in a browser. No install, no build step, no
camera required.

> **Want the standalone, browser-native static HTML page so you can tinker with it and run it locally?** It lives on the `gh-pages` branch — [view it on GitHub](https://github.com/charlesreid1/dc34sim/blob/gh-pages/index.html), or grab just that branch with:
>
> ```
> git clone -b gh-pages --single-branch https://github.com/charlesreid1/dc34sim.git
> python3 -m http.server 8000
> 
> **Then visit `http://localhost:8000`.**

## Tabs

Five tabs, each a different lens on the same
genetics core:

- **popsim** - evolve a whole population of badges. Click any cell to inspect it.
- **vim gene** - edit one badge's genome (a diploid: 2 x 9 bytes) and mint QRs a real badge will accept.
- **skeet** - take one diploid and enumerate every gamete it can produce (all `2^5 = 32`), then pick one and seal it into a QR for a real badge.
- **borg** - same as popsim, but every mating uses one frozen partner, ROOT.
- **badgecestry** - decode one badge's ancestry composition with a two-stage hidden Markov model.

Any specimen inspector can send its selected badge to **vim gene**, **skeet**,
or **badgecestry**, so the tabs compose into a pipeline: evolve a population,
pull out an individual, edit it, enumerate its gametes, seal one, then hand
the QR to a physical badge.

## Skeet in one paragraph

Skeet is meiosis-as-a-truth-table. One diploid in, all 32 possible gametes
out (5 linkage groups, fair independent coins, exact enumeration - no
sampling). Each row shows the 5-bit coin pattern, 9 locus cells colored by
source haploid so the linkage blocks are visible, and a `+ Pick` button.
Picking a gamete reveals a QR gene-exchange panel that mirrors vim gene's
three-phase flow: paste your badge's nonce, seal the picked gamete under
that nonce, and hand back a QR your badge will accept as a mate. See the
[`skeet.md`](skeet.md) and the [FAQ's Skeet section](faq.md#skeet).

## Docs

- [`faq.md`](faq.md) - user-facing FAQ, keyed by tab.
- [`synthetic-population-genetics.md`](synthetic-population-genetics.md) - the genetics core (haploids, loci, mutation, meiosis, mating).
- [`gene-exchange.md`](gene-exchange.md) - the three-phase QR protocol for talking to a physical badge.
- [`borg-genetics.md`](borg-genetics.md) - the ROOT-partner colony variant.
- [`badgecestry.md`](badgecestry.md) - the two-stage HMM used for ancestry decoding.
- [`skeet.md`](skeet.md) - the gamete-enumeration tab: genetics primer + physical-badge exchange.

## Related repositories

- [`bunnie/dc34-api`](https://github.com/bunnie/dc34-api) - genetics core: `Haploid`, `Diploid`, `BadgeType`, `phenotype`, `meiosis`, `mutate`.
- [`bunnie/dc34-vault`](https://github.com/bunnie/dc34-vault) - exchange state machine and the k0 oracle.
- [`bunnie/dc34-console`](https://github.com/bunnie/dc34-console) - LED renderer.
- [`nastea1/dc34-gamete`](https://github.com/nastea1/dc34-gamete) - QR wire format, published k0, and the reference workbench.
- [U.S. Patent 12,626,778](https://patents.justia.com/patent/12626778) - "Accelerated hidden Markov models for genotype analysis".

## Acknowledgements

- Huge thanks to [Bunnie Huang](https://www.bunniestudios.com/) - creator of the DC34 badge and chip - and to the entire DEF CON 34 badge team, without whom none of this would exist.
- **badgecestry** - the two-stage HMM ancestry decoder draws on techniques described in [Ancestry.com](https://www.ancestry.com/)'s [U.S. Patent 12,626,778](https://patents.justia.com/patent/12626778), "Accelerated hidden Markov models for genotype analysis".

## License

See [`LICENSE`](LICENSE).
