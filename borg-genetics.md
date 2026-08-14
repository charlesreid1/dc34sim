# BORG GENETICS

What if an organism acquired ROOT access in the DEF CON badge colony?
Once ROOT hacks the mating process, nobody in the colony will ever have another mate.
Every mating event always involves ROOT. Without mutations, the genome of
ROOT will eventually become the genome of the entire colony. The irresistable
pull of genetic gravity.

Can a colony of badges mutate fast enough to overcome the pull of ROOT? 
What would escape look like? How does the shooting star phenome distribution
change, if ROOT has it?

## What this is

This is a supplement to `synthetic-population-genetics.md`. Read that first. Everything in §§0-6 of that
document — genome layout, meiosis, phenotype expression (including the `nonlin` asymmetry), mutation, and the
Gray-code trick — still applies.

We just add one twist.

The **borg** tab in [`index.html`](https://github.com/charlesreid1/dc34sim/blob/gh-pages/index.html) runs the same genetic simulator, but replaces panmictic random
pairing with a fixed second partner in every mating: **ROOT**. ROOT is itself a diploid — two chromosomes,
same nine loci as anyone else — and every mating event in the borg regime is `(some colony member, ROOT)`.
Every colony member is a diploid too; that never changes. What changes is the pairing rule and the fact
that ROOT's two chromosomes are frozen for the entire run: never mutated, never overwritten, never replaced.

The whole regime is a special case of `synthetic-population-genetics.md` §7, with the mating scheduler
collapsed to a constant function (partner = ROOT, always) and a small handful of side conditions turned off.
But the *behavior* under that collapse is qualitatively different enough to deserve its own page.

---

## 1. What Actually Changes

Both the colony members and ROOT are diploids — nine loci × two chromosomes each. Every mating produces a
child diploid whose two chromosomes are: **slot 0** = a fresh meiotic gamete from the colony member's own
diploid, **slot 1** = a fresh meiotic gamete from ROOT's diploid. Slot 0 replaces that colony member in
place; ROOT is not replaced at all. So a generation is `N` mating events, one per colony member, each
producing one child diploid that overwrites its parent.

That's it. The vocabulary from popsim (`_egg`, `_sperm` scratch buffers) still shows up in the shipped code
because the same primitives are reused, but here the names are just labels for "slot 0 gamete" and "slot 1
gamete" — no biological sex, no gendered roles.

Compared to the panmictic sim (`synthetic-population-genetics.md` §7):

| step                             | popsim (panmictic)                            | borg (ROOT-paired)                                  |
|----------------------------------|-----------------------------------------------|-----------------------------------------------------|
| partner selection                | Fisher-Yates shuffle each generation          | none — every member is paired with ROOT             |
| slot-0 gamete source             | `meiosis(recipient diploid)`                  | `meiosis(self diploid)` — the member's own diploid  |
| slot-0 gamete mutation           | `mutate(...)` at baseRate (elev if inbreed)   | `mutate(...)` at baseRate                           |
| slot-1 gamete source             | `meiosis(donor diploid)`                      | `meiosis(rootDiploid)` — from the fixed ROOT        |
| slot-1 gamete mutation           | `mutate(...)` at baseRate (+inbreed pass)     | **not applied** — the ROOT gamete is written raw    |
| inbreeding pass (same-type)      | extra mutation on both slots at floor rate    | **disabled entirely** (no `inbreedingMode` control) |
| child badge type                 | inherited from slot-0 (recipient) member      | inherited from the colony member (unchanged)        |
| slot layout of the child diploid | `[slot 0, slot 1]` (matches `replace_gene()`) | `[slot 0, slot 1]` — unchanged                      |
| ROOT diploid itself              | n/a                                           | **frozen** — never mutates, never recombines        |

Four things worth considering:

1. **ROOT is diploid, and so is every colony member.** Nobody in this simulation is haploid. The word
   "haploid" only ever applies to a *gamete* — the transient 9-byte output of `meiosis()` that lives on the
   scratch buffer for the duration of one mating and is then written into a slot of some diploid. Every
   individual — ROOT included — is 2 × 9 bytes at rest. What makes ROOT "frozen" is that its 18 bytes are
   never mutated and never overwritten; the colony members' 18 bytes are overwritten every generation.
   Meiosis still runs on ROOT every mating, so a mating draws one of `2^5 = 32` possible ROOT gametes (five
   independent linkage-group coin flips), not the same 9 bytes every time.

2. **The ROOT-side gamete is not mutated.** In popsim, both gametes going into the child diploid get a
   `mutate()` pass (§5.1 of the popsim doc), and on same-type matings the inbreeding branch adds another
   pass. In borg, only the colony member's own gamete is mutated. The ROOT gamete is written into slot 1
   raw. All generation-over-generation variability on the slot-1 side comes entirely from the 5-coin meiotic
   reshuffle of ROOT's two chromosomes. This is what keeps ROOT alleles recognizable in the population even
   after hundreds of generations: they can be reshuffled but not eroded.

3. **No inbreeding pass.** The popsim inbreeding branch fires only when the two mating individuals
   have equal badge types. ROOT doesn't participate in that comparison — the borg regime skips the branch
   unconditionally, so there is no `Elevated`-floor elevation, no double mutation, no elevated slot-0
   mutation. The UI has no inbreeding toggle on the borg tab for exactly this reason.

4. **Slot-0 mutation is the only source of new alleles.** Every generation the slot-0 side does one meiotic
   reshuffle of the colony member's own two chromosomes, followed by one `mutate()` pass at `baseRate`.
   That is the entire entropy budget for the run. If you set `mutationRate=None`, slot 0 becomes a pure
   shuffle of the previous generation's slot-0 pool, which drifts under classic Wright-Fisher dynamics —
   except that slot 1 is *always* pinned to a ROOT gamete, so half of every child diploid is drawn from a
   two-allele distribution regardless.

Everything else (linkage groups, Gray-code mutation, additive-dominant phenotypes, the
`nonlin` asymmetry, which still fires on the `[slot 0, slot 1]` child in slot order) is
identical to popsim.

---

## 2. Why The Default Rate Is Higher

The borg tab defaults to **Radioactive** (54.7% per-locus mutation), not Baseline. This is a policy choice,
not a code change: `RATE_VALUE` and `RATE_BITS` are the same tables as popsim, but the default `<select>` in
the borg control panel starts one step further up the ladder.

The reason is entropy budget. In popsim there are three independent `mutate()` calls per mating event:
slot-1 gamete mutation, slot-0 gamete mutation, and (on same-type) the inbreeding-elevated extra passes on
both slots. In borg there is exactly **one** — slot-0 gamete mutation. The slot-1 (ROOT) gamete is written
raw. So the mutational "budget" per generation drops by a factor of two to three at the same rate label, and
drift toward the ROOT phenotype accelerates because there is less noise to push the population away from what
ROOT keeps delivering. Radioactive restores enough per-locus churn on the slot-0 side to keep the
allele-frequency spectra visibly evolving instead of collapsing onto ROOT-flavored equilibria in a handful of
generations.

If you want a quiet, "ROOT takes full control" run, set the rate to Baseline or None. If you want to see how
much slot-0 mutation is needed to *escape* ROOT's pull on a given locus, sweep the rate and watch how far the
allele-frequency spectra drift from ROOT's own genotype. The `nonlin` phenotype is a particularly satisfying
one to watch because the shipped asymmetric expression reads slot-0 `chaser` + slot-1 `nonlin` — and slot 1
is always a ROOT gamete here, so every displayed `nonlin` in the population inherits a fixed additive
contribution from ROOT regardless of what the colony member's own alleles say.

---

## 3. Where ROOT Comes From

ROOT is snapshotted from the **vim gene** tab at reset time. The borg tab has no ROOT editor of its own; when
you press `reset` on the borg controls, `initPopulation()` reads the current 2×9 bytes out of the vim gene
editor (`vimGene.chrom0` / `vimGene.chrom1`) and hands them to `BorgPopulation` as `params.rootDiploid`. The
constructor takes an immutable `Uint8Array` copy — subsequent edits in vim gene don't reach into a running
borg population. Reset the borg tab to install a new ROOT.

This lets you construct a specific ROOT genotype in vim gene (e.g. a rare shooting-star Uber chaser, or a
deliberately gap-filling `hue_base=25` allele), then watch what happens to a 10,000-member population when
every child gets one of ROOT's gametes as its slot-1 chromosome.

The initial population is generated the normal way (§2 of the popsim doc): each member gets a badge type
drawn from the mix preset, and its two chromosomes are freshly generated from that type's ranges. ROOT does
not participate in that step — ROOT only shows up at mating time in `step()`. So generation 0 is a "clean"
per-type population, and each subsequent generation replaces the slot-1 chromosome of every member with a
fresh ROOT gamete.

---

## 4. Population Genetics Under A Frozen Partner

Some quick predictions to test against the sim:

1. **Slot-1 is fixed, immediately.** After one generation, every colony member's slot 1 is a fresh ROOT
   gamete. The slot-1 allele frequency spectrum, per locus, is exactly `0.5 * ROOT.chrom0[locus] + 0.5 *
   ROOT.chrom1[locus]` (two possible values per locus in groups B, D, E; a correlated pair of two possible
   values per locus in groups A and C).

2. **Slot-0 is drift plus mutation.** The slot-0 side never sees an outside allele, so its `N` haploid
   slots (one per colony member) are just the previous generation's slot-0 pool reshuffled and mutated.
   Under `mutationRate=None`, this is textbook Wright-Fisher drift with `N_e = N`. Under Baseline+, the
   mutation kernel pushes it toward a rate-dependent stationary distribution.

3. **Additive-dominant loci converge on ROOT fast.** `sat`, `cd_dir`, `chaser` are all `sat_add`-dominant.
   Because every member gets one of ROOT's alleles at those loci in slot 1, any member whose slot-0 allele
   plus ROOT's slot-1 allele exceeds 255 saturates to 255. If ROOT carries a high allele at any of these, the
   whole population's phenotype at that locus goes near-255 immediately. If ROOT carries a low allele, the
   population *can* still saturate but only if the slot-0 side drifts high — a good indirect readout of how
   much slot-0 entropy is present.

4. **`nonlin` phenotype has a hard floor from ROOT.** The shipped asymmetric expression reads slot-0
   `chaser` + slot-1 `nonlin`. Slot 1 is always a ROOT gamete, so `slot-1.nonlin` is one of ROOT's two nonlin
   alleles every mating. `slot-0.chaser` is whatever the colony member's own slot-0 gamete brought. So
   displayed `nonlin` across the population is `sat_add(slot-0.chaser, ROOT_gamete.nonlin)` — a fixed
   additive contribution from ROOT overlaid on the population's own slot-0 chaser distribution. Under
   symmetric mode, it becomes `sat_add(slot-0.nonlin, ROOT_gamete.nonlin)` — the slot-0 nonlin drift shows
   up cleanly against a fixed ROOT baseline. Diffing the two modes here is the crispest possible instance of
   the "quantify the nonlin bug" experiment (`synthetic-population-genetics.md` §9.7), because the slot-1
   side is a constant instead of a moving target.

5. **Badge-type distribution is invariant.** Every child inherits the colony member's own badge type
   (`_nextTypes[i] = this.types[i]`). ROOT is not part of the badge-type accounting, so the type histogram is
   frozen at the initial mix. No demographic drift, no colonization by ROOT's type. This is a deliberate
   choice — badge type in the shipped code is a hardware property, not genetic — but it means the
   type-composition chart in the borg tab is decorative rather than diagnostic. Watch the allele spectra
   instead.

---

## 5. What ROOT Is Not

- **Not a fitness function.** ROOT does not KEEP or REVERT anyone. Every colony member is unconditionally
  replaced every generation. Selection is still exogenous and, in the current implementation, absent — this
  is neutral drift with a frozen slot-1 distribution. Layering a selection policy on top would be a
  one-function extension (see `synthetic-population-genetics.md` §8), and would compose cleanly with the borg
  mating regime because the two concerns don't touch each other.

- **Not the shipped badge's behavior.** No DC34 badge implements this. ROOT is a simulator-only construct
  that answers the counterfactual "what if there were one fixed genotype that every badge in the room
  invariably mated with?" It's the badge-genetics equivalent of a founder effect turned into a mating policy.

- **Not immune to the `nonlin` asymmetry.** ROOT chromosome order matters in exactly the same way as popsim
  chromosome order matters. `rootDiploid[0]` (chrom0) and `rootDiploid[1]` (chrom1) are distinguishable at
  meiosis time (the E-group coin picks between them), but the `nonlin` phenotype only ever reads
  `slot-1.nonlin` from the child diploid, and slot 1 is always a ROOT gamete. So ROOT's two `nonlin` alleles
  appear in the population's `nonlin` phenotype with 50/50 frequency, and the `chaser`-leak that produces
  the asymmetry comes from the colony member's *own* slot-0, not from ROOT. If you swapped the slot
  assignment so ROOT wrote slot 0 and the colony member wrote slot 1, the picture would flip: ROOT's
  chaser would leak into every displayed nonlin. That's not implemented, but it would be a two-line edit in
  `BorgPopulation.step()`.

---

## 6. Implementation Delta

The whole regime is one class extending `Population` with a rewritten `step()`. Everything else is inherited.

Code below quoted verbatim from [`index.html`](https://github.com/charlesreid1/dc34sim/blob/gh-pages/index.html). Note that the shipped source names its two
scratch buffers `_egg` and `_sperm` — labels for "the gamete that ends up in slot 0" and "the gamete that ends up in slot 1".
There is no biological sex in the model.

```javascript
class BorgPopulation extends Population {
  constructor(N, seed, params) {
    super(N, seed, params);
    // Immutable copy — the ROOT diploid is fixed at construction time.
    this.rootDiploid = new Uint8Array(params.rootDiploid);
  }

  step() {
    const baseRate = this.params.mutationRate;
    for (let i = 0; i < this.N; i++) {
      // Slot-0 gamete: from the colony member's own diploid, mutated at baseRate.
      meiosis(this.pop, i * 2 * N_LOCI, this._egg, 0, this.rng);
      mutate(this._egg, 0, baseRate, this.rng);
      // Slot-1 gamete: fresh meiotic draw from ROOT. Not mutated.
      meiosis(this.rootDiploid, 0, this._sperm, 0, this.rng);
      // Child = [slot 0, slot 1], same slot order as popsim.
      const childOff = i * 2 * N_LOCI;
      for (let li = 0; li < N_LOCI; li++) {
        this._nextPop[childOff + li]          = this._egg[li];
        this._nextPop[childOff + N_LOCI + li] = this._sperm[li];
      }
      this._nextTypes[i] = this.types[i];
    }
    // Swap buffers, recompute phenotypes.
    const tp = this.pop; this.pop = this._nextPop; this._nextPop = tp;
    const tt = this.types; this.types = this._nextTypes; this._nextTypes = tt;
    this._computePhenotypes();
    this.generation++;
  }
}
```

Compare to `Population.step()`: gone are the Fisher-Yates shuffle, the reciprocal `_mate(iA, iB) + _mate(iB,
iA)` loop, the `inbreedOn` branch, the `maxRate(...)` floor calculation, the slot-1 mutation, and the
inbreeding-elevated slot-0 mutation. What's left is the minimal viable mating: one meiosis per side, one
mutation on the slot-0 gamete, one slot-ordered write of the child diploid.

The rest of the borg tab (histograms, specimen inspector, ROOT panel, seed input, mix preset) is plumbing
reuse — the `makeView()` factory lets popsim and borg share the same render + inspect + histogram code,
differing only in element ids and the `inbreedingMode: null` config that hides the inbreeding control.

---

## 7. What To Measure (Borg Edition)

Beyond the popsim measurements in `synthetic-population-genetics.md` §9, ROOT enables a few borg-specific
experiments:

1. **Distance from ROOT, per locus, over time.** For each locus, plot `|mean(pop[:, 1, locus]) - mean(ROOT[:,
   locus])|` (slot-1 side) and `|mean(pop[:, 0, locus]) - mean(ROOT[:, locus])|` (slot-0 side). Slot 1 snaps
   to ROOT's mean in one generation and stays pinned. Slot 0 shows the drift-vs-mutation battle explicitly.

2. **Mutation rate needed to escape ROOT's phenotype.** Sweep `mutationRate` from None to Apocalyptic and
   measure the fraction of the population whose displayed phenotype differs from ROOT's own phenotype by more
   than some threshold. This gives you the "escape velocity" for a given locus and quantifies how dominant a
   frozen slot-1 genotype actually is.

3. **ROOT's `nonlin` fingerprint.** Because slot-1 `nonlin` is always a ROOT allele, the population's `nonlin`
   phenotype distribution has a hard lower bound set by ROOT. Comparing the shipped asymmetric expression to
   the symmetric one under identical ROOT and RNG seeds isolates the `nonlin`-bug effect against a constant
   slot-1 baseline — cleaner than the popsim version of the same experiment, where the slot-1 side is itself
   moving.

4. **What happens to the initial per-type mix.** Because slot 1 is always a ROOT gamete and the population
   members' badge types never change, you can watch how e.g. a Human member's *displayed phenotype* moves
   toward or away from ROOT's phenotype as a function of ROOT's own badge type. Set ROOT to Uber (rare
   chaser, sparse nonlin) with a Human-heavy mix and count generations to phenotypic assimilation. That's
   your quantitative answer to "how fast does the ROOT genome take over a mixed population."

---

## 8. Source Provenance

Everything in this document is behavior of `class BorgPopulation` in [`index.html`](https://github.com/charlesreid1/dc34sim/blob/gh-pages/index.html) (§1
implementation, ~line 1371) plus the borg-specific pieces of `initPopulation()` (~line 2218) and the borg
view configuration (~line 2166). The genetics primitives (`meiosis`, `mutate`, `phenotype`, `haploidFromType`,
Gray-code LUTs, rate tables) are unchanged from popsim and are the browser port of the Rust genetics core
described in `synthetic-population-genetics.md` §12.

Anything in this document labeled as a *prediction* or an *experiment* is not implemented in the sim itself —
those are the natural next-step measurements the sim makes cheap.
