# SYNTHETIC POPULATION GENETICS

## Related repositories

- [`bunnie/dc34-api`](https://github.com/bunnie/dc34-api) - genetics core: `Haploid`, `Diploid`, `BadgeType`, `MutationRate`, `phenotype`, `meiosis`, `mutate`, `gray_encode`, `gray_decode`.
- [`bunnie/dc34-vault`](https://github.com/bunnie/dc34-vault) - exchange state machine and the k0 oracle (`sha256(k0)[..8]`).
- [`bunnie/dc34-console`](https://github.com/bunnie/dc34-console) - LED renderer.
- [`nastea1/dc34-gamete`](https://github.com/nastea1/dc34-gamete) - QR wire format (`PROTOCOL.md`), published `k0`, and a browser reference implementation (`gamete-workbench.html`).

## What The Badge Is Actually Doing

The DEF CON 34 badge runs a population-genetics simulator: diploid genome, meiosis with independent assortment,
mutation with a tunable rate, an assortative-mating penalty (the "inbreeding" pass), and a partitioned starting
population where each badge type has its own allele ranges. The physical badge is a distributed node running one
instance of this sim; the mate-choice function is two people pointing their badges at each other.

Stripped of the firmware layer, the mathematical core is a few hundred lines of vectorized array code. This
document pulls it out so it can be re-run at 10,000+ individuals with no hardware in the loop, and lets us ask
questions the physical conference can't answer at any reasonable statistical power:

- How fast does allele diversity accumulate under Baseline vs Apocalyptic?
- Under what conditions does the rare shooting star phenotype (`lin < 88`) fix, drift, or just vanish?
- Does the asymmetric `nonlin` phenotype expression actually change equilibrium behavior, or does it wash out in the population average?
- Does the assortative-mating penalty actually diverge isolated subpopulations, or does mutation swamp it?

---

## 0. Why This Works As A Model

The badge system is delightfully clean for a homebrew genetics toy. The badge team was certainly showing mercy by
building something this tractable:

- **Fixed genome length.** Nine one-byte loci. No indels, no crossover inside a gene, no ploidy changes.
- **Discrete allele space.** Each locus is a `u8`, exactly 256 possible alleles. The whole haploid space has cardinality `256^9 ≈ 4.7 × 10^21`, which sounds huge but decomposes cleanly per locus so it doesn't matter.
- **Independent assortment.** Meiosis picks each locus's parent haploid independently (see §3). No linkage, no intra-locus recombination.
- **Pure functions.** `phenotype()`, `meiosis()`, `mutate()` are all deterministic given their RNG stream. Trivially vectorizable.
- **No selection baked in.** The badge has no fitness function. Selection is entirely exogenous, i.e. humans pressing KEEP or REVERT. Which means the sim's selection layer is a *policy* we get to pick, and the same core can run drift-only, truncation-selection, mate-choice, whatever regime we want to study.

The three things that usually make a pop-gen simulator gnarly (variable genome length, linkage, continuous trait
spaces) are all just... absent. By construction. That's why 10k individuals is comfortable on a laptop instead of
"reserve a cluster."

The **k0** key is what makes the exchange verifiable between real badges. In simulation it drops out entirely. We
don't need to encrypt gametes to move them between in-memory individuals. The wire format itself - two phases,
AES-256-GCM-SIV under `k0`, base45 QRs, byte-15 badge-type semantics - is specified in
[`nastea1/dc34-gamete/PROTOCOL.md`](https://github.com/nastea1/dc34-gamete/blob/main/PROTOCOL.md); the vim gene
tab's QR panel implements that wire format so a simulated diploid can mint gametes real badges would accept.

---

## 1. The Genome

A **Haploid** is nine `u8` loci, in this order (matches `dc34-api/src/lib.rs:231`):

| # | locus         | role                                       |
|---|---------------|--------------------------------------------|
| 0 | `cd_period`   | number of animation periods (0..=6 valid)  |
| 1 | `cd_rate`     | animation speed                            |
| 2 | `cd_dir`      | animation direction bias                   |
| 3 | `sat`         | saturation                                 |
| 4 | `hue_ratedir` | hue-cycling speed and direction            |
| 5 | `hue_base`    | low end of hue interval                    |
| 6 | `hue_bound`   | high end of hue interval                   |
| 7 | `chaser`      | chaser animation (renderer calls it `lin`) |
| 8 | `nonlin`      | nonlinear brightness / gamma rolloff       |

A **Diploid** is an ordered pair `(haplo0, haplo1)` of Haploids. Order matters, and it matters in a way that's
going to bite us in §4, so remember that.

In vectorized form this is naturally an `(N, 2, 9)` `uint8` array for a population of `N` individuals: axis 0 =
individual, axis 1 = haploid slot, axis 2 = locus. (The shipped browser sim, `index.html`,
flattens this to a single `Uint8Array` of length `N * 2 * 9`; a NumPy reference implementation is in §10.)

```
population.shape == (N, 2, 9)
population.dtype == uint8
```

Every operation described in this document is a pure array transform over that shape. No per-individual loops
required, in either the shipped JS or the NumPy reference in §10.

---

## 2. Badge Types Are Priors, Not Species

The eight badge types are NOT phenotypic categories. This is the thing to internalize. 
They are **per-locus prior distributions** used only at fresh-genome generation. 

Once a genome exists, the badge type is metadata that tags along; it affects the inbreeding penalty (§5) 
but does not constrain future allele values. A Human's offspring can end up with an Uber-native chaser allele,
and there is nothing in the code that says no.

The ranges, transcribed from `dc34-api/src/lib.rs:102-173`:

| type       | hue         | sat        | chaser     | nonlin    | cd_dir    | cd_period_max |
|------------|-------------|------------|------------|-----------|-----------|---------------|
| Goon       | 0..=20 \*   | 160..=255  | 90..=255   | 0..=255   | 0..=255   | 4             |
| Community  | 32..=80     | 32..=160   | 90..=255   | 0..=255   | 0..=255   | 2             |
| Village    | 80..=128    | 32..=160   | 90..=255   | 0..=255   | 0..=45    | 4             |
| Human      | 128..=160   | 32..=255   | 90..=255   | 0..=255   | 0..=255   | 5             |
| Other      | 160..=192   | 16..=255   | 0..=255    | 0..=90    | 0..=255   | 6             |
| CtfContest | 192..=220   | 16..=255   | 90..=255   | 0..=90    | 0..=255   | 6             |
| Uber       | 220..=255 † | 130..=255  | 0..=45     | 0..=44    | 0..=45    | 3             |
| None       | 128..=160   | 32..=255   | 90..=255   | 0..=255   | 0..=255   | 4             |

Little gotchas from the source that are easy to miss:

- \* Goon's `hue_base` is **forced to 0** regardless of the range roll. No democracy there.
- † Uber's `hue_bound` is **forced to 255** regardless of the range roll. Uber gets the whole rainbow whether they want it or not.
- `hue_bound` is drawn from `hue_base..=range.end()`, so it's always ≥ `hue_base` at generation. Mutation can later violate that ordering; `phenotype()` fixes it at expression time.
- `cd_period` at generation is drawn from `0..=cd_period_max`. Mutation is post-processed with `% 7`, so it stays in `0..=6` afterward too.
- `cd_rate` and `hue_ratedir` are unconstrained `u8` at generation for every badge type. The wild west.

The 21-31 hue gap falls out immediately: no type's `hue_range()` covers that interval, so no fresh badge starts
there. This is the "unused hue slice" shown on the chart. Crossbreeding and mutation can still reach it, and
watching the gap fill in is one of the more satisfying things to plot.

For simulation, a badge type reduces to a small dict-of-ranges. The population starts with a chosen mixture, e.g.
70% Human, 10% Goon, 5% each of Village/Community/Other/CtfContest, 1% Uber, mirroring whatever conference
demographics you care about.

---

## 3. Meiosis: Independent Assortment, But With Linkage Groups

From `Diploid::meiosis()` in `dc34-api/src/lib.rs:357-374`. The linkage structure here is subtle and it's the one
place a naive re-implementation will probably be silently wrong.

- `cd_period`, `cd_rate`, `cd_dir` all inherit from the **same** parent haploid (one shared coin flip).
- `sat` gets its **own** coin flip.
- `hue_ratedir`, `hue_base`, `hue_bound` all inherit from the **same** parent (one shared coin flip).
- `chaser` gets its **own** coin flip.
- `nonlin` gets its **own** coin flip.

So the nine loci form **five linkage groups**:

```
group A: (cd_period, cd_rate, cd_dir)
group B: (sat)
group C: (hue_ratedir, hue_base, hue_bound)
group D: (chaser)
group E: (nonlin)
```

Within a group, the two loci co-segregate perfectly (linkage = 1). Across groups, segregation is independent
(linkage = 0). There's no partial crossover, no recombination hotspots, none of that biology-textbook stuff. Five
coins, five groups.

Vectorized: for a population of N mating pairs, sample five independent `uniform{0,1}` masks of length N, one per
group, and use each mask to select whichever haploid slot contributes that group's alleles.

```javascript
// pop: flat Uint8Array of length N*2*N_LOCI (individual, haploid slot, locus).
// out: flat Uint8Array of length N*N_LOCI, one gamete per individual.
const GROUPS = [[0,1,2], [3], [4,5,6], [7], [8]];  // five linkage groups

for (let i = 0; i < N; i++) {
  for (let g = 0; g < GROUPS.length; g++) {
    const pick = rng.coin();                       // 0 or 1, one coin per group
    const haploOff = i * 2 * N_LOCI + pick * N_LOCI;
    for (const locus of GROUPS[g]) {
      out[i * N_LOCI + locus] = pop[haploOff + locus];
    }
  }
}
```

There is no true "chromosome" object in the badge - the diploid is just two flat 9-byte haploids. The linkage structure is expressed purely by which
loci share a coin flip in `meiosis()`. Any simulator that draws one coin per locus (the naive "each locus
segregates independently" version) will overstate diversity in groups A and C. Easy mistake to make, hard to notice
from the output.

---

## 4. Phenotype: The Diploid --> Haploid Expression Map

From `Diploid::phenotype()` in `dc34-api/src/lib.rs:317-347`. Let `a = haplo0[locus]`, `b = haplo1[locus]`, all
arithmetic on `u8` with the saturating-add semantics noted:

| locus         | expression                                            | flavor                      |
|---------------|-------------------------------------------------------|-----------------------------|
| `cd_period`   | `min(6, (a + b) / 2)`                                 | mean, capped                |
| `cd_rate`     | `(a + b) / 2` in u16                                  | true mean, no wrap          |
| `cd_dir`      | `sat_add(a, b)`                                       | additive dominance          |
| `sat`         | `sat_add(a, b)`                                       | additive dominance          |
| `hue_ratedir` | `(2 + (14 - min(14, sat_add(a, b)))) % 14`            | inverse-add mod 14          |
| `hue_base`    | `min(a, b)`                                           | wider interval dominant     |
| `hue_bound`   | `max(a, b)` then `max(hue_bound, hue_base)`           | wider interval dominant     |
| `chaser`      | `sat_add(a, b)`                                       | additive dominance          |
| `nonlin`      | **`sat_add(haplo0.chaser, haplo1.nonlin)`**           | **asymmetric, see §4.1**    |

Where `sat_add(x, y) = min(255, x + y)`.

The structural stuff worth noting:

1. **Additive dominance is the default.** Six of the nine phenotypes get shoved toward 255 by any nonzero allele on
   either haploid. Under random mating with uniform priors, `sat`, `cd_dir`, `chaser` equilibrate very close to
   saturated within a handful of generations, because `P(sat_add(a,b) < 255) = P(a + b < 256)` and both `a` and `b`
   drift upward as saturation accumulates in the gene pool.

2. **`hue_base` = min, `hue_bound` = max.** Under drift, the phenotype interval `[min, max]` widens monotonically
   because min-drift shrinks and max-drift grows. Mutation can push them back but the phenotype itself has a strong
   "wider is dominant" bias. So populations tend to display broader hue rainbows over time even without anyone
   selecting for it.

3. **`hue_ratedir` is peculiar.** The inner clamp caps `a + b` at 14, then subtracts from 14, then adds 2, then
   mods 14. The output space is `{0, 1, 2, ..., 13}` and the mapping is not injective, many `(a, b)` pairs produce
   the same phenotype. Empirical mapping in the simulator takes about six lines.

4. **The `nonlin` line is asymmetric.** Which brings us to...

### 4.1 The nonlin Anomaly: Shouldn't Work, Ships Anyway

The line at `dc34-api/src/lib.rs:342`:

```rust
nonlin: self.0[0].chaser.saturating_add(self.0[1].nonlin),
```

The recipient's chaser allele shows up in the donor's nonlin phenotype. For simulation, the important facts are:

- **Haploid order matters.** Slot 0 is the recipient's egg, slot 1 is the donor's sperm (`replace_gene(egg,
  sperm)` builds `Diploid([egg, sperm])` in that order). So the *displayed* `nonlin` depends on which badge
  received the exchange, not just on the two contributing genomes. Same two badges, opposite roles, different
  nonlin.
- **The recipient's chaser allele leaks into the donor's nonlin phenotype.** This is not a modeling error on my
  part - it is the shipped code. Bug? Feature? You decide!
- Empirically, roughly 78% of parent pairs produce a different `phenotype().nonlin` if you reverse the haploid
  order. That's not a rounding-error asymmetry, that's structural.

At the population level `nonlin` and `chaser` are **phenotypically coupled** in a way meiosis does not model, and
the coupling is directional (egg-side chaser only). Any equilibrium analysis of `nonlin` that assumes symmetric
additive combination is going to be off by a lot.

For a "corrected model" run, replace slot-0 `chaser` with slot-1 `nonlin` (or slot-0 `nonlin`; both are defensible
symmetric readings). Running both variants side by side with identical RNG seeds and diffing the equilibria is one
of the cleanest ways to quantify the bug's actual population-level effect. The sim can make answering that question a breeze.

---

## 5. Mutation

From `mutate()`, `mutation_func()`, and the Gray-code helpers in `dc34-api/src/lib.rs:537-584`, plus `MutationRate` in `:178-225`.

**Invariant:** every mating exchange has up to three mutation events: (1) the sender ALWAYS mutates its sperm before sending, at the sender's own `final_rate`; (2) the receiver ALWAYS mutates its egg, at plain `rate` for cross-type or elevated `rate` for same-type; (3) if same-type (inbreeding), the sperm gets an EXTRA mutation pass on the receiver side at the elevated rate. The sender-side sperm mutation is the one that operates on the QR data anyone else can observe.

### 5.1 How Often A Locus Mutates

Each of the nine loci gets its own independent Bernoulli roll:

```
p_mutate(rate) = rate_value / 256
```

Where `rate_value` is the numeric byte of the rate:

| rate         | value | per-locus P(mutate) |
|--------------|-------|---------------------|
| None         | 0     | 0                   |
| Baseline     | 64    | 0.250               |
| Elevated     | 100   | 0.391               |
| Radioactive  | 140   | 0.547               |
| Apocalyptic  | 240   | 0.938               |

Rolls are per-locus, so the number of loci mutated in one pass is `Binomial(9, p_mutate(rate))`. At Apocalyptic
you're mutating nearly every locus every pass.

### 5.2 How Big A Mutation Is (or: why Gray codes are ingenious here)

If a locus mutates, the mutation is a **Gray-code bit flip**:

```
new_allele = gray_decode(gray_encode(old) ^ (bits << shift))
```

Where:

- `bits` is the rate-dependent mask: `Baseline=0x01, Elevated=0x03, Radioactive=0x07, Apocalyptic=0x1F`.
- `shift ~ Uniform{0..=7}`, resampled per locus per pass.

Semantics: encode the allele into its Gray code, XOR in a run of `popcount(bits)` contiguous bits at a random 8-bit
position (with wrap-around at bit 8 discarded because `shift ≤ 7`, so the top bit of the flip mask can fall off),
decode back.

**Why Gray code matters here.** In standard binary, flipping bit *k* of a `u8`
changes the value by `±2^k` regardless of what the current value is. So a bit-3 flip is always ±8, boring, uniform.
In Gray code, a bit flip in Gray-space is a *monotone-neighborhood* move in ordinary space: adjacent Gray codes
differ by exactly one bit, so a one-bit Gray flip is a small ordinary-space step (often, though not always, `±1` or
a nearby jump). This means Baseline (`bits=0x01`, one-bit flip in Gray) produces mostly *small* value changes,
whereas Apocalyptic (`bits=0x1F`, five contiguous bits flipped in Gray) produces genuinely large, structurally
disruptive changes. The rate knob isn't just "how often" but "how far."

For simulation, precompute the Gray transforms as 256-entry LUTs and vectorize the mutation as a bit-XOR on the
whole population gamete array. Sample one shift per locus per individual as `uniform{0..7}`, form the mask as
`(bits << shift)` truncated to `u8` (matching the Rust `u8` overflow behavior), and XOR into the Gray-encoded
allele.

Post-mutation, one locus needs special handling: `cd_period` is `mutation_func(gene, bits) % 7`, keeping it in
`0..=6`. Every other locus keeps the raw mutated `u8`.

### 5.3 The Mutation Meter Is UX, Not Math

`mutation_param` (`dc34-vault/src/main.rs:302`) is a physical-input integrator: side-rocker reversals add +2,
accelerometer orientation changes add +8, redraws decay it (faster above 160).
`MutationRate::from_param(mutation_param)` bins it into the four levels.

In simulation this whole subsystem collapses to: **the mutation rate is a scalar per exchange, chosen by policy.**
There is no wall clock, no physical input, no decay. You can model whatever tempo you want explicitly as a policy
on the mating scheduler. Or, you can just pass in whatever rate you want to study.

This is also where you'd bolt on any policy the actual physical system can't express: type-conditional
attraction (Goons preferentially seek Speakers), refractory periods after mating, love triangles, whatever.
The mating scheduler is the natural home for that kind of thing, see §8.

### 5.4 lock_rate() and final_rate

Two rates matter in the shipped code (`dc34-vault/src/config.rs:306, 342, 357`):

- **`mutation_rate`** - the live value derived from `mutation_param`, changing continuously.
- **`final_rate`** - the value captured at the moment the QR is displayed (`lock_rate()`).

Whenever a badge produces a gamete for a QR (`get_padded_gamete()` at `config.rs:341`), the sender **always** runs `mutate(gamete, final_rate)` before handing it over. This sender-side pass is unconditional — it happens on every gamete production regardless of who the recipient is, whether it's same-type, or whether inbreeding applies. The `final_rate` here is the sender's own locked rate at QR display time (see below), NOT anything the recipient chose.

The snapshot semantic exists so the visible level at the moment of exchange is what governs the mutation, not whatever the meter has decayed to by the time the QR actually gets scanned. Scanning takes real time, and this prevents people from cheesing it by locking high and letting the meter decay before the exchange.

In simulation there is no meter to decay, so `final_rate` = whatever rate the policy chose for this mating event.
It only matters if you want to model the human "lock high, exchange later" trick.

---

## 6. The Inbreeding Pass

From `dc34-vault/src/main.rs:704-732`.

On top of the unconditional sender-side pass in §5.4, when the incoming badge type equals the receiver's own badge type, the sperm gets a **second, extra** mutation pass on the receiving side, at a rate that is `max(inbreeding_floor, donor_final_rate)`:

- Human's inbreeding floor is `Elevated` (100/256).
- Every other type's inbreeding floor is `Baseline` (64/256).
- `None` counts as Baseline in this table too.

The receiver's egg is then also mutated at that same rate (see `get_egg(rate)` at `config.rs:353-362`: it uses `rate.unwrap_or(self.final_rate).max(self.final_rate)`).

So a same-type Human/Human exchange at Baseline effectively looks like:

```
donor sperm       : mutate(rate=Baseline)   [on the donor side]
recipient inbreed : mutate(rate=Elevated)   [on the recipient side, extra pass on sperm]
recipient egg     : mutate(rate=Elevated)   [egg gets the inbreeding rate too]
```

Whereas a cross-type Human/Goon Baseline exchange is just:

```
donor sperm       : mutate(rate=Baseline)
recipient egg     : mutate(rate=Baseline)
```

**No inbreeding pass on the sperm, and no elevation on the egg.** Cross-type is cheap; same-type is expensive.

In pop-gen language this is an **assortative-mating mutation penalty**: mating within a subpopulation costs extra
mutational load. The design intent (per the source comment) is to force isolated populations to diverge faster. In
simulation it's a one-flag toggle, and the effect on within-type variance is one of the cleanest experiments to run.

---

## 7. The Full Mating Transform

Putting §3, §5, §6, §4 together, a single mating event `(egg_parent, sperm_parent, rate) → child` is:

```
1. sperm = meiosis(sperm_parent)          # §3
   mutate(sperm, rate)                    # §5 -- ALWAYS, on the sender side, every gamete
2. egg = meiosis(egg_parent)              # §3
3. if sperm_parent.type == egg_parent.type:
       inbreed_rate = max(rate,
                          Elevated if egg_parent.type == Human else Baseline)
       mutate(sperm, inbreed_rate)        # §6, extra pass on the sperm
       mutate(egg,   inbreed_rate)        # §6, egg gets the elevated rate
   else:
       mutate(egg, rate)                  # egg uses the plain rate
4. child.diploid = Diploid(egg, sperm)    # slot 0 = egg, slot 1 = sperm
5. child.phenotype = phenotype(child.diploid)   # §4
```

The `sperm` gets **two** mutation passes on same-type mating (once during gamete formation at the donor, once on receipt at the recipient's inbreeding pass). Not a bug, it's the shipped behavior, and it explains why same-type Apocalyptic exchanges are so much more disruptive than cross-type ones. Every step here is a whole-population array op given a mating schedule.

---

## 8. Selection Is Whatever You Say It Is

The badge implements **no fitness function.** The only decision is KEEP vs REVERT (`dc34-vault/src/config.rs`,
`prior_gene` slot), and it's made by a human squinting at the previewed phenotype. So in the sim, selection is a
*policy* you plug in. Some useful ones:

- **Neutral drift.** KEEP always. Every child replaces the parent. The null hypothesis, always run this first.
- **Truncation on a trait.** E.g. KEEP iff `phenotype.chaser < 88` (chasing the shooting-star variant), or iff
  `hue_base ∈ [21, 31]` (colonizing the hue gap), or iff `sat == 255`.
- **Weighted selection.** KEEP with probability proportional to some scalar function of phenotype (e.g. "prefer
  wider hue intervals").
- **Human-in-the-loop analog.** KEEP with probability that depends on Euclidean distance in phenotype space from a
  target, plus a temperature. Models a human with imperfect but directional taste.

The generality is the whole point. Because selection is not baked into the core, the same simulator answers "how
fast does drift alone fix an Uber chaser allele in a Human population" and "how fast does directional selection at
rate p fix it" with only the policy swapped. Free experiment.

Similarly, the **mating scheduler** is a policy:

- **Panmictic random.** Every generation, shuffle and pair everyone. Easy, unrealistic.
- **Type-biased.** Weight pairing probability by badge-type similarity (models people naturally clustering by badge
  color).
- **Geographic.** Assign 2D positions, pair by proximity, walk them.
- **badge.sex remote.** A small fraction of pairings are drawn from the global population regardless of position.  The tinder for badges.

At scale (10k individuals or more), these policies all start to produce *different equilibria*.

---

## 9. What To Measure

Some quantities that fall naturally out of the model:

1. **Allele frequency spectrum per locus, per generation.** For each locus, the empirical histogram over `0..=255`
   in the population's `2N` haploid slots. Under neutral drift it drifts; under selection it skews.

2. **Effective diversity, per locus.** `H = 1 - Σ p_i^2` over allele frequencies. Rate of loss ≈ `1/(2N)` per
   generation under neutral drift, faster under bottlenecks. Classic pop-gen number, cheap to compute.

3. **Fixation and loss rates.** Track how often a locus goes to a single allele across the population, and how
   quickly. This is where the inbreeding pass should really show up.

4. **Cross-type introgression.** Fraction of Human-badge haploids carrying an allele in the Uber-native chaser
   range `0..=45`. Answers "how long does it take Uber chaser genes to spread into the Human population under given
   mixing assumptions."

5. **P(phenotype.chaser < 88) per generation.** The shooting-star variant. Compare to the pure-source calculation
   (~99.7% Uber, ~6% Other, 0% for the other six at gen 0) and watch it evolve.

6. **21-31 hue-gap occupancy.** Fraction of population with `phenotype.hue_base ∈ [21, 31]` OR interval covering
   that range. Answers "does the gap get colonized under mixing / mutation, and how fast."

7. **`nonlin` bug delta.** Run the population twice with identical RNG seeds, once with the shipped asymmetric
   expression, once with a symmetric version. Compare the population-level `nonlin` distributions. This
   *quantifies* the "very likely a bug" hypothesis in a way source reading alone cannot.

8. **Same-type vs cross-type mating variance.** With the inbreeding pass on and off. Should demonstrate the design
   claim ("adds more diversity more quickly for populations that are isolated").

None of these require anything the badge itself doesn't already do. They only require running it 10,000 times in
parallel and looking at the aggregate.

---

## 10. Reference Implementation Sketch

The shipped simulator is JavaScript, in [`index.html`](https://github.com/charlesreid1/dc34sim/blob/gh-pages/index.html), using a single flat `Uint8Array` for the
population. The sketch below is a NumPy reference of the same math - a second implementation is useful for
cross-checking, and NumPy makes the array-transform shape of each step easier to read than the flattened JS
version. Every function is a pure array transform.

```python
import numpy as np

N_LOCI = 9
LOCUS = dict(cd_period=0, cd_rate=1, cd_dir=2, sat=3,
             hue_ratedir=4, hue_base=5, hue_bound=6,
             chaser=7, nonlin=8)

GROUPS = [(0, 1, 2), (3,), (4, 5, 6), (7,), (8,)]

RATE_VALUES = dict(None_=0, Baseline=64, Elevated=100,
                   Radioactive=140, Apocalyptic=240)
RATE_BITS   = dict(None_=0, Baseline=0x01, Elevated=0x03,
                   Radioactive=0x07, Apocalyptic=0x1F)

# Gray-code lookups, precomputed once.
_g = np.arange(256, dtype=np.uint8)
GRAY_ENCODE = (_g ^ (_g >> 1)).astype(np.uint8)
GRAY_DECODE = np.empty(256, dtype=np.uint8)
GRAY_DECODE[GRAY_ENCODE] = _g

def meiosis(pop, rng):
    N = pop.shape[0]
    picks = rng.integers(0, 2, size=(N, len(GROUPS)), dtype=np.uint8)
    out = np.empty((N, N_LOCI), dtype=np.uint8)
    idx = np.arange(N)
    for g, loci in enumerate(GROUPS):
        for locus in loci:
            out[:, locus] = pop[idx, picks[:, g], locus]
    return out

def mutate(gametes, rate_name, rng):
    val  = RATE_VALUES[rate_name]
    bits = RATE_BITS[rate_name]
    if val == 0:
        return gametes
    N = gametes.shape[0]
    # per-locus roll, per-locus shift
    rolls   = rng.integers(0, 256, size=(N, N_LOCI), dtype=np.uint8) < val
    shifts  = rng.integers(0, 8,   size=(N, N_LOCI), dtype=np.uint8)
    masks   = ((bits << shifts) & 0xFF).astype(np.uint8)
    encoded = GRAY_ENCODE[gametes]
    flipped = encoded ^ (masks * rolls.astype(np.uint8))
    mutated = GRAY_DECODE[flipped]
    # cd_period wraps mod 7 post-mutation
    mutated[:, LOCUS['cd_period']] %= 7
    return np.where(rolls, mutated, gametes)

def phenotype(pop):
    a, b = pop[:, 0, :].astype(np.int16), pop[:, 1, :].astype(np.int16)
    out = np.empty(pop.shape[0::2], dtype=np.uint8)  # (N, 9)
    out[:, LOCUS['cd_period']]   = np.minimum(6, (a[:, 0] + b[:, 0]) // 2)
    out[:, LOCUS['cd_rate']]     = ((a[:, 1] + b[:, 1]) // 2).astype(np.uint8)
    out[:, LOCUS['cd_dir']]      = np.minimum(255, a[:, 2] + b[:, 2])
    out[:, LOCUS['sat']]         = np.minimum(255, a[:, 3] + b[:, 3])
    hr = np.minimum(14, a[:, 4] + b[:, 4])
    out[:, LOCUS['hue_ratedir']] = ((2 + (14 - hr)) % 14).astype(np.uint8)
    hb = np.minimum(a[:, 5], b[:, 5])
    ub = np.maximum(a[:, 6], b[:, 6])
    out[:, LOCUS['hue_base']]    = hb.astype(np.uint8)
    out[:, LOCUS['hue_bound']]   = np.maximum(ub, hb).astype(np.uint8)
    out[:, LOCUS['chaser']]      = np.minimum(255, a[:, 7] + b[:, 7])
    # asymmetric nonlin: slot0 chaser + slot1 nonlin. This is the bug.
    out[:, LOCUS['nonlin']]      = np.minimum(255, a[:, 7] + b[:, 8])
    return out
```

`mate()` composes `meiosis`, `mutate`, and the inbreeding branch against a pair of parent arrays and their
badge-type arrays. `step()` runs one generation given a mating scheduler and a selection policy.

At 10,000 individuals, one generation is a few milliseconds. Running 1000 generations across 20 seeded replicates
is a coffee break.

---

## 11. But Why Are We Doing This

Others have performed a source audit of the badge. That's useful - it tells us what the code does, and it can
hand-calculate the fresh-generation `lin < 88` probabilities per badge type. What it can't tell us is what a
*population* of these things does over time. And the physical conference can't tell us either, at any reasonable
statistical power: the population is too small, the mating is too sparse, and the timeline is just one weekend
where only a fraction of attendees will be doing the badge-mating dance.

A Monte Carlo re-implementation lets us:

- **Quantify hypotheses.** "The nonlin bug matters" becomes "the nonlin bug shifts the population median displayed nonlin from X to Y at generation 100 under panmictic Baseline mating."
- **Falsify strategies.** Anyone can use this tool to propose a hypothesis, build a strategy, and test it. The simulator can tell us the *expected number of generations* to fixation under various mutation rates, mating rates, and population fractions. The result may confirm the hypothesis, or it may show that the target variant is out of reach on any realistic conference timeline.
- **Explore counterfactuals.** What if the inbreeding pass didn't exist? What if Uber was 5% of the population instead of <1%? What if `hue_range()` had no gap? Each of these is a one-line change in the simulator.

None of this replaces what the badges are actually doing at DEF CON. It runs *alongside* source audits and field
observations, answering the kind of questions you can only answer at scale, with 10,000+ badges.

---

## 12. Source Provenance

The mathematical core in this document is transcribed from:

- `dc34-api/src/lib.rs` - `Haploid`, `Diploid`, `BadgeType`, `MutationRate`, `Haploid::from_type`, `Diploid::meiosis`, `Diploid::phenotype`, `mutate`, `mutation_func`, `gray_encode`, `gray_decode`.
- `dc34-vault/src/config.rs` - `final_rate`, `lock_rate`, `get_padded_gamete`, `get_egg`, `replace_gene`.
- `dc34-vault/src/main.rs` - `mutation_param` accumulation and decay; the inbreeding-detection and rate-elevation branch around line 704.

All ranges, expressions, and probabilities are lifted directly from those files. Anything labeled *hypothesis* or *policy* in this document (the selection regimes, the mating schedulers) is deliberately not in the badge source. Those are extensions the simulator is free to make because the physical badge imposes no such constraints on the mathematical core.
